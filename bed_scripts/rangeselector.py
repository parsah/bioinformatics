'''
Given a BED file and a set of genomic FASTA files, each genomic region in the 
BED file is paired with an equal-length genomic segment, matching its GC and 
repeat content. Selecting such sequences therefore helps extract sequences 
representative of a BED file which match its intrinsic properties; serving as
a baseline or control-set.
'''

import random
import os
import sys
from Bio import SeqIO # for reading FASTA files.
from Bio.SeqUtils import GC

BED_FILE = sys.argv[1] # first argument
GENOME_DIR = sys.argv[2] # folder containing FASTA sequences

class GCCounter():
    ''' 
    A useful class which efficiently computes the GC percentage of a
    user-provided input sequence as well as its hard-masked repeat
    percentage. The user provides a maximally-accepted GC percentage (basegc)
    and a cushion at-which to accept GC percentages (delta). All sliding
    windows with GC percentages within this delta are therefore accepted.
    '''
    def __init__(self, seq, win, basegc, delta, top_num):
        self.seq = seq # input query sequence
        self.win = win # sliding window size
        self.top_num = top_num # save the number of GC-matching sequences.
        self.basegc = basegc # the baseline GC percentage
        self.delta = delta # the buffer or cushion at-which a window is accepted.
        self.upper_bound = self.basegc + self.delta # upper, lower bounds of GC %s.
        self.lower_bound = self.basegc - self.delta
    
    def count_gc(self):
        prior_gc = 0 # define initial GC count
        prior_seq = '' # define empty prior sequence
        hits = []
        for i in range(len(self.seq)):
            if i == 0: # for the first sliding window, compute GC content
                prior_seq = self.seq[i: self.win]
                prior_gc = prior_seq.count('G') + prior_seq.count('C')
                gc_perc = round((prior_gc / self.win) * 100, 2)
                if self.lower_bound <= gc_perc <= self.upper_bound:
                    #print(prior_seq+'\t'+str(gc_perc))
                    hits.append(prior_seq)
            else:
                curr_seq = self.seq[i: i+self.win]
                if len(curr_seq) != self.win: # do not pull pre-mature windows
                    break
                
                # if the prior sequence's head is G or C, decrease count
                # since this base is not part of current window.
                if prior_seq[0] == 'G' or prior_seq[0] == 'C':
                    prior_gc -= 1
                
                # if the current sequence's tail (new addition) contains
                # a G or C, increase count since this base is within window.
                if curr_seq[-1] == 'G' or curr_seq[-1] == 'C':
                    prior_gc += 1
                
                # update reference
                prior_seq = curr_seq
                gc_perc = round((prior_gc / self.win) * 100, 2)
                if self.lower_bound <= gc_perc <= self.upper_bound:
                    #print(curr_seq+'\t'+str(gc_perc))
                    hits.append(curr_seq)
            
            # if max number of sequences reached, break
            if len(hits) == self.top_num:
                break
        return hits

def run_selector():
    ''' 
    Executes the main algorithm that takes a BED entry and extracts a 
    corresponding genomic sequence which matches its GC and repeat content.
    '''
    
    # create set of valid FASTA files before any analysis is even performed.
    fasta_entries = {}
    for f in os.listdir(GENOME_DIR):
        if f.endswith('.fasta') and 'chrM' not in f:
            record = SeqIO.read(GENOME_DIR + '/' + f, 'fasta')
            fasta_entries[record.description] = record.seq
            #print('added', record.description, len(record.seq))
    
    for bed_line in open(BED_FILE):
        bed_line = bed_line.strip().split('\t')
        # parse the BED entry.
        chrom, start, end = bed_line[0], int(bed_line[1]), int(bed_line[2])
        #print('\nAnalyzing BED entry:', chrom, start, end)
        
        # pull-out the sequence referencing the the BED chromosome.
        record = fasta_entries[chrom]
        bed_sequence = record[start: end] # slice segment

        # choose a random chromosome and extract a suitable segment from it.
        keys = list(fasta_entries.keys())
        keys.remove(chrom) # remove the current ID to avoid duplicates.
        rand_chrom_name = random.choice(keys)
        
        gc_perc = round(GC(bed_sequence), 2)
        #print('=> BED GC %:', gc_perc, '; querying', rand_chrom_name, '...')
        #print(keys)
        rand_chrom_record = fasta_entries[rand_chrom_name]
        
        # perform GC-counting operations.
        counter = GCCounter(seq = str(rand_chrom_record), win = len(bed_sequence), 
                        basegc = gc_perc, delta = 0.2, 
                        top_num = 10)
        hits = counter.count_gc() # get all the matching segments for the BED entry.
        
        # write results; only one hit is randomly chosen.
        print('>matched.'+chrom+'.'+str(start)+'.'+str(end) + '\n' +random.choice(hits) + '\n')

if __name__ == '__main__':
    try:
        run_selector()
    except IOError as e:
        print(e)
    except KeyboardInterrupt:
        print()
