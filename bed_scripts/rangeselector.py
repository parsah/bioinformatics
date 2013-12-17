'''
Given a BED file and a set of genomic FASTA files, each genomic region in the 
BED file is paired with an equal-length genomic segment, matching its GC and 
repeat content. Selecting such sequences therefore helps extract sequences 
representative of a BED file which match its intrinsic properties; serving as
a baseline or control-set.
'''

import argparse
import random
import os
from Bio import SeqIO # for reading FASTA files.

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

class BEDSequence():
    def __init__(self, seq):
        self.seq = seq
        
    def gc_perc(self):
        gc_count = float(self.seq.count('G') + self.seq.count('C'))
        return round(gc_count / self.get_len() * 100, 2)
    
    def get_seq(self):
        return self.seq
    
    def get_len(self):
        return len(self.seq)

class RangeSelector():
    ''' 
    Given genomic ranges defined in a user-provided BED file, this class
    selects corresponding genomic sequences having equivalent GC and
    repetitive properties as the genomic range. These resultant sequences
    are therefore equivalent in intrinsic properties and suitable for 
    contrast-based analyses.
    '''
    def __init__(self, args):
        self.args = args
        self.fasta_folder = args['folder']
        self.bed = self.args['bed']
        self.outhandle = open(self.args['o'], 'w')        
    
    def run_selector(self):
        ''' 
        Executes the main algorithm that takes a BED entry and extracts a 
        corresponding genomic sequence which matches its GC and repeat content.
        '''
        
        # create set of valid FASTA files before any analysis is even performed.
        fasta_files = [f for f in os.listdir(self.fasta_folder)
                       if f.endswith('.fasta')]
        
        for bed_line in open(self.bed):
            bed_line = bed_line.strip().split('\t')
            # parse the BED entry.
            chrom, start, end = bed_line[0], int(bed_line[1]), int(bed_line[2])
            print('Analyzing BED entry:',chrom, start, end)
            
            # pull-out the sequence referencing the the BED chromosome.
            record = SeqIO.read(self.fasta_folder + '/' + chrom+'.fasta', 'fasta')
            bed_seq = BEDSequence(seq=record.seq[start: end])
            
            rand_chrom_name = random.choice(fasta_files)
            print('=> BED GC %:', bed_seq.gc_perc(), '; querying', rand_chrom_name, '...')
            rand_chrom_rec = SeqIO.read(self.fasta_folder + '/' + rand_chrom_name, 'fasta')
            
            # perform GC-counting operations.
            gcc = GCCounter(seq = str(rand_chrom_rec.seq), win = bed_seq.get_len(), 
                            basegc = bed_seq.gc_perc(), delta = self.args['d'], 
                            top_num = self.args['t'])
            hits = gcc.count_gc()
            
            self.outhandle.write('>matched.'+chrom+'.'+str(start)+'.'+str(end) + '\n' +random.choice(hits) + '\n')
            self.outhandle.flush()
            
            print()
        self.outhandle.close()

if __name__ == '__main__':
    desc = 'Script to match BED entries with feature-specific genomic regions.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-bed', metavar='FILE', help='BED file [req]',
                        required=True)
    parser.add_argument('-folder', metavar='DIR', help='Hard-masked folder of FASTA files [req]',
                        required=True)
    parser.add_argument('-d', metavar='FLOAT', help='Difference between BED GC percentage and FASTA [0.1]',
                        default=0.1, type=float)
    parser.add_argument('-t', metavar='INT', help='Fetch the top N matching sequences [10]',
                        default=10, type=int)
    parser.add_argument('-o', metavar='FILE', help='FASTA output file [./sequences.fasta]',
                        default='./sequences.fasta')
    args = vars(parser.parse_args())
    try:
        rs = RangeSelector(args)
        rs.run_selector()
    except IOError as e:
        print(e)
    except KeyboardInterrupt:
        print()
