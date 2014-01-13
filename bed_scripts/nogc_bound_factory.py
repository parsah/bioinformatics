'''
A script which pulls-out sequences from a genome that match a given BED entry
length. This script is a variant of the GC-bound script that only pulls
corresponding genomic sequences within a user-defined GC percentage.
'''

import random
import os
import sys
from Bio import SeqIO # for reading FASTA files.

BED_FILE = sys.argv[1] # first argument
GENOME_DIR = sys.argv[2] # folder containing FASTA sequences

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
            fasta_entries[record.description] = {'seq': record.seq, 'len': len(record.seq)}
            
    for bed_line in open(BED_FILE):
        bed_line = bed_line.strip().split('\t')
        # parse the BED entry.
        chrom, start, end = bed_line[0], int(bed_line[1]), int(bed_line[2])
        #print('\nAnalyzing BED entry:', chrom, start, end)
        bed_len = end - start
        
        # choose a random chromosome and extract a suitable segment from it.
        keys = list(fasta_entries.keys())
        keys.remove(chrom) # remove the current ID to avoid duplicates.
        rand_chrom_name = random.choice(keys)
        rand_chrom_record = fasta_entries[rand_chrom_name]
        pos = random.randrange(0, rand_chrom_record['len']) # find random position
        seq = rand_chrom_record['seq'][pos: pos + bed_len] # extract seq
        
        # write results; only one hit is randomly chosen.
        print('>matched|'+chrom+':'+str(start)+'-'+str(end) + '\n' + seq)

if __name__ == '__main__':
    try:
        run_selector()
    except IOError as e:
        print(e)
    except KeyboardInterrupt:
        print()