'''
Given a BED file and a set of genomic FASTA files, each genomic region in the 
BED file is paired with an equal-length genomic segment, matching its GC and 
repeat content. Selecting such sequences therefore helps extract sequences 
representative of a BED file which match its intrinsic properties; serving as
a baseline or control-set.
'''

import argparse
import os
import random
from Bio import SeqIO # for reading FASTA files.

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
        self.fasta_folder = args['folder']
        self.fasta_files = set() # references FASTA files within the folder.
        self.bed = args['bed']
        self.outfile = args['o']
        
    def bed_check(self):
        ''' 
        Verify that each BED entry maps to a FASTA file in the user-provided
        FASTA folder.
        '''
        self.fasta_files = set([os.path.splitext(f)[0]
                            for f in os.listdir(self.fasta_folder) 
                            if f.endswith('fasta')])
        for i in open(self.bed): # iterate over BED, check region matches a chromosome.
            i = i.strip().split('\t')
            chromosome = i[0] # chromosome column is the first column
            if chromosome not in self.fasta_files:
                raise IOError(chromosome + " is not an entry in the BED file [error].")

    def run_selector(self):
        ''' 
        Executes the main algorithm that takes a BED entry and extracts a 
        corresponding genomic sequence which matches its GC and repeat content.
        '''
        for bed_line in open(self.bed):
            bed_line = bed_line.strip().split('\t')
            chrom, start, end = bed_line[0], int(bed_line[1]), int(bed_line[2])
            print(chrom, start, end)
            record = SeqIO.read(self.fasta_folder + '/' + chrom+'.fasta', 'fasta')
            bed_seq = BEDSequence(seq=record.seq[start: end])
            print(bed_seq.gc_perc(), bed_seq.get_seq())

if __name__ == '__main__':
    desc = 'Script to match BED entries with feature-specific genomic regions.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-bed', metavar='FILE', help='BED file [req]',
                        required=True)
    parser.add_argument('-folder', metavar='DIR', help='Soft or hard-masked folder of FASTA files [req]',
                        required=True)
    parser.add_argument('-o', metavar='FILE', help='FASTA output file [./sequences.fasta]',
                        default='./sequences.fasta')
    args = vars(parser.parse_args())
    try:
        rs = RangeSelector(args)
        rs.bed_check()
        rs.run_selector()
    except IOError as e:
        print(e)
