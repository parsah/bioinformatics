''' 
Trivial script to count GC percentage of sequences in a FASTA file.
'''

import argparse
import os
from parser import *

def get_gc_content(sequence):
    ''' 
    Trivial function to compute GC percentage of a DNA string.
    @param sequence: DNA sequence
    @return: float referencing the sequence's GC percentage.
    '''
    sequence = str(sequence).upper()
    num_a = sequence.count('A')
    num_t = sequence.count('T')
    num_g = sequence.count('G')
    num_c = sequence.count('C')
    total = num_a + num_t + num_c + num_g
    perc = float(num_g + num_c) / total * 100
    return perc # return sequence GC percentage

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', metavar='FASTA', required=True,
                        help='Input FASTA file')
    args = vars(parser.parse_args())
    seqs = parse_fasta(fname = args['in'])
    print(os.path.basename(args['in']))
    for s in seqs:
        gc_perc = get_gc_content(sequence=seqs[s])
        print(gc_perc)
    
