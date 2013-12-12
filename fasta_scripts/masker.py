'''
Script which removes FASTA entries with user-provided invalid characters.
'''

import argparse
import re
from Bio import SeqIO

def display_seq(header, seq):
    ''' 
    Prints sequence header and sequence information in FASTA format
    @param header: FASTA header
    @param seq: FASTA sequence
    '''
    print('>' + header + '\n' + seq)

def run(fname, is_filter):
    entries = SeqIO.parse(fname, 'fasta')
    for i in entries:
        seq = str(i.seq).upper()
        num_a, num_t, num_n = seq.count('A'), seq.count('T'), seq.count('N')
        num_g, num_c = seq.count('G'), seq.count('C')
        if (num_a + num_c + num_g + num_t + num_n) == len(seq):
            display_seq(header = i.description, seq = seq) # clean sequence
        elif is_filter:
            corrected_seq = re.sub('[^ATGCN]', 'N', seq) # filter bases
            display_seq(header = i.description, seq = corrected_seq)

if __name__ == '__main__':
    d = """
    FASTA filter script. Non-ATGCN bases are by-default filtered. The --mask
    option replaces such bases with an N so that the sequence is kept.
    """
    parser = argparse.ArgumentParser(description=d)
    parser.add_argument('-i', required=True, metavar='FASTA',
                        help='FASTA file [na]')
    parser.add_argument('--mask', action='store_const', const=True,
                        default=False, help='Mask non-ATGCN bases as N [false]')
    args = vars(parser.parse_args())
    run(fname = args['i'], is_filter = args['mask'])
    
