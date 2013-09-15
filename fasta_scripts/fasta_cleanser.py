'''
Script which removed FASTA entries with user-provided invalid characters.
'''

import argparse
from Bio import SeqIO

def run(fname, chars):
    entries = SeqIO.parse(fname, 'fasta')
    for i in entries:
        is_clean = True
        for c in chars:
            if c in i.seq:
                is_clean = False
        if is_clean in i.description:
            print('>' + i.description + '\n' + str(i.seq))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, metavar='FASTA',
                        help='FASTA file [na]')
    parser.add_argument('-chars', required=True, metavar='', nargs='+',
                        help='List of characters to filter [na]')
    args = vars(parser.parse_args())
    run(fname = args['i'], chars = args['chars'])
    