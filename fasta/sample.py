'''
Helpful script to sample a user-provided FASTA file using the uniform
distribution. Such sampling takes place with replacement.
'''

import argparse
import random
from Bio import SeqIO


def sample(f, n):
    '''
    Sample (with replacement) a user-provided FASTA file.
    @param f: User-provided FASTA file.
    @param n: Number of samples from the FASTA file to perform.
    '''
    entries = list(SeqIO.parse(f, 'fasta'))
    for seqnum in range(n):
        loc = round(random.uniform(0, len(entries) - 1))
        entry = entries[loc]  # get index of randomly-selected FASTA entry
        header = '>' + str(seqnum + 1) + '-' + entry.description  # header
        print(header + '\n' + str(entry.seq))  # print-out entire entry

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', metavar='FASTA', required=True,
                        help='User-provided FASTA file [req]')
    parser.add_argument('-n', metavar='INT', type=int, default=10,
                        help='Number of samples to identify [10]')
    args = vars(parser.parse_args())
    sample(f=args['f'], n=args['n'])
