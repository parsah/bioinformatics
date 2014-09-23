"""
Given a user-provided FASTA file, this script segregates all FASTA
entries by their sequence length. Resultant groupings are subsequently
saved in their own filename. This functionality is ideal for instances
whereby sequences within a file must all be by the same length.
"""

import argparse
from Bio import SeqIO


def partition(f):
    """
    Segregates FASTA entries by their sequence length and saves each
    respective segregation to its own concordant file.
    :param f: FASTA filename.
    """
    records = list(SeqIO.parse(open(f), 'fasta'))
    unique_lens = set([len(record.seq) for record in records if len(record.seq) % 100 == 0])
    groupings = {unique_len: [] for unique_len in unique_lens}  # key => length, value => FASTA records.

    # save length-specific sequences to their own data-structure.
    for record in records:
        record_len = len(record.seq)
        if record_len in groupings:
            groupings[record_len].append(record)

    # next, save length-specific sequences to their own respective file.
    for length in groupings:
        out_handle = open(f + '.len.' + str(length), 'w')
        SeqIO.write(groupings[length], out_handle, 'fasta')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', metavar='FASTA', required=True,
                        help='User-provided FASTA file [req]')
    args = vars(parser.parse_args())
    partition(f=args['f'])
