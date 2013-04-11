'''
A simple hash-based script which only produces a set of unique FASTA entries.
'''

import argparse
import hashlib
import sys
from Bio import SeqIO

def process_seqs(fname):
    ''' 
    For each sequence, compute a hash and test if this hash exists in our
    dictionary. If so, that respective sequence is a duplicate and therefore
    not represented as output.
    '''
    hashes = set()
    records = SeqIO.parse(fname, 'fasta')
    for record in records:
        sha512 = hashlib.sha512()
        header, sequence = record.description, str(record.seq)
        as_bytes = sequence.encode(encoding='utf-8')
        sha512.update(as_bytes)
        hex_repr = sha512.hexdigest()
        if hex_repr not in hashes: # if sequence is unique, send to stdout.
            hashes.add(hex_repr)
            sys.stdout.write('>' + header + '\n' + sequence+'\n')
            sys.stdout.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', metavar='FILE', required=True,
                        help='Input FASTA file [na]')
    args = vars(parser.parse_args())
    process_seqs(fname=args['in'])