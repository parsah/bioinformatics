''' 
A very useful script to filter FASTA sequences if a specified filter-string
is provided. Sequences satisfying this criteria will be kept.
'''

import sys
import argparse
from Bio import SeqIO

def filter_seqs(fname, pattern):
    # Read FASTA files and see if pattern is found in their headers.
    fasta_records = SeqIO.parse(fname, 'fasta')
    for record in fasta_records:
        if pattern in record.description: # print if pattern in the header.
            sys.stdout.write(">"+record.description+"\n"+\
                             str(record.seq)+"\n")
            sys.stdout.flush()
    sys.stdout.close()
    

# A simple script to remove FASTA entries containing a specific substring.
parser = argparse.ArgumentParser()
parser.add_argument('-in', metavar='FASTA', required=True,
                    help='Input FASTA file [na]')
parser.add_argument('-pattern', metavar='STR', required=True,
                    help='Substring for which to keep FASTA entries [na]')
args = vars(parser.parse_args())
filter_seqs(fname=args['in'], pattern=args['pattern'])
