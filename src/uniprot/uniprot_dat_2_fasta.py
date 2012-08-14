
# Simple script to convert a Swissprot / Uniprot file to a fasta file.
# Ambigious annotations ('uncharacterized', 'hypothetical', etc) are not parsed

from Bio import SeqIO
import argparse

# instantiate parser; only swissprot .dat file is required
parser = argparse.ArgumentParser(description='Simple uniprot (.dat) to fasta parser')
parser.add_argument('-in', help='Input .dat file', metavar='', required=True)
args = vars(parser.parse_args())

# parse input file
entries = SeqIO.parse(open(args['in']), format='swiss')
for entry in entries: # iterate over each entry and determine if annotation is ambigious
	if 'uncharacterized' in entry.description\
		or 'Uncharacterized' in entry.description\
		or 'Hypothetical' in entry.description\
		or 'hypothetical' in entry.description\
		or 'unknown' in entry.description\
		or 'Unknown' in entry.description:
		pass
	else: # if annotation is specific, send to standard-out
		print '>'+entry.id+'|'+entry.name+'\n'+entry.seq
