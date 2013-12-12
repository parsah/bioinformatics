''' 
Removes FASTA entries if their headers are found in a user-provided list of
undesired headers.
'''

import argparse
from parser import *

def main():
	desc = 'Filters a fasta file given a list of required fasta entries'
	p = argparse.ArgumentParser(description=desc)
	p.add_argument('-in', help='Input FASTA file [na]', 
				metavar='FILE',required=True)
	p.add_argument('-list', help='List of FASTA headers to keep [na]', 
				metavar='FILE', required=True)
	args = vars(p.parse_args())
	seqs = parse_fasta(fname=args['in']) # parse master fasta file
	seqs_keeping = parse_list(fname=args['list']) # parse file of desired entries
	get_desired_seqs(total_seqs=seqs, to_keep=seqs_keeping) # run analysis

def get_desired_seqs(total_seqs, to_keep):
	''' 
	Runs analysis and prints all sequences which are not found in the list to
	standard-out.
	'''
	for seq in to_keep:
		if seq not in total_seqs: # if a desired seq is not in the master, exit
			print(seq, 'is not in the FASTA file')
			break
		else: # if the desired seq is found, print-out to stdout stream
			print('>'+seq+'\n'+total_seqs[seq])

def parse_list(fname):
	''' 
	Simple function to parse a user-provided file of FASTA headers.
	@return: list of FASTA headers.
	'''
	seq_to_keep = []
	for line in open(fname):
		line = line.strip()
		seq_to_keep.append(line)
	return seq_to_keep

if __name__ == '__main__':
	main()