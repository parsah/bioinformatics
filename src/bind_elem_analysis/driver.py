'''
A script which takes a series of inputs:
* A user-provided fasta file of input sequences
* A reference genome
Given the set of inputs above, each is mapped against the reference genome 
(using BLAST), and the location is then derived. Upsteam and downstream of this
hit-site is then sliced out. These two set of slices are then compared against
a dictionary of binding-elements (BEs) to derive which BEs exist and which do
not.
'''

import argparse


if __name__ == '__main__':
	str_epilog = 'Ex: python driver.py Fasta_file Ref_genome Db -u 2000 -d 200'
	str_desc = 'Drives the ability to find binding-elements (BEs) in a set of'+\
		' sequences. Given a an populated input fasta file, each sequence is'+\
		' mapped against a reference such that up-stream and down-stream regions'+\
		' are analyzed for the presence of BEs. Such BEs are from a user-provided'+\
		' dictionary.'
	p = argparse.ArgumentParser(description=str_desc, epilog=str_epilog)
	p.add_argument('input_fasta', help='Input fasta file; fasta [required]')
	p.add_argument('ref', help='Reference genome; fasta [required]')
	p.add_argument('db', help='Binding Element DB [required]')
	p.add_argument('-u', '-up', help='up-stream bp', default=2000, type=int)
	p.add_argument('-d', '-down', help='down-stream bp', default=200, type=int)
	
	args = p.parse_args()
	print args.input_fasta, args.ref, args.db
	print args.u
	