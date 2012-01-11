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
import os
import time


def print_params(fasta_file, ref_genome, be_db, bp_up, bp_down):
	print 'Binding Element Search Script - [',time.asctime(),']' 
	print 'Fasta file:', fasta_file
	print 'Reference genome:', ref_genome
	print 'Binding Element (BE) DB:', be_db
	print 'up-stream bp:', bp_up
	print 'down-stream bp:', bp_down
	print '*'*40
	return is_param_valid(fasta_file, ref_genome, be_db, bp_up, bp_down)
	
	
def is_param_valid(fasta_file, ref_genome, be_db, bp_up, bp_down):
	print 'Checking input parameters:'
	is_valid = os.path.exists(fasta_file)
	# check to see if fasta file is found
	if is_valid is False:
		print '[ERROR: Fasta file (',fasta_file,') is not found]'
	
	# check to see if reference genome is found
	is_valid = os.path.exists(ref_genome)
	if is_valid is False:
		print '[ERROR: Reference genome (',ref_genome,') is not found]'
	
	# check to see if binding-element database is found
	is_valid = os.path.exists(be_db)
	if is_valid is False:
		print '[ERROR: Binding Element DB (',be_db,') is not found]'
	
	# stdout status of validity
	if is_valid:
		print 'VALID: Input parameters all valid. Proceeding with analysis.'	
	else:
		print 'INVALID: Input parameters invalid. Please refer to message above.'	
	return is_valid

if __name__ == '__main__':
	str_epilog = 'Ex: python driver.py Fasta_file Ref_genome Be_db -u 2000 -d 200'
	str_desc = 'Drives the ability to find binding-elements (BEs) in a set of'+\
		' sequences. Given a an populated input fasta file, each sequence is'+\
		' mapped against a reference such that up-stream and down-stream regions'+\
		' are analyzed for the presence of BEs. Such BEs are from a user-provided'+\
		' dictionary.'
	p = argparse.ArgumentParser(description=str_desc, epilog=str_epilog)
	p.add_argument('-i','-input', help='Input fasta file; fasta [required]',required=True)
	p.add_argument('-r','-ref', help='Reference genome; fasta [required]', required=True)
	p.add_argument('-b','-be', help='Binding Element DB [required]', required=True)
	p.add_argument('-u', '-up', help='up-stream bp <def: 2000>', default=2000, type=int)
	p.add_argument('-d', '-down', help='down-stream bp <def: 200>', default=200, type=int)
	
	args = p.parse_args()
	fasta_file, ref, be_db, bp_up, bp_down = args.i, args.r, args.b, args.u, args.d
	input_valid = print_params(fasta_file, ref, be_db, bp_up, bp_down)
	if input_valid:
		pass
	
	
	