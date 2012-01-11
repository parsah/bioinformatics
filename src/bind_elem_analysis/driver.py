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
import tempfile
import time

# represents all the files and states pertaining to the analysis
class Analysis_Bucket():
	def __init__(self, fasta_file, ref, be_db, bp_up, bp_down):
		self.fasta_file = fasta_file
		self.ref_genome = ref
		self.be_db = be_db
		self.bp_up = bp_up
		self.bp_down = bp_down
		self.output_temp_dir = tempfile.mkdtemp(dir=os.getcwd())
	
	@staticmethod
	def new_segment_line():
		return '*'*40
		
	def __str__(self):
		return 'Binding Element Search Script - ['+str(time.asctime())+']\n'+\
			'Fasta file: '+str(self.fasta_file)+'\n'+\
			'Reference genome: '+str(self.ref_genome)+'\n'+\
			'Binding Element (BE) DB: '+str(self.be_db)+'\n'+\
			'up-stream bp: '+str(self.bp_up)+'\n'+\
			'down-stream bp: '+str(self.bp_down)+'\n'+\
			Analysis_Bucket.new_segment_line()
	
	
	def is_param_valid(self):
		print '### Checking validity of input files: ###'
		is_valid = os.path.exists(self.fasta_file)
		# check to see if fasta file is found
		if is_valid is False:
			print '[ERROR: Fasta file (',self.fasta_file,') is not found]'
		
		# check to see if reference genome is found
		is_valid = os.path.exists(self.ref_genome)
		if is_valid is False:
			print '[ERROR: Reference genome (',self.ref_genome,') is not found]'
		
		# check to see if binding-element database is found
		is_valid = os.path.exists(self.be_db)
		if is_valid is False:
			print '[ERROR: Binding Element DB (',self.be_db,') is not found]'
		
		# stdout status of validity
		if is_valid:
			print 'VALID: Input parameters all valid. Proceeding with analysis.'	
		else:
			print 'INVALID: Input parameters invalid. Please refer to message above.'	
		return is_valid

	def parse_fasta(self, fasta_file):
		print Analysis_Bucket.new_segment_line()
		print '### Parsing fasta file and beginning analysis ###'
		entries = {}
		for line in open(self.fasta_file):
			line = line.strip()
			if len(line) == 0:
				break
			else:
				if '>' in line:
					header = line
					entries[header] = ''
				else:
					entries[header]+=line
		
		print '#/fasta entries parsed:',len(entries)
		os.rmdir(self.output_temp_dir)

	# for each fasta entry, blast against reference genome
	def run_analysis(self):
		# first, parse the fasta file
		self.parse_fasta(self.fasta_file)

if __name__ == '__main__':
	str_epilog = 'Ex: python driver.py Fasta_file Ref_genome Be_db -u 2000 -d 200'
	str_desc = 'Drives the ability to find binding-elements (BEs) in a set of'+\
		' sequences. Given a an populated input fasta file, each sequence is'+\
		' mapped against a reference such that up-stream and down-stream regions'+\
		' are analyzed for the presence of BEs. Such BEs are from a user-provided'+\
		' dictionary.'
	p = argparse.ArgumentParser(description=str_desc, epilog=str_epilog)
	p.add_argument('-i','-input', help='Input fasta file; fasta [required]',required=True)
	p.add_argument('-r','-ref', help='Reference genome BLAST DB; fasta [required]', required=True)
	p.add_argument('-b','-be', help='Binding Element DB [required]', required=True)
	p.add_argument('-u', '-up', help='up-stream bp <def: 2000>', default=2000, type=int)
	p.add_argument('-d', '-down', help='down-stream bp <def: 200>', default=200, type=int)
	
	args = p.parse_args()
	analysis_obj = Analysis_Bucket(args.i, args.r, args.b, args.u, args.d)
	print analysis_obj
	if analysis_obj.is_param_valid(): # if input is valid, process or stop
		analysis_obj.run_analysis()
	
	
	