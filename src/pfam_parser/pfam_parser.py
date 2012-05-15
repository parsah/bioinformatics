
''' 
Given a tab-delimited text file with PFAM accessions, this script parses this 
file, extracts the PFAM accession and retrieves this respective accessions'
PFAm functional annotation.
'''

from optparse import OptionParser
import os.path

def parse_blast(options):
	in_file = options.in_file # parse-in the user-provided input file
	if (os.path.isfile(in_file)):
		handle = open(in_file) # read-in the user-provided file
		for line in handle:
			line = line.strip().split('\t')
			glyma, pfam_block = str(line[0]), str(line[2])
			print glyma+'\t'+pfam_block
		
	else:
		print 'Please provide a valid input file'

def main():
	parser = OptionParser()
	parser.add_option("-i", "--in_file", help="Input BLAST file; tab-delimited")
	
	(options, args) = parser.parse_args()
	parse_blast(options)
	
if __name__ == '__main__':
	main()

