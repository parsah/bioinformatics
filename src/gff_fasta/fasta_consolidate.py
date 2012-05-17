'''
A simple function to concat a multi-line fasta line into a single line
'''

import optparse
import os

# function to parse multi-line fasta and yield only one single line each
def parse_file(filename):
	pass

def main():
	# create argument to parse user input file
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", help="Multi-line fasta file")
	(options, args) = parser.parse_args()
	
	# test to see if filename is valid
	filename = options.input
	if os.path.isfile(filename):
		parse_file(filename)
	

if __name__ == '__main__':
	main()