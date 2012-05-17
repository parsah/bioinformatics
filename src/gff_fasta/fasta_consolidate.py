'''
A simple function to concat a multi-line fasta line into a single line
'''

import optparse
import os

# function to parse multi-line fasta and yield only one single line each
def parse_file(in_file, out_file):
	results = {}
	header = ''
	for line in open(in_file):
		line = str(line.strip())
		if line.startswith('>'):
			header = line
			results[header] = ''
		else:
			results[header] += line
	if len(results) > 0:
		# save the consolidated results
		out_handle = open(out_file, 'w')
		for i in results:
			out_handle.write(i+'\t'+ str(results[i])+'\n')
			out_handle.flush()
		out_handle.close()
		print len(results), ' entries saved [OK]'
	else:
		print 'Either file is not fasta or it does not contain any entries'

def main():
	# create argument to parse user input file
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", help="Multi-line fasta file")
	parser.add_option("-o", "--output", help="Consolidated output in_file")
	(options, args) = parser.parse_args()
	
	# test to see if in_file is valid
	in_file, out_file = options.input, options.output
	if os.path.isfile(in_file):
		parse_file(in_file, out_file)
	

if __name__ == '__main__':
	main()