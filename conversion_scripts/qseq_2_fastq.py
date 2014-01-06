'''
Useful script to convert an Illumina QSEQ file into FASTA
'''

import argparse

# Drives the parsing of the qseq file and prints stdout
def parse(input_file):
	for line in open(input_file):
		line = line.strip().split()
		header = '@'+line[0]+':'+line[2]+':'+line[3]+':'+line[4]+':'+line[5]+'#'+line[6]
		seq = line[8]
		comment = '+'
		qual_str = line[9]
		print(header+'\n'+seq+'\n' + comment + '\n' + qual_str)

# Simple script to convert a qseq file to fastq; conversion sent to stdout	
def main():
	try:
		description = 'qseq -> fastq: easy and quick conversion script'
		parser = argparse.ArgumentParser(description)
		parser.add_argument('-input', help='Input QSEQ file', 
						metavar='', required=True)
		
		# parse args and run script
		args = vars(parser.parse_args())
		if not args['input']:
			parser.print_help()
		else:
			parse(args['input'])
	except IOError as e:
		print(e, 'is an invalid file')
	
if __name__ == '__main__':
	main()