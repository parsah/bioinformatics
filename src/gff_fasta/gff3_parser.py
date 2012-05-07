
import optparse
import os


# simple function to determine if filename is valid
def is_valid_file(filename):
	return os.path.isfile(filename)

# parse the GFF3 filename into a simple hash
def parse_gff(gff_filename):
	print 'parsing gff: ', gff_filename

# parse the fasta file which will be mapped against the GFF3 annotation
def parse_fasta(fasta_filename):
	hash_fasta = {}
	header = '' # fasta header
	for line in open(fasta_filename):
		line = line.strip() # remove its trailing whitespace
		if '>' in line:
			header = line[1:] # add each fasta entry to the hash
			hash_fasta[header] = ''
		else:
			hash_fasta[header]+= line
	print 'fasta parsing [OK], #/entries: ', len(hash_fasta)
	return hash_fasta

if __name__ == '__main__':
	parser = optparse.OptionParser()
	# pass-in a GFF3 as well as a fasta file
	parser.add_option("-g", "--gff3", help="GFFv3 input file")
	parser.add_option("-f", "--fasta", help="input fasta file")
	(options, args) = parser.parse_args()
	
	if (is_valid_file(options.gff3) and is_valid_file(options.fasta)):
		parse_fasta(options.fasta)
		parse_gff(options.gff3)
	else:
		print ''
	