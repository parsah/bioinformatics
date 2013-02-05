
import argparse, os

# Trivial function to parse a fasta file
def parse_fasta(fname):
	if os.path.exists(fname):
		seqs, header = {}, ''
		for line in open(fname):
			line = line.strip()
			if line.startswith('>'):
				header = line[1:] # remove the fasta header
				seqs[header] = ''
			else:
				seqs[header]+=line
		print(len(seqs), 'entries parsed')
		return seqs
	else:
		raise IOError('Input file:' + fname + ' is an invalid file')

# If a sequence header contains a semi-colon, that segment becomes a new 
# FASTA entry. This feature is useful in instances in-which you have alternate
# isoforms but they share the same sequence. Doing so makes management and
# readability easier and each isoform is its own entry.
def split(seqs):
	for header in seqs:
		variants = header.split('|') # eg. >gene|isoform1;isoform2;isoform3
		base = variants[0] # the actual 'gene' is called base
		variants = variants[1].split(';') # everything beyond is a variant
		for variant in variants: # for each variant, create a new entry
			h = '>' + base + '|' + variant +'\n' + seqs[header]
			print(h)
		
if __name__ == '__main__':
	desc = 'Splits FASTA entries if its header contains a semi-colon'
	p = argparse.ArgumentParser(description=desc)
	p.add_argument('-in', help='Input FASTA file [na]', 
				metavar='FILE',required=True)
	try:
		args = vars(p.parse_args())	
		seqs = parse_fasta(fname=args['in']) # parse master fasta file
		split(seqs) # split FASTA headers on whether a semicolon is found
	except IOError as e:
		print(e)