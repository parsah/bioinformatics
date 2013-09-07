
import argparse

def main():
	desc = 'Filters a fasta file given a list of required fasta entries'
	p = argparse.ArgumentParser(description=desc)
	p.add_argument('-in', help='Input fasta file [na]', 
				metavar='FILE',required=True)
	p.add_argument('-list', help='List of fasta headers to keep [na]', 
				metavar='FILE', required=True)
	args = vars(p.parse_args())
	seqs = parse_fasta(fname=args['in']) # parse master fasta file
	seqs_keeping = parse_list(fname=args['list']) # parse file of desired entries
	get_desired_seqs(total_seqs=seqs, to_keep=seqs_keeping) # run analysis

# Run analysis; sends sequences to standard-out
def get_desired_seqs(total_seqs, to_keep):
	for seq in to_keep:
		if seq not in total_seqs: # if a desired seq is not in the master, exit
			print(seq, 'is not in the fasta file')
			break
		else: # if the desired seq is found, print-out to stdout stream
			print('>'+seq+'\n'+total_seqs[seq])

# Parse a list of user-provided fasta headers; entries to keep
def parse_list(fname):
	seq_to_keep = []
	for line in open(fname):
		line = line.strip()
		seq_to_keep.append(line)
	return seq_to_keep
	
# Trivial function to parse a fasta file
def parse_fasta(fname):
	try:
		seqs, header = {}, ''
		for line in open(fname):
			line = line.strip()
			if line.startswith('>'):
				header = line[1:] # remove the fasta header
				seqs[header] = ''
			else:
				seqs[header]+=line
		return seqs
	except IOError:
		print('[ERROR] input file:', fname, 'not found')

if __name__ == '__main__':
	main()