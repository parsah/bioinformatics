
from Bio.Seq import Seq
import re, argparse

def parse_fasta(infile):
	seqs = {}
	header = ''
	for line in open(infile):
		line = line.strip()
		if '>' in line:
			
			header = line.replace('|', '_')
			seqs[header] = ''
		else:
			seqs[header]+=line
	print('#/fasta entries:', len(seqs),'\n')
	return seqs

def get_ORFs(sequence, title, cutoff):
	sequence = Seq(sequence).upper()
	dict_frame_and_protein = {}
	
	# Get Open Reading Frames on the forward strand (+)
	for frame in range(3):
		seq_in_frame = sequence[frame:]
		protein_seq = seq_in_frame.translate(to_stop=False)
		dict_frame_and_protein['frame_'+str(frame+1)] = str(protein_seq)

	# Get Open Reading Frames on the r strand (-)
	sequence = sequence.reverse_complement()
	for frame in range(3):
		seq_in_frame = sequence[frame:]
		protein_seq = seq_in_frame.translate(to_stop=False)
		dict_frame_and_protein['frame_-'+str(frame+1)] = str(protein_seq)
	return _fetch_orfs_with_met(dict_frame_and_protein, title, cutoff)
	
def _fetch_orfs_with_met(dict_frame_and_protein, title, cutoff):
	# for all the 6x open-reading frames, you do not want to simply blast each
	# reading frame because there are many MET start-sites within this large
	# sequence. What you want to do is split this large reading frame sequence
	# into sub-sequences whereby each sub-sequence starts with a MET start-site.
	
	max_orf_len, longest_orf, info = -1, '', '' # represents most likely ORF
	for frame_number in dict_frame_and_protein:
		frame_protein_seq, all_orfs = dict_frame_and_protein[frame_number], {}
		matcher = re.findall('M[A-Z]*\*', frame_protein_seq) # find M sites
		if len(matcher) >= 1:
			for counter, each_orf in enumerate(matcher):
				all_orfs[frame_number+'_suborf_'+str(counter+1)] = each_orf
		
		for orf in all_orfs: # for all ORFs, find the longest reading frame
			if len(all_orfs[orf]) >= max_orf_len:
				max_orf_len = len(all_orfs[orf])
				longest_orf = all_orfs[orf]
				info = title+'_'+orf # set information regarding it
	if len(longest_orf) >= cutoff: # the longest ORF must pass cutoff-length
		return info+'\n'+longest_orf
	else:
		return ''
	
def _get_protein_length(protein_seq):
	raw_protein_length = len(protein_seq)
	num_stop_codon = str(protein_seq).count('*')
	return raw_protein_length - num_stop_codon

def main():
	desc = 'Helpful script to find the longest ORF from amongst all six ORFs'
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument('-in', help='Input fasta file', metavar='FILE')
	parser.add_argument('-cutoff', help='ORF length cutoff [10]', 
					default=10, type=int, metavar='INT')
	args = vars(parser.parse_args())
	if args['in']:
		entries = parse_fasta(args['in'])
		for entry in entries:
			print(get_ORFs(sequence=entries[entry], title=entry, cutoff=args['cutoff']))
		

if __name__ == '__main__':
	main()
