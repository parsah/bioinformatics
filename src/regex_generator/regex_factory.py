
import argparse, collections

''' 
Given a variable number of equal-length raw DNA strings, this module generates
a regular expression to model the provided input sequences.
'''

# Checks to validate whether input sequences are equal in length
def __check_sequences(list_seq):
	if len(list_seq) < 2:
		print '[ERROR] At least 2 list_seq must be provided'
	else:
		# if all strings have same length, a list-set will be of size 1
		all_length = len(set([len(seq) for seq in list_seq]))
		if all_length == 1:
			return True # all sequences are the same length
		else:
			print '[ERROR] All sequences must be of the same length'
			return False

# Create mappings of what char is found at what sequence-index
def generate_regex(sequences):
	dict_positions = collections.defaultdict(set)
	for each_seq in sequences:
		for index, each_char in enumerate(each_seq):
			dict_positions[index].add(each_char)
			
	# Once mapping are performed, represent maps as a regex
	regex_str = ''
	for index in dict_positions:
		index_combinations = ''.join(dict_positions[index])
		if len(dict_positions[index]) > 1:
			index_combinations = '[' + index_combinations + ']'
		regex_str += index_combinations # concat index's combinations into one
	
	print regex_str	

def main():
	desc = 'Regular Expression generator given equally-sized strings'
	usage='%(prog)s sequences [seq1, seq2, ...]'
	parser = argparse.ArgumentParser(usage, description=desc)
	parser.add_argument('sequences', type=str, nargs='+',
					help='List of equally-sized raw DNA sequences')
	args = vars(parser.parse_args())
	
	# test if > 1 sequence is provided and all sequences are equal in length
	if __check_sequences(args['sequences']):
		generate_regex(args['sequences'])

if __name__ == '__main__':
	main()