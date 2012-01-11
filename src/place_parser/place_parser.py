
""" 
PLACE parser - a simple for parsing PLACE (plant cis-element DB) ASCII files. 
"""

def parse_place_file(filename):
	curr_block = [] # represents each sub-section in this ascii file
	num_blocks = 0
	for line in open(filename):
		line = line.strip()
		if '//' not in line: # create a chunk until delimiter is reached
			curr_block.append(line)
		else:
			num_blocks+=1
			process_block(num_blocks, curr_block)
			curr_block = []
	print '#/blocks parsed:', num_blocks

def process_block(block_num, block):
	accn, desc_block, organism = '', '', ''
	for each_line in block:
		header, desc = each_line[0:3].strip(), each_line[3:].strip()
		if header == 'AC':
			accn = desc
		if header == 'DE':
			desc_block = desc_block + ' ' + desc
		if header == 'OS':
			organism = organism + ' '+desc
	sequence = block[-1].strip()
	result = [accn, desc_block.strip(), organism.strip(), dna_to_regex(sequence)]
	print '\t'.join(result)

def dna_to_regex(sequence):
	# for easy regex parsing, convert ambigious nucleotides into their
	# respective regular-expression equivalent.
	sequence = sequence.replace('R', '[AG]') # R == A or G
	sequence = sequence.replace('Y', '[TC]') # Y == T or C; pyrimidines
	sequence = sequence.replace('K', '[GT]') # K == G or T; keto
	sequence = sequence.replace('M', '[AC]') # M == A or C; amino
	sequence = sequence.replace('S', '[GC]') # S == G or C; 3x H+ bonds
	sequence = sequence.replace('W', '[AT]') # W == A or T; 2x H+ bonds
	sequence = sequence.replace('B', '[CGT]') # B == not A
	sequence = sequence.replace('D', '[AGT]') # D == not C
	sequence = sequence.replace('H', '[ACT]') # H == not G
	sequence = sequence.replace('V', '[ACG]') # V == not T
	sequence = sequence.replace('N', '[ATGCU]') # N == any character

	valid_chars = set(['A', 'T', 'G', 'C', 'U', '[', ']'])
	for character in sequence:
		if character not in valid_chars:
			raise RuntimeError('error with seq:', sequence)
	return sequence
		
if __name__ == '__main__':
	parse_place_file('/home/bioinfo/Desktop/place.seq')