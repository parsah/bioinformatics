'''
A simple function to concat a multi-line fasta line into a single line
'''

import optparse
import os

N = 2000 # default chunk size; #/sequences per file

# function to parse multi-line fasta and yield only one single line each
def parse_file(in_file, out_dir, show_stats):
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
		print len(results), ' entries saved [OK]'
		if show_stats:
			for entry in results:
				print entry[1:]+'\t', len(results[entry])
		else:
			partition(results, out_dir)
	else:
		print 'Either file is not fasta or it does not contain any entries'

def partition(results, out_dir):
	# obtain keys of the dict for easy slicing
	l = results.keys()
	chunks = [l[i:i+N] for i in range(0, len(l), N)]
	
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
	# iterate over each fragment and save it to its respective file
	for frag_num, fragment in enumerate(chunks):
		
		# create output file to represent each fragment
		out_handle = open(out_dir+'/fragment_' +str(frag_num)+'.fasta', 'w')
		
		# iterate over each sequence mapping in that fragment, save its seq
		for a_key in fragment:
			out_handle.write(a_key+'\n' + results[a_key]+'\n')
			out_handle.flush()
		out_handle.close()	
	
	print len(chunks), 'fragments'


def main():
	# create argument to parse user input file
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", help="Multi-line fasta file")
	parser.add_option("-o", "--output", help="Output directory")
	parser.add_option('--stat', action='store_true', default=False,
				help='Only stdout length of fasta entry [false]')
	(options, args) = parser.parse_args()
	
	# test to see if in_file is valid
	in_file, out_dir, show_stat = options.input, options.output, options.stat
	if os.path.isfile(in_file):
		parse_file(in_file, out_dir, show_stat)
	

if __name__ == '__main__':
	main()