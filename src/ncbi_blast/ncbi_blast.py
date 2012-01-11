
""" 
Simple NCBI blast implementation.
"""

import time
from optparse import OptionParser
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# BLAST-SPECIFIC VARIABLES
PROGRAM = 'blastp'
DATABASE = 'nr'
EXPECT = 1e-5
USE_HIT_NUMBER = 0

# KEYWORDS FOR FILTERING HITS
AMBIGIOUS_KEYWORDS = ['hypothetical', 'putative', 'unknown', 'unnamed', 'predicted']
NUM_ENTRIES = 0

def parse_results(records):
	# iterate through each hit and print its information
	results = []
	for each_record in records:
		alignments = each_record.alignments
		if alignments: # if there are hits
			for each_alignment in alignments:
				hsps = each_alignment.hsps
				for each_hsp in hsps:
					title = each_alignment.title # the title has delimiters
					title = str(title).split('|') # which can be separated.
					# slide accession number and desc information 
					accession, description = title[3].strip(), title[4].strip()
					eval = float(each_hsp.expect)
					score = int(each_hsp.score)
					if eval <= EXPECT:
						results.append([accession, description, score, eval])
	
	out_str = get_best_hit(results)
	return out_str

def get_best_hit(results):
	new_list = [] # dedicate this array to store hits without ambigious chars
	for each_blast_hit in results:
		accession, desc, score, eval = each_blast_hit
		num_times_not_found = 0
		for each_kword in AMBIGIOUS_KEYWORDS:
			if each_kword not in desc.lower():
				num_times_not_found+=1
		
		if num_times_not_found == len(AMBIGIOUS_KEYWORDS):
			new_list.append([accession, desc, score, eval])
#			print "\tPASSED - no ambig", accession, desc, score, eval, num_times_not_found
	
	new_list.sort(key=lambda x: x[-2], reverse=True)
	
	out_str = ''
	if new_list:
		best_hit = new_list[0]
		accession, desc, score, eval = best_hit
		out_str = accession+'\t'+desc+'\t'+str(score)+'\t'+str(eval)
	else:
		out_str = '-\t-\t-\t-'

	
	return out_str

def get_num_entries():
	global NUM_ENTRIES
	handle = open(IN_FILENAME)
	for each_line in handle:
		each_line = each_line.strip()
		if '>' in each_line:
			NUM_ENTRIES+=1

def begin_blast(input_filename):
	get_num_entries()
	print '\nSTARTING BLAST - ' + time.asctime()+'\tID:'+str(UNIQUE_ID)+'\n'
	in_handle = open(input_filename, 'r') # read fasta file and parse
	in_handle = SeqIO.parse(in_handle, 'fasta') 
	
	output_handle = open(OUT_FILENAME, 'w')
	blast_xml_filename = 'blast_results_'+str(UNIQUE_ID)+'.xml'
	for counter, each_fasta in enumerate(in_handle): # for each fasta entry ...
		header = str(each_fasta.description) # to-string its header
		seq = str(each_fasta.seq) # to-string the sequence and BLAST it
		print '['+str(counter+1)+'/'+str(NUM_ENTRIES)+'] - '+header+'... ',
		blast_result = NCBIWWW.qblast(PROGRAM, DATABASE, seq, EXPECT, descriptions=5)
		blast_result = blast_result.read() # read the object contents to string
		
		# save the output string to a local file for parsing
		save_file = open(blast_xml_filename, 'w')
		save_file.write(blast_result) # save the blast results to file (XML)
		save_file.close()
		
		# next, parse the xml results and display on-screen
		results_handle = open(blast_xml_filename, 'r')
		records = NCBIXML.parse(results_handle)
		result_string = parse_results(records)
		output_handle.write(header+'\t'+result_string)
		output_handle.flush()
		output_handle.write('\n')
		output_handle.flush()
		print 'OK'
	output_handle.close()
	print '\nBLAST COMPLETE - ' + time.asctime()
	
if __name__ == '__main__':
	global IN_FILENAME, OUT_FILENAME, UNIQUE_ID
	parser = OptionParser()
	parser.add_option("-i", "--in_file", help="input fasta file")
	parser.add_option("-o", "--out_file", help="output txt file")
	parser.add_option("-v", "--value", help="a unique value")
	
	(options, args) = parser.parse_args()
	
	if not options.in_file or not options.out_file or not options.value:
		print 'values must be provided for each argument. see -h for info'
		
	else:
		IN_FILENAME = str(options.in_file)
		OUT_FILENAME = str(options.out_file)
		UNIQUE_ID = str(options.value)
		
		if '.txt' in OUT_FILENAME:	
			OUT_FILENAME = OUT_FILENAME.replace('.txt', '_'+UNIQUE_ID+'.txt')
			print 'input:', IN_FILENAME
			print 'output:', OUT_FILENAME
			print 'unique ID:', UNIQUE_ID
			begin_blast(IN_FILENAME)
		