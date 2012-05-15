
""" 
A simple multiprocess NCBI blast implementation; similar in logic to
ncbi_blast.py
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import multiprocessing
from multiprocessing import Pool
import time, os

# BLAST VARIABLES
DB = 'nr'
PROG = 'blastp'
E_VAL = 1e-5

# NUMBER OF PROCESSES
NUM_PROCESSES = multiprocessing.cpu_count()

# PROG-SPECIFIC VARIABLES
AMBIGIOUS_KEYWORDS = ['hypothetical', 'putative', 'unknown', 'unnamed', 'predicted']
ALL_HITS = []

NUM_COMPLETED = 0
NUM_RECORDS = 0 # assigned while file is read

# input fasta file to iteratively run blast on
IN_FILE = 'seqs.fasta'
# each successful blast will have its XML stored in a directory. This directory
# will then be parsed to fetch corresponding hits/results.
OUT_XML_DIR = 'blast_dir/'

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
					if eval <= E_VAL:
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

def blast_seq(fasta):
	try:
		output = NCBIWWW.qblast(program=PROG, database=DB, sequence=fasta.seq, expect=E_VAL, descriptions=2)
		xml_output = output.read() # xml output is an xml file
		
		xml_filename = OUT_XML_DIR+fasta.id+'.xml'
		xml_handle = open(xml_filename, 'w')
		xml_handle.write(xml_output)
		xml_handle.flush()
		xml_handle.close()
			
		print fasta.id
		return (fasta.id, xml_output) # fasta header followed by xml string
	except Exception:
		print '\nERROR:', fasta.id,' - ', fasta.seq, '- BLASTING AGAIN\n'
		blast_seq(fasta)
		
def get_number_fasta_records():
	global NUM_RECORDS
	handle = open(IN_FILE, 'r')
	for each_line in handle:
		each_line = each_line.strip()
		if '>' in each_line:
			NUM_RECORDS+=1

def init_blast():
	global ALL_HITS
	get_number_fasta_records()
	handle = open(IN_FILE, 'r')
	fasta_handle = SeqIO.parse(handle, 'fasta') # read-in fasta file
	
	records = [record for record in fasta_handle]
	
	pool = Pool(processes=NUM_PROCESSES) # create a pool for workers
	print '#/processes created:',NUM_PROCESSES, '-',time.asctime()
	ALL_HITS = pool.map(blast_seq, records)

def process_hits():
	global ALL_HITS, OUT_XML_DIR
	out_handle = open(OUT_XML_DIR, 'w')
	for each_result in ALL_HITS:
		header, xml_string = each_result # last index is the xml string
		
		xml_filename = 'blast_output.xml'
		xml_handle = open(xml_filename, 'w')
		xml_handle.write(xml_string)
		xml_handle.flush()
		xml_handle.close()
		
		records = NCBIXML.parse(open(xml_filename))
		out_string = header+'\t'+parse_results(records)
		out_handle.write(out_string+'\n')
		out_handle.flush()
	out_handle.close()

		
if __name__ == '__main__':
	try:
		os.mkdir(OUT_XML_DIR) # create the output folder
		init_blast() # use multiple processes to run blast
#		process_hits() # results are then fetched and saved to the txt file
		print 'analysis complete', time.asctime()
		
	except KeyboardInterrupt:
		print 'BLAST interrupted by user'