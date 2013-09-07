
import argparse, multiprocessing, time, os
from Bio.Blast import NCBIXML
import string

# words to filter from annotations 
AMBIGIOUS_KEYWORDS = ['hypothetical', 'putative', 'unknown', 'unnamed', 'predicted', 
					'uncharacterized', 'Predicted']

# points to ncbi bin folder
BLAST_EVALUE = ''

DB_FOLDER, DB_NAME = '', ''
INPUT_FILE, NUM_PROCESSES = '', 0
NUM_ENTRIES, CURR_COUNTER = 0,0
NUM_WITH_HIT, NUM_WITHOUT_HIT = 0, 0

def callback(return_from_blast):
	fasta_id, fasta_file, xml_file = return_from_blast
	global CURR_COUNTER
	CURR_COUNTER+=1
	output_str = process_blast_xml(fasta_id, xml_file)
	
	print('[',CURR_COUNTER, '/', NUM_ENTRIES,']', fasta_id,'\t', output_str)
	OUT_HANDLE.write(fasta_id+'\t'+output_str+'\n')
	OUT_HANDLE.flush()
	
	# at the end, remove the fasta and xml file
	os.remove(fasta_file)
	os.remove(xml_file)

def create_fasta_file(fasta_id, fasta_seq):
	filename = DB_FOLDER+'/'+fasta_id+'.fasta'
	out_handle = open(filename, 'w')
	out_handle.write('>'+fasta_id[1:]+'\n'+str(fasta_seq)+'\n')
	out_handle.flush()
	out_handle.close()
	return filename

def onlyascii(char):
	if char not in string.printable: 
		return 'X'
	else:
		return char

def sanitize_xml(xml_fname): 
	curr_str = open(xml_fname, 'r').read()
	new_str = ''
	for achar in curr_str:
		new_str += onlyascii(achar)
	
	handle = open(xml_fname,'w')
	handle.write(new_str)
	handle.flush()
	handle.close()
	
def process_blast_xml(fasta_id, xml_filename):
	global NUM_WITH_HIT, NUM_WITHOUT_HIT
	sanitize_xml(xml_filename)
	records = NCBIXML.parse(open(xml_filename))
	all_hits = [] # represents all the hits per fasta id
	for blast_records in records:
		for alignment in blast_records.alignments:
			for hsp in alignment.hsps:
				# ALTER THIS BASED ON HOW MANY DELIMITERS ARE IN THE BLAST DB
				# BELOW IS GOOD IFF TITLE IS: 
				# >gi|18411468|ref|NP_567197.1| a protein [A. thaliana]
				#					-- accn --  ------- desc ---------
				
				title = alignment.title.split('|')
				desc = title[-1]
				accn = title[-2]
				all_hits.append([str(accn), str(desc), float(hsp.expect), float(hsp.score)])
	
	# sort list so hits at different hsps are in-order
	all_hits = sorted(all_hits, key= lambda x: x[2])
	only_valid_hits = []
	for hit in all_hits: # for each hit...
		accn, desc, expect, score = hit
		is_ambigious_hit = False
		for ambig_word in AMBIGIOUS_KEYWORDS:
			if ambig_word in desc: # ... if an ambigious word is in the hit,
				is_ambigious_hit = True # flag and filter when all hits are done 
		
		# VALUES WHICH ARE POTENTIAL TOP-HITS GET STORED FOR SELECTION
		if not is_ambigious_hit:
			# to-string numerics for easy joining to string
			only_valid_hits.append([accn, desc, str(expect), str(score)])
	
	top_hit = ''
	if len(only_valid_hits) > 0:
		top_hit = '\t'.join(only_valid_hits[0])
		NUM_WITH_HIT+=1
	else:
		top_hit = '-\t-\t-\t-'
		NUM_WITHOUT_HIT+=1
	return top_hit.strip()

def run_blast(fasta_id, fasta_seq):
	fasta_file = create_fasta_file(fasta_id, fasta_seq)
	xml_file = DB_FOLDER+'/'+fasta_id+'.xml'
	cmd = BLAST_PROG+" -db " + DB_FOLDER+'/'+DB_NAME +\
		' -query ' + fasta_file + ' -evalue ' +str(BLAST_EVALUE) +\
		' -outfmt 5 ' + ' -out ' + xml_file + ' -num_alignments 10'
	os.system(cmd) # execute the local-blast command
	return fasta_id, fasta_file, xml_file
	
def parse_fasta():
	try:
		global NUM_ENTRIES
		seqs = {}
		header = ''
		for line in open(INPUT_FILE):
			line = line.strip()
			if '>' in line:
				
				header = line.replace('|', '_')
				seqs[header] = ''
			else:
				seqs[header]+=line
		NUM_ENTRIES = len(seqs)
		print('#/fasta entries:', NUM_ENTRIES)
		return seqs
	except IOError:
		print('[ERROR] input file:', INPUT_FILE, 'not found')
	

def execute_program():
	start_time = time.time()
	global DB_NAME, DB_FOLDER
	DB_NAME = os.path.basename(DB_FOLDER)
	DB_FOLDER = os.path.dirname(DB_FOLDER)
	
	# configure input file and blast db so output file can be saved correctly
	base_infile = os.path.basename(INPUT_FILE)
	base_infile = base_infile[0:base_infile.index('.')]
	base_db = os.path.basename(DB_NAME)
	base_db = base_db[0:base_db.index('.')]
	out_file = DB_FOLDER+'/'+base_infile+'__vs__'+base_db+'__'+'.txt'
	global OUT_HANDLE
	OUT_HANDLE = open(out_file, 'w')
	
	pool = multiprocessing.Pool(processes=NUM_PROCESSES)
	print(BLAST_PROG, 'LOCAL BLAST\t', DB_NAME)
	print('='*20,'#/worker processes:', NUM_PROCESSES,'='*20)
	dict_ids_seqs = parse_fasta()
	
	for fasta_id in dict_ids_seqs:
		seq = dict_ids_seqs[fasta_id]
		pool.apply_async(run_blast, args=(fasta_id.replace('>', ''), seq), callback=callback)
	pool.close()
	pool.join()
	pool.terminate()
	
	OUT_HANDLE.close()
	print('\n','='*20, BLAST_PROG,'of',NUM_ENTRIES,'seqs complete','='*20)
	print('total time: %.4f sec' % (time.time() - start_time),\
		'[ HIT:',NUM_WITH_HIT,'/',NUM_ENTRIES,';','NO HIT:',NUM_WITHOUT_HIT,'/',NUM_ENTRIES,']')
	print('output saved to [',out_file,']\n')

if __name__ == '__main__':
	in_file, out_file = '', '' 
	desc = 'A simple multiprocess local-blast python2.7 script'
	epilog = 'By Parsa Hosseini [August 2011]'
	parser = argparse.ArgumentParser(description=desc, epilog=epilog)
	parser.add_argument('-in', action='store', help='input fasta file', dest='in_file')
	parser.add_argument('-prog', action='store',help='blast program',dest='prog')
	parser.add_argument('-evalue', action='store', help='blast evalue', default='1e-5',dest='evalue')
	parser.add_argument('-num_proc', action='store', help='#/processes', default=8,dest='num_proc',type=int)
	parser.add_argument('-blast_db', action='store', help='db folder', dest='db_folder')
	
	parser = parser.parse_args()
	BLAST_EVALUE, BLAST_PROG, NUM_PROCESSES = parser.evalue, parser.prog, parser.num_proc
	INPUT_FILE = parser.in_file
	DB_FOLDER = parser.db_folder
	if all([BLAST_EVALUE, BLAST_PROG, NUM_PROCESSES,INPUT_FILE,\
			DB_FOLDER]):
		print(BLAST_EVALUE, os.path.isfile(BLAST_PROG), BLAST_PROG)
		execute_program()
	else:
		print('[ERROR] >=1 parameters are required:')
		print('-prog:', BLAST_PROG, '\n-evalue:', BLAST_EVALUE, '\n-in:', INPUT_FILE,\
			'\n-blast_db:', DB_FOLDER,'\n-num_proc:',NUM_PROCESSES)
	