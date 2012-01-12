'''
A script which takes a series of inputs:
* A user-provided fasta file of input sequences
* A reference genome
Given the set of inputs above, each is mapped against the reference genome 
(using BLAST), and the location is then derived. Upsteam and downstream of this
hit-site is then sliced out. These two set of slices are then compared against
a dictionary of binding-elements (BEs) to derive which BEs exist and which do
not.
'''

import argparse
from Bio.Blast import NCBIXML
import multiprocessing
import os
import tempfile
import time

# words to filter from annotations 
AMBIGIOUS_KEYWORDS = ['hypothetical', 'putative', 'unknown', 'unnamed', 'predicted', 
					'uncharacterized']
HITS = {}
CURR_BLAST_HIT = 0


# represents all the files and states pertaining to the analysis
class Analysis_Bucket():
	def __init__(self, fasta_file, ref, be_db, bp_up, bp_down, evalue, num_proc, executable):
		self.fasta_file = fasta_file
		self.ref_genome = ref
		self.be_db = be_db
		self.bp_up = bp_up
		self.bp_down = bp_down
		self.evalue = evalue
		self.num_proc = num_proc
		self.blast_exec = executable
		self.blast_exec_dir = os.path.dirname(self.blast_exec)+'/'
		self.output_temp_dir = tempfile.mkdtemp(dir=os.getcwd())
		self.num_entries = 0
	
	@staticmethod
	def new_segment_line():
		return '\n'+'*'*40
		
	def __str__(self):
		return 'Input parameters:\n'+\
			'Binding Element Search Script - ['+str(time.asctime())+']\n'+\
			'Fasta file: '+str(self.fasta_file)+'\n'+\
			'Reference genome: '+str(self.ref_genome)+'\n'+\
			'Binding Element (BE) DB: '+str(self.be_db)+'\n'+\
			'up-stream bp: '+str(self.bp_up)+'\n'+\
			'down-stream bp: '+str(self.bp_down)+'\n'+\
			'BLAST e-value: '+str(self.evalue)+'\n'+\
			'BLAST location: '+str(self.blast_exec)+'\n'+\
			'#/processes: '+str(self.num_proc)+'\n'+\
			Analysis_Bucket.new_segment_line()
	
	
	def is_param_valid(self):
		print '--- Checking validity of input files: ---'
		is_valid = os.path.exists(self.fasta_file)
		# check to see if fasta file is found
		if is_valid is False:
			print '[ERROR: Fasta file (',self.fasta_file,') is not found]'
		
		# check to see if reference genome is found
		is_valid = os.path.exists(self.ref_genome)
		if is_valid is False:
			print '[ERROR: Reference genome (',self.ref_genome,') is not found]'
		
		# check to see if binding-element database is found
		is_valid = os.path.exists(self.be_db)
		if is_valid is False:
			print '[ERROR: Binding Element DB (',self.be_db,') is not found]'
		
		# stdout status of validity
		if is_valid:
			print 'VALID: Input parameters all valid. Continuing-on with analysis.'	
		else:
			print 'INVALID: Input parameters invalid. Please refer to message above.'	
		return is_valid

	def parse_fasta(self, fasta_file):
		print Analysis_Bucket.new_segment_line()
		print '--- Parsing fasta file and beginning analysis ---'
		entries = {}
		for line in open(self.fasta_file):
			line = line.strip()
			if len(line) == 0:
				break
			else:
				if '>' in line:
					header = line
					entries[header] = ''
				else:
					entries[header]+=line
		
		self.num_entries = len(entries)
		print '#/fasta entries parsed:', self.num_entries
		return entries
	
	def process_blast_results(self):
		global HITS
		print '\nParsing BLAST results:'
		for counter in HITS:
			top_hit = HITS[counter]
			start, end, chrom = top_hit['start'], top_hit['end'], top_hit['chrom']
			accn, evalue = top_hit['accn'], top_hit['e_val']
			score, fasta = top_hit['score'], top_hit['fasta']
			if end < start: # reverse strand (-)
				print counter, ')' ,fasta, 'maps to chrom.', chrom, 'from',\
					start, '->', end, '[-]\t', start-end,'bps.'
			else: # forward strand (+)
				print counter, ')' ,fasta, 'maps to chrom.', chrom, 'from',\
					start, '->', end, '[+]', end-start, 'bps.'	
					
	def ref_genome_split(self):
		filename, header = '', ''
		print self.output_temp_dir
		for line in open(self.ref_genome):
			line = line.strip()
			if '>' in line:
				line = line.replace('>', '')
				header = line
				filename = open(self.output_temp_dir+'/'+line, 'w+')
				filename.write('>'+header+'\n')
				print 'added', header
			else:
				filename.write(line+'\n')
				filename.flush()
			
# for each fasta entry, blast against reference genome
def run_analysis(obj_analysis):
	# first, parse the fasta file
	entries = obj_analysis.parse_fasta(obj_analysis.fasta_file)
	pool = multiprocessing.Pool(processes=obj_analysis.num_proc)
	print '#/processes created:', obj_analysis.num_proc
	for each_fasta in entries:
		seq = entries[each_fasta]
		header = each_fasta.replace('>', '')
		pool.apply_async(func=run_blast, args=(header, seq, obj_analysis), callback=callback_stdout)
		
	pool.close()
	pool.join()
	pool.terminate()

def callback_stdout(return_from_blast):
	global CURR_BLAST_HIT, HITS
	CURR_BLAST_HIT+=1
	fasta_id = return_from_blast['fasta_id']
	obj_analysis = return_from_blast['curr_obj']
	xml_filename = return_from_blast['xml_file']
	fasta_filename = return_from_blast['fasta_file']
	
	blast_results = process_blast_xml(fasta_id, xml_filename)
	# store blast result as hash; pointed by its sequence ID number.
	HITS[CURR_BLAST_HIT] = blast_results
	print CURR_BLAST_HIT,'/', obj_analysis.num_entries,'[OK]' # progress output
	
	# delete xml and fasta filenames
	os.remove(xml_filename)
	os.remove(fasta_filename)

def process_blast_xml(fasta_id, xml_filename):
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
				desc = desc.split(' ')[-1]
				accn = title[-2]
				all_hits.append([str(accn), str(desc), float(hsp.expect),
					float(hsp.score), int(hsp.sbjct_start), int(hsp.sbjct_end), fasta_id])
	
	# sort list so hits at different hsps are in-order
	all_hits = sorted(all_hits, key= lambda x: x[2])[0]
	all_hits = dict(zip(['accn', 'chrom', 'e_val', 'score', 'start', 'end', 'fasta'], all_hits))
	
	return all_hits # return hit with best e-value

def run_blast(fasta_id, fasta_seq, obj_analysis):
	fasta_id = fasta_id.replace('|', '_')
	blast_exec_dir = obj_analysis.blast_exec_dir
	fasta_file = create_fasta_file(fasta_id, fasta_seq, obj_analysis)
	xml_file = obj_analysis.output_temp_dir+'/'+fasta_id+'.xml'
	cmd = blast_exec_dir+"blastn -db " + obj_analysis.ref_genome +\
		' -query ' + fasta_file + ' -evalue ' +str(obj_analysis.evalue) +\
		' -outfmt 5 ' + ' -out ' + xml_file + ' -num_alignments 3'
	os.system(cmd) # execute the local-blast command
	return {'fasta_id': fasta_id, 'curr_obj': obj_analysis, 'xml_file':xml_file, 'fasta_file': fasta_file}

def create_fasta_file(fasta_id, fasta_seq, obj_analysis):
	filename = obj_analysis.output_temp_dir+'/'+fasta_id+'.fasta'
	out_handle = open(filename, 'w')
	out_handle.write('>'+fasta_id[1:]+'\n'+str(fasta_seq)+'\n')
	out_handle.flush()
	out_handle.close()
	return filename

if __name__ == '__main__':
	str_epilog = 'Ex: python driver.py Fasta_file Ref_genome Be_db -u 2000 -d 200'
	str_desc = 'Drives the ability to find binding-elements (BEs) in a set of'+\
		' sequences. Given a an populated input fasta file, each sequence is'+\
		' mapped against a reference such that up-stream and down-stream regions'+\
		' are analyzed for the presence of BEs. Such BEs are from a user-provided'+\
		' dictionary.'
	p = argparse.ArgumentParser(description=str_desc, epilog=str_epilog)
	# required params
	p.add_argument('-i','-input', help='Input fasta file; fasta [REQUIRED]',required=True)
	p.add_argument('-r','-ref', help='Reference genome BLAST DB; fasta [REQUIRED]', required=True)
	p.add_argument('-b','-be', help='Binding Element DB [REQUIRED]', required=True)
	p.add_argument('-x','-blastn', help='BLASTN executable [REQUIRED]',required=True)
	
	# optional params
	p.add_argument('-e','-evalue', help='BLAST e-value <def: 1e5>', default=1e5)
	p.add_argument('-n','-num_proc', help='#/processes <def: 2>', default=2,type=int)	
	p.add_argument('-u', '-up', help='up-stream bp <def: 2000>', default=2000, type=int)
	p.add_argument('-d', '-down', help='down-stream bp <def: 200>', default=200, type=int)
	
	args = p.parse_args()
	analysis_obj = Analysis_Bucket(fasta_file=args.i, ref=args.r, be_db=args.b, 
			bp_up=args.u, bp_down=args.d, evalue=args.e, num_proc=args.n, 
			executable=args.x)
	print analysis_obj
	if analysis_obj.is_param_valid(): # if input is valid, process or stop
		run_analysis(analysis_obj)
	
	
	# split genome into reference chromosomes; easy parsing
	analysis_obj.ref_genome_split()
	
	# when complete, process results and map to corresponding chromosome
	analysis_obj.process_blast_results()
	
	os.rmdir(analysis_obj.output_temp_dir)
	