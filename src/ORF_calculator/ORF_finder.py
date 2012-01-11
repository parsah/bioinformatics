
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML

def get_ORF(sequence, title, blast_ncbi=False):
	sequence = Seq(sequence).upper()
	dict_frame_and_protein = {}
	
	# Get Open Reading Frames on the forward strand (+)
	for frame in range(3):
		seq_in_frame = sequence[frame:]
		frame_number_header = 'Frame #'+str(frame+1)
		num_spaces = (len(frame_number_header)+2) *' '
		protein_seq = seq_in_frame.translate(to_stop=False)
		dict_frame_and_protein['frame_'+str(frame+1)] = str(protein_seq)
		
	# Get Open Reading Frames on the r strand (-)
	sequence = sequence.reverse_complement()
	for frame in range(3):
		seq_in_frame = sequence[frame:]
		frame_number_header = 'Frame #-'+str(frame+1)
		num_spaces = (len(frame_number_header)+2) *' '
		protein_seq = seq_in_frame.translate(to_stop=False)
		dict_frame_and_protein['frame_-'+str(frame+1)] = str(protein_seq)
		
	out_str = output_results(dict_frame_and_protein, title,blast_ncbi)
	return out_str
	
def output_results(dict_frame_and_protein, title, blast_ncbi):
	out_str = ''
	for each_orf in dict_frame_and_protein:
		out_str += each_orf+'\t'+dict_frame_and_protein[str(each_orf)]
	sorted_dict = sorted(dict_frame_and_protein.items(), key=lambda x: get_protein_length(x[1]), reverse=True)
	for i in sorted_dict:
		frame_number, protein_seq = i
		out_str+= str(i[0])+'\t'+ str(get_protein_length(i[1]))
	collection_orfs = fetch_orfs_with_met(dict_frame_and_protein, title)	
	if blast_ncbi:
		ncbi_blast(protein_seq)
	
	return collection_orfs

def ncbi_blast(protein_seq):
	result_handle = NCBIWWW.qblast('blastp', 'nr', protein_seq, expect=1e-5)
	save_file = open("my_blast.xml", "w")
	save_file.write(result_handle.read())
	save_file.close()
	
	result_handle = open("my_blast.xml")
	blast_records = NCBIXML.parse(result_handle)
	
	num_hits = 0
	for each_record in blast_records:
		for alignment in each_record.alignments:
			for hsp in alignment.hsps:
				title = alignment.title
				title = title.split('|')[4].strip()
				title = title.split('>')[0].strip()
				out_str = title+'\t'+str(alignment.length)+'\t'+str(hsp.score)+'\t'+str(hsp.expect)
				num_hits+=1
				if num_hits == 0:
					print '\t'+out_str
					
def fetch_orfs_with_met(dict_frame_and_protein, title):
	# for all the 6x open-reading frames, you do not want to simply blast each
	# reading frame because there are many MET start-sites within this large
	# sequence. What you want to do is split this large reading frame sequence
	# into sub-sequences whereby each sub-sequence starts with a MET start-site.
	entire_collection = ''
	for frame_number in dict_frame_and_protein:
		frame_protein_seq = dict_frame_and_protein[frame_number]
		# split the larger sequence into wherever a MET start-site is found.
		potential_orfs = frame_protein_seq.split('M')
		
		# for each site which starts with a MET, slice the sequence to the 
		# first occurance of a stop site. This slice will represent a candidate
		# open reading frame within the larger open-reading frame.
		
		counter_sub_orf = 0
		for each_potential_orf in potential_orfs:
			index_first_stop = each_potential_orf.find('*')
			sub_orf = 'M'+each_potential_orf[0:index_first_stop]+'*'
			
			if len(sub_orf) >= 41:
				counter_sub_orf+=1
				header = '>'+title+'_'+frame_number+'_'+'sub_orf_'+str(counter_sub_orf)+'_'+str(len(sub_orf))+'bp'
				sub_orf
				
				fasta_str = header+'\n'+sub_orf+'\n'
				
				entire_collection = entire_collection+fasta_str
	return entire_collection			
	
def get_protein_length(protein_seq):
	raw_protein_length = len(protein_seq)
	num_stop_codon = str(protein_seq).count('*')
	return raw_protein_length - num_stop_codon
