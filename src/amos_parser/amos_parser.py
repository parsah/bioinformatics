
import linecache

IN_DIR = 'DIR CONTAINING AMOS AFG FILE'

def get_line_from_file(line_number):
	contig_id = linecache.getline(IN_DIR+'velvet_asm.afg', line_number+2)
	contig_id = contig_id.replace('iid:', '').strip()
	return contig_id

def parse():
	current_contig = None
	if_for_src = 'src:'
	contig_id_and_seqs = {}
	
	handle = open(IN_DIR+'velvet_asm.afg', 'r')
	for line_number, each_line in enumerate(handle):
		each_line = each_line.strip()
		if '{CTG' in each_line:
			current_contig = get_line_from_file(line_number)
			contig_id_and_seqs[current_contig] = []
#			print '\n',current_contig, 'added to dict'
		
		if current_contig:
#			print 'Contig #:'+current_contig+'\t'+str(line_number)+'\t'+each_line
			if if_for_src in each_line:
				src_mapping_2_contig = int(each_line.replace(if_for_src, ''))
				contig_id_and_seqs[current_contig].append(src_mapping_2_contig)
	
	out_writer = open(IN_DIR+'count_reads_per_contig.txt', 'w')
	header = 'Ctg #\t#/reads'
	out_writer.write(header+'\n')
	out_writer.flush()
	print header
	for i in contig_id_and_seqs:
		print i, '\t', len(contig_id_and_seqs[i])
		out_writer.write(str(i)+'\t'+str(len(contig_id_and_seqs[i]))+'\n')
		out_writer.flush()
	out_writer.close()
	print 'analysis done'

if __name__ == '__main__':
	parse()