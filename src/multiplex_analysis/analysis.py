
import os, random, subprocess, time

# information for running alignment against a reference
BWA_EXEC = 'bwa-0.5.9/bwa'
# REF points to reference genome fasta file
REF = '/unmasked_genome/ref.fasta'

# annotation-based information for building binary tree
# variable 'ANNOTATIONS' will be split based on indices of the chromosome, start
# and end
ANNOTATIONS = '/Annotations/annots.txt'
ANNOTATION_HEADERS = ''
IDX_CHROM, IDX_START, IDX_END = 3, 5, 6

class Chromosome_Binary_Search_Tree():
	def __init__(self, applicable_chrom):
		self.chromosome = applicable_chrom
		self.numNodes = 0
		self.rootNode = None
		self.tag_count_results = []
		
	def insert(self, start, end, chrom, annotation):
		new_node = Binary_Search_Tree_Node(start, end, chrom, annotation)
		if self.rootNode is None:
			self.rootNode = new_node
		else:
			curr_node = self.rootNode
			while True:
				if start < curr_node.start:
					if curr_node.leftNode is None:
						curr_node.leftNode = new_node
						new_node.parentNode = curr_node
						break
					curr_node = curr_node.leftNode
				else:
					if curr_node.rightNode is None:
						curr_node.rightNode = new_node
						new_node.parentNode = curr_node
						break
					curr_node = curr_node.rightNode
		self.numNodes+=1
	
	def find_node(self, node, index):
		if not node:
			return None
		if index < node.start:
			return self.find_node(node.leftNode, index)
		elif index > node.end:
			return self.find_node(node.rightNode, index)
		else:
			node.tag_count+=1
			return str(node.start)+'_'+str(node.end)

	def infix(self, node):
		if node.leftNode:
			self.infix(node.leftNode)
		self.tag_count_results.append([node.tag_count, node.annotation])
		if node.rightNode:
			self.infix(node.rightNode)
	
class Binary_Search_Tree_Node():
	def __init__(self, start, end, chromosome, annotation):
		self.start =  start
		self.end = end
		self.chromosome = chromosome
		self.annotation = annotation
		self.tag_count = 0
		self.rightNode = None
		self.leftNode = None
		self.parentNode = None
	
	def __str__(self):
		return str(self.start)+' '+str(self.end)+' '+str(self.chromosome)

class Analysis():
	def __init__(self, index, num_entries, fasta_file):
		self.index = index
		self.in_file = fasta_file
		self.num_entries = num_entries
		self.seqs = None # #/seqs in a specific dataset
		self.out_dir = None
		self.str_num_seqs = '' # formatted, stringified ver of len(self.seqs)
		
	# count the number of non-empty entries in the file
	@staticmethod
	def get_entries(in_file):
		list_datasets = [i.strip() for i in open(in_file) if len(i) > 1]
		return list_datasets
	
	def analyze(self):
		# phase 1 - parsing input
		# for each dataset, parse its seqs be-it fasta or fastq
		in_file = self.in_file.strip()
		print '#'*4, self.index+1, '/', self.num_entries,'-', self.in_file,'#'*4
		a_dset_extension = os.path.splitext(self.in_file)[-1]
		if a_dset_extension == '.fasta':
			self.parse_fasta(fasta = self.in_file)
		else:
			print '\tfastq' # implement fastq parser
		
		# phase 2 - read mapping
		# first, make a folder per dataset and store all reads in there
		start_time, self.out_dir = time.time(), os.path.dirname(in_file)+'/output/'
		if not os.path.isdir(self.out_dir):
			os.mkdir(self.out_dir)
		sam_file = self.out_dir+'reads.sam'
		bwa_cmd = BWA_EXEC+' bwasw '+ REF+' '+self.in_file +' -t 8 -f ' + sam_file
		print 'mapping', self.str_num_seqs, 'seqs to genome...',
		
		fnull = open(os.devnull, 'w') # suppress stdout output
		subprocess.call(bwa_cmd, shell = True, stdout = fnull, stderr = fnull)
		fnull.close() # close stream since cmd execution is now complete
		total_time = format(time.time() - start_time,'.2f')
		print '[OK]\n', '[',total_time,'s]\n'
		
		# phase 3 - process sam file, extract reads which map to the 
		# reference once, multiple times or neither (i.e. potential SR seqs)
		# finally, save such reads to file
		seqs_map, seqs_not_map = self.process_sam_file(sam_file)
		out_handle_reads_map2_genome = open(self.out_dir+'reads_map_2_genome.fasta', 'w')
		out_handle_reads_not_map2_genome = open(self.out_dir+'reads_not_map_2_genome.fasta', 'w')
		out_str = ''
		for read in seqs_map:
			score, chrom = seqs_map[read]['score'], seqs_map[read]['chrom']
			index = seqs_map[read]['idx']
			out_str = '>'+read+'|'+chrom+'|'+str(score)+'|'+str(index)+\
				'\n'+self.seqs[read]
			out_handle_reads_map2_genome.write(out_str+'\n')
			out_handle_reads_map2_genome.flush()
		for read in seqs_not_map:
			out_str = '>'+read+'\n'+seqs_not_map[read]
			out_handle_reads_not_map2_genome.write(out_str+'\n')
			out_handle_reads_not_map2_genome.flush()
		
		out_handle_reads_map2_genome.close()
		out_handle_reads_not_map2_genome.close()
		return seqs_map, seqs_not_map
	
	def process_sam_file(self, sam_input):
		print 'processing SAM file, extracting mapped/non-mapped reads...',
		seqs_map, seqs_not_map, start_time = {}, {}, time.time()
		num_no_map, num_map_many_chrom = 0, 0
		for line in open(sam_input):
			line = line.strip().split('\t')
			if len(line) >= 11: # 11, 16 == #/columns which do not/do hit ref.
				read_id, sequence, chrom = line[0], line[9], line[2]
				if chrom == '*': # if the read did not hit the reference
					seqs_not_map[read_id] = sequence
					num_no_map+=1
				else: # otherwise, the read mapped to the chromosome
					aln_sc, idx = int(line[11].strip().split(':')[-1]), line[3]
					# if the read only hit 1 unique location, store it ...
					if read_id not in seqs_map:
						seqs_map[read_id] = {'idx':idx, 'score':aln_sc,'chrom':chrom}
					else: # otherwise, take the better alignment and replace
						# the alignment which is already in-place. 
						old_score = seqs_map[read_id]['score']
						if aln_sc > old_score:
							seqs_map[read_id]['idx'] = idx
							seqs_map[read_id]['score'] = aln_sc
							seqs_map[read_id]['chrom'] = chrom
						num_map_many_chrom+=1
		print '[OK]'
		# display #/reads which map to reference
		perc_map = '[ '+format(float(len(seqs_map))/len(self.seqs)*100, '.2f')+'%]'
		num_map = format(len(seqs_map), ',d')+' / '+self.str_num_seqs
		# display #/reads which map to multiple regions on reference
		perc_multi_map = '[ '+format(float(num_map_many_chrom)/len(self.seqs)*100, '.2f')+'%]'
		num_multi_map = format(num_map_many_chrom, ',d')+' / '+self.str_num_seqs
		# display #/reads which do not map to reference
		perc_not_map = '[ '+format(float(len(seqs_not_map))/len(self.seqs)*100, '.2f')+'%]'
		num_not_map = format(len(seqs_not_map), ',d')+' / '+self.str_num_seqs
		
		print '#/seq mapping to reference (hit 1x):', num_map, perc_map
		print '#/seq not mapping to reference:', num_not_map, perc_not_map
		print '#/seq mapping to reference (> 1x):', num_multi_map, perc_multi_map
		print '[',format(time.time() - start_time, '.2f'),'s]'
		return seqs_map, seqs_not_map
		
	def parse_fasta(self, fasta):
		print 'parsing FASTA reads...',
		header, all_seqs, start_time = '', {}, time.time()
		for line in open(fasta):
			line = line.strip()
			if '>' in line:
				header = line[1:]
				all_seqs[header] = ''
			else:
				all_seqs[header]+= line
		total_time = format(time.time() - start_time,'.2f')
		print '[OK]\n', '[',total_time,'s]', 
		print format(len(all_seqs), ',d'),'seqs -', os.path.basename(fasta),'\n'
		self.seqs = all_seqs
		self.str_num_seqs = format(len(self.seqs), ',d')
	
	def tree_parse_annotations(self):
		# parse all the annotations and store them in a simple bin. search tree
		# each chromosome in the annotations-file (see IDX_CHROM) represents
		# its own tree, with annotations in that tree each representing either
		# a L or R child.
		start_time = time.time()
		print '\nparsing annotations...',
		counter, dict_all_nodes_per_chrom = 0, {}
		for line in open(ANNOTATIONS):
			line = line.strip().split('\t')
			if counter == 0:
				global ANNOTATION_HEADERS
				ANNOTATION_HEADERS = '\t'.join(line) # goes to header of BST
			else:
				stringified_line = '\t'.join(line)
				chrom, start, end = line[IDX_CHROM], line[IDX_START], line[IDX_END]
				
				if not chrom in dict_all_nodes_per_chrom:
					dict_all_nodes_per_chrom[chrom] = []
				value = {'chrom': chrom, 'start': int(start), 'end':int(end), 
						'annot':stringified_line}
				
				dict_all_nodes_per_chrom[chrom].append(value)
			counter+=1
		
		# next, shuffle each chromosomes nodes to randomize the list and make
		# it fit for storing as a binary-search tree
		for chrom in dict_all_nodes_per_chrom:
			random.shuffle(dict_all_nodes_per_chrom[chrom])
		total_time = format(time.time() - start_time, '.2f')
		print '[OK]', len(dict_all_nodes_per_chrom), 'chromosome(s) parsed.'
		print '[',total_time,'s]'
		return dict_all_nodes_per_chrom

	def tree_generate(self, dict_all_nodes_per_chrom):
		dict_binary_trees = {}
		print '\ngenerating BST per chromosome...',
		# generate a 'Chromosome_Binary_Search_Tree' per chromosome
		for chrom in dict_all_nodes_per_chrom:
			tree = Chromosome_Binary_Search_Tree(chrom)
			# for each tree, add all its nodes to the tree
			for each_entry in dict_all_nodes_per_chrom[chrom]:
				node_start, node_end = each_entry['start'], each_entry['end']
				node_chrom, node_annot = each_entry['chrom'], each_entry['annot']
				tree.insert(node_start, node_end, node_chrom, node_annot)
			dict_binary_trees[chrom] = tree
		print '[OK]'
		return dict_binary_trees
	
	def perform_tag_counting(self, trees, mapped_seqs):
		# perform tag-counting given a set of all the reads mapping to the
		# respective chromosome
		handle_out_map_2_genome_and_transcript = open(self.out_dir+\
				'reads_map_2_genome_and_transcript.fasta', 'w')
		handle_out_map_2_genome_not_transcript = open(self.out_dir+\
				'reads_map_2_genome_not_transcript.fasta', 'w')
		num_map_2_tscript, num_not_map_2_tscript, start_time = 0, 0, time.time()
		print 'tag-counts', format(len(mapped_seqs), ',d'), 'reads...',
		for mapped_read in mapped_seqs:
			chrom = mapped_seqs[mapped_read]['chrom']
			index = mapped_seqs[mapped_read]['idx']
			score = mapped_seqs[mapped_read]['score']
			applicable_tree = trees[chrom]
			# either the read maps to a transcript or it doesn't.
			hit = applicable_tree.find_node(applicable_tree.rootNode, int(index))
			
			out_str = '>'+mapped_read+'|'+chrom+'|'+str(score)+'|'+\
				str(index)+'|'+str(hit)+'\n'+self.seqs[mapped_read]
			if hit: # if the read maps to a transcript
				handle_out_map_2_genome_and_transcript.write(out_str+'\n')
				handle_out_map_2_genome_and_transcript.flush()
				num_map_2_tscript+=1
			else: # the read maps to the genome but not to a transcript
				handle_out_map_2_genome_not_transcript.write(out_str+'\n')
				handle_out_map_2_genome_not_transcript.flush()
				num_not_map_2_tscript+=1
		
		handle_out_map_2_genome_not_transcript.close()
		handle_out_map_2_genome_and_transcript.close()
		print '[OK]\n[', format(time.time() - start_time, '.2f'), 's]'
		
		print '#/seq mapping to genome and transcript:',\
			format(num_map_2_tscript, ',d'),'/', format(len(mapped_seqs),',d')
		print '#/seq mapping to genome but not to transcript:',\
			format(num_not_map_2_tscript,',d'),'/', format(len(mapped_seqs),',d')
		self.get_tag_counts(trees)
		
	def get_tag_counts(self, trees):
		# all tag_counts are stored as individual trees.
		# this list stores them as one complete list for easy parsing to file.
		all_tag_counts = []
		for each_chrom in trees:
			tree = trees[each_chrom]
			tree.infix(tree.rootNode) # saves tree results as list
			all_tag_counts.extend(tree.tag_count_results)
		out_handle_tag_counts = open(self.out_dir+'tag_counts.txt', 'w')
		out_handle_tag_counts.write('Tag_count\t'+ANNOTATION_HEADERS+'\n')
		out_handle_tag_counts.flush()
		for i in all_tag_counts:
			tag_count, annotation = str(i[0]), i[1]
			out_handle_tag_counts.write(tag_count+'\t'+annotation+'\n')
			out_handle_tag_counts.flush()
		out_handle_tag_counts.close()
		print '\n[ANALYSIS COMPLETE]\n'
		
if __name__ == '__main__':
	# provide a file which points to the fasta files. fastq not yet implemented
	list_entries = Analysis.get_entries('fasta_files.txt')
	for counter, each_fasta in enumerate(list_entries):
		analysis_runner = Analysis(counter, len(list_entries), each_fasta)
		dict_seqs_map, dict_seqs_not_map = analysis_runner.analyze()
		dict_all_nodes_per_chrom = analysis_runner.tree_parse_annotations()
		trees = analysis_runner.tree_generate(dict_all_nodes_per_chrom)
		analysis_runner.perform_tag_counting(trees, dict_seqs_map)
		