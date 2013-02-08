import argparse, os
from collections import Counter

# A wrapper for BLASTP results. The file creating such an object is 2x columns
# large; the first being the query sequence while the latter is the homolog.
# These hits are tab-delimited.
class BLASTPResults():
	def __init__(self, fname):
		self.fname = fname
		self._is_file_valid() # determine if the file exists or not
		self.data = {} # key => source header, value => homolog
		self._parse() # parse user-provided file
	
	# Since the filename is large, only get the basename; much shorter
	def get_base_name(self):
		basename = os.path.basename(self.fname)
		index_dot = basename.rfind('.')
		return basename[0: index_dot] # do not return file extension
		
	# Determine if the query filename is valid or not
	def _is_file_valid(self):
		if os.path.exists(self.fname):
			return True
		else:
			raise OSError(self.fname +' is not a valid file [ERROR]')
	
	# Parse the user-provided file; is tab-delimited and 2x columns wide
	def _parse(self):
		for line in open(self.fname):
			key, value = line.strip().split('\t')
			if key in self.data: # if key already exists, exit
				raise IndexError(key + ' already exists; found >1 [ERROR]')
			else:
				self.data[key] = value # implies key is unique
		print(len(self.data), 'hits parsed for', self.get_base_name(), '[OK]')
	
	# Get BLASTP results
	def get_results(self):
		return self.data	
		
# Represents genes which are mutually conserved between two organisms.
# Two input files comprise this gene-set: a query and target. The query
# is such that the first column is the query organism gene, and its second
# column is its BLASTP homology. Similarly, the target sequence is the
# reciprocal. For a gene to be mutually conserved, its homology must be
# the same given the query and target.  
class MutuallyConservedGeneSet():
	def __init__(self, query, target):
		self.hits_query = BLASTPResults(query) # query BLASTP results
		self.hits_target = BLASTPResults(target) # target BLASTP results
	
	# Create handle to save output to
	def create_output_file(self):
		out = open(self.hits_query.get_base_name() + '_vs_' +\
			self.hits_target.get_base_name()+'_conserved.tab', 'w')
		out.write('Query\tTarget\n')
		out.flush()
		return out
	
	# Given parsed query and target hits, merge them to see if the homolog
	# of one is also the homolog of another. This reciprocal identity is
	# crucial for deducing magnitude of conservation.
	# For example: Suppose the query has the following:
	# Glyma15g23761.1	AT3G17450.1
	# If AT3G17450.1 also maps to Glyma15g23761.1 in the target-set, this
	# would render a conserved alignment as there is reciprocal homology.
	# The parameter 'whole_isoform' means entire accession is kept, otherwise
	# the accession is trimmed upto the first dot. Preceeding this is the gene.
	def join(self, whole_isoform):
		out = self.create_output_file() # create output file
		results = {}
		for query_key in self.hits_query.get_results():
			# Use value to see if it is a key in the target; implies conservation
			query_val = self.hits_query.get_results()[query_key]
			recip = '-' # reciprocal sequence
			if query_val in self.hits_target.get_results():
				# next, joining is based on the entire isoform, keep the 
				# accession (i.e. Glyma02g37560.1), otherwise trim upto the
				# last dot (i.e. Glyma02g37560) as this is the gene accession.
				if whole_isoform:
					recip = self.hits_target.get_results()[query_val] # isoform
				else:
					recip = isoform_to_gene(self.hits_target.get_results()[query_val])
					query_key = isoform_to_gene(query_key) # gene-query
					query_val = isoform_to_gene(query_val) # gene-homolog
			if query_key == recip:
				results[query_key] = query_val # save the conserved region and save later
				#print('Query:', query_key, 'Target:', query_val, 'Reciprocal:', recip)
		
		# A dictionary of query -> target hits are now present such that the 
		# reverse is true. Some targets though are isoforms and 
		# gene-representation produces duplicate entries.
		# eg. AAA.1 and AAA.2 yield AAA AAA, therefore duplicates.
		counts = Counter(results.values()) # target => frequency
		for query in results:
			if counts[results[query]] == 1:
				print('Query: '+query+'\tTarget: '+results[query])
				out.write(query+'\t'+results[query] + '\n')
				out.flush()
		print('Analysis complete given', self.hits_query.get_base_name(),
			'and', self.hits_target.get_base_name())
		
# Models which genes are conserved across all user-provided analysis. Each
# analysis is the pairwise alignment between a query and target organism.
class MultiSpeciesConservedSet():
	def __init__(self, files):
		self.files = list(set(files)) # remove duplicate filenames
		self.cons_geneset = {} # references genes conserved amongst all files
	
	# Find all the genes which are found within all the files
	def union(self):
		gene_set, no_hit_str = {}, None
		for file_no, fname in enumerate(self.files):
			for line in open(fname).readlines()[1:]: # skip first line; header
				line = line.strip().split('\t')
				gene, target = line[0], line[1]
				if gene not in gene_set:
					# create array however long the number of files are
					gene_set[gene] = [no_hit_str] * len(self.files)
				# set the target to the index of its file
				gene_set[gene][file_no] = target
		# remove entries which are not conserved across all species
		for i in gene_set:
			if gene_set[i].count(no_hit_str) == 0:
				self.cons_geneset[i] = gene_set[i]
					
	# Iterate through each file and see if the gene is conserved in other files.
	def output_conserved(self):
		out = open('conserved.tab', 'w')
		for counter, i in enumerate(self.cons_geneset):
			s = str(counter+1)+'\t' + i +'\t'+ '\t'.join(self.cons_geneset[i])
			out.write(s + '\n')
			out.flush()
		out.close()
		print('Analysis complete,', len(self.cons_geneset), 'conserved regions')

# Given an isoform accession, trim end-dot to yield gene name
def isoform_to_gene(accession):
	rindex_dot = accession.rfind('.') # find right-most dot (.)
	return accession[: rindex_dot] # slice upto the last-most dot

# Parses user-provided arguments
def create_parser():
	desc = 'Identify conserved genes across multi-species BLASTP hits'
	p = argparse.ArgumentParser(description = desc)
	p.add_argument('-query', metavar='', default=None,
		help='Query BLASTP results; 2 columns, tab-delimited [na]')
	p.add_argument('-target', metavar='', default=None,
		help='Target BLASTP results; 2 columns, tab-delimited [na]')
	p.add_argument('-merge', metavar='', nargs='+',
		help='Finds conserved genes across multiple GeneSet objects')
	p.add_argument('--whole_isoform', action='store_true', default=False,
					help='Perform conservation using isoform only [false]')	
	return vars(p.parse_args()) # hash of user-provided parameter arguments

if __name__ == '__main__':
	try:
		args = create_parser() # parse args and create a conserved gene-set
		if args['query'] and args['target']: # only if inputs have been provided
			gs = MutuallyConservedGeneSet(query=args['query'], target=args['target'])
			gs.join(whole_isoform=args['whole_isoform'])
		elif args['merge']: # if merging of inputs wishes to be performed
			conserved_set = MultiSpeciesConservedSet(files=args['merge'])
			conserved_set.union() # union of genes given all files 
			conserved_set.output_conserved() # output conserved genes to file
			
	except OSError as e:
		print(e)