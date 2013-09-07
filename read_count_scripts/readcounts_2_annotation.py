""" 
This script takes a query and baseline count-set and maps it to a
master annotation dataset. The purpose of this script is to deduce
differentially-expressed transcripts.
"""

import argparse, csv, math, os

class CountBuild():
	""" 
	A 'CountBuild' object is the product of computing how many reads mapped to 
	a SAM reference accession. This is computed using the command: 
	cat <samfile> | awk 'print {$3}' | sort | uniq -c
	Output is saved to a text file, able to be parsed as a CountSet object. 
	"""
	
	def __init__(self, fname, is_baseline):
		self.fname = fname # filename of this build
		self.numreads = 0 # total number of reads mapping to the build
		self.counts = {} # key => transcript ID, value => read count (int)
		self.basename = os.path.basename(fname) # get sole filename
		self.is_baseline = is_baseline # is this CountBuild a control
		
	# Parse a user-provided count-build; column 1 is an integer count, while
	# column 2 is the transcript accession
	def count_frequency(self):
		for line in open(self.fname):
			line = line.strip().split(' ')
			accn, count = line[1], int(line[0])
			self.counts[accn] = count # each accession references a count
			self.numreads += count # increment read-count
		print(self, len(self.counts), 'entries;', self.numreads, 'reads')

	# The current object is represented as the file basename
	def __str__(self):
		state = 'Baseline' if self.is_baseline else 'Query'
		return state + ' [' + os.path.basename(self.fname) + ']\t'

# Wraps a query and control CountBuild, producing a CountSet object; an object
# which references read-counts per transcript across both CountBuilds.
class CountSetFactory():
	""" 
	Performs all the logic which merges a control and query CountBuild,
	followed by addition of a user-provided set of annotations.
	"""
	
	def __init__(self, query, baseline):
		self.query = query
		self.baseline = baseline
		self.counts = {} # references counts given two CountBuild objects
		self.headers = ['Accession', 'Baseline_Num', 'Query_Num',
			'Baseline_RPKM', 'Query_RPKM', 'Log2_Query/Control_RPKM']
	
	# Merges both the query and baseline datasets so that transcript abundance
	# derivations can be performed easily.
	def merge(self):
		in_both, in_baseline, in_query = self._partition()
		
		# merge accessions found in both the query and control
		for a in in_both:
			count_query, count_baseline = self.query.counts[a], self.baseline.counts[a]
			self.counts[a] = {'query': count_query, 'baseline':count_baseline}		
		# if accession only in the baseline, set count of 1 to query accession
		for a in in_baseline:
			count_query, count_baseline = 1, self.baseline.counts[a]
			self.counts[a] = {'query': count_query, 'baseline':count_baseline}
		# if accession only in the query, set count of 1 to baseline accession
		for a in in_query:
			count_query, count_baseline = self.query.counts[a], 1
			self.counts[a] = {'query': count_query, 'baseline':count_baseline}
	
	# Partitions accession per CountBuild so that you have sets representing
	# accessions mapping to the query and baseline exclusively, and to both.
	def _partition(self):
		# create sets of both transcript IDs per countbuild for easy mapping
		accns_query = set([i for i in self.query.counts.keys()])
		accns_baseline = set([i for i in self.baseline.counts.keys()])
		
		# obtain accessions found in both CountBuilds as well as those which
		# are exclusive to each CountBuild object; 3x sets will be returned
		set_both = accns_query.intersection(accns_baseline) # accns in both
		set_only_query = accns_query.difference(accns_baseline) # only in query
		set_only_baseline = accns_baseline.difference(accns_query) # only in baseline
		
		# accessions are either found in both builds or exclusive to one
		return set_both, set_only_baseline, set_only_query
		
	# A CountSetFactory object can query its keys (transcript IDs) against a
	# user-provided annotation object (of type GenomeAnnotations)
	def add_annotations(self, annot, count_min, rpkm_min):
		# for each accession in the CountSetFactory object, get its annotation
		in_annot, not_in_annot = {}, []
		self.headers.extend(annot.headers) # add annotation headers 
		for accn in self.counts:
			if accn not in annot.annotations:
				not_in_annot.append(accn) # list of unannotated transcripts
			else:
				line = annot.annotations[accn] # get annotation for transcript
				start, end, accn =  int(line[annot.col_start]), int(line[annot.col_end]), line[annot.col_ids]
				accn_len = end - start # length of the transcript; for RPKM
				
				# compute rpkm for both query and baseline accessions
				rpkm_query = rpkm(c=self.counts[accn]['query'], 
							n=self.query.numreads, l=accn_len)
				rpkm_baseline = rpkm(c=self.counts[accn]['baseline'], 
							n=self.baseline.numreads, l=accn_len)
				log2_rpkm = math.log(rpkm_query/rpkm_baseline, 2) # log2 rpkm
				entry = [accn, str(self.counts[accn]['baseline']), 
					str(self.counts[accn]['query']),
					rpkm_baseline, rpkm_query, log2_rpkm] + line
				if abs(log2_rpkm) >= rpkm_min and\
                    self.counts[accn]['query'] >= count_min and\
                    self.counts[accn]['baseline'] >= count_min:
					in_annot[accn] = entry # transcript ID is found in annotations
		self.counts = in_annot
		
	# saves results to a file
	def write_countset(self):
		outfname = self.query.basename+'_vs_'+self.baseline.basename+'.txt'
		handle = open(outfname, 'w') # output file handle
		write('\t'.join(self.headers) + '\n', handle)
		for output in self.counts:
			outline = [str(i) for i in self.counts[output]] # cast to string
			write('\t'.join(outline) + '\n', handle)
		handle.close()
		
# Models a user-provided functional annotation
class GenomeAnnotations():
	def __init__(self, fname):
		self.fname = fname # annotation filename; provided CSV file
		self.annotations = {} # key => transcript ID, value => full annotation
		self.headers = [] # headers from the annotations file
		self.col_ids = 0 # columns representing annotation elements are unknown 
		self.col_start = 0
		self.col_end = 0
	
	# Sets arguments pertaining to the annotation such as start and end indices
	def set_args(self, args):
		self.col_ids = args['ids']
		self.col_start = args['start']
		self.col_end = args['end']
		
	# Parse the user-provided CSV file; use in-built module for such parsing.
	def count_frequency(self):
		linenum = 0
		for l in csv.reader(open(self.fname)):
			if linenum != 0: # the first line is the GFF3 comment
				# Represents GFF columns; zero-indexing
				self.annotations[l[self.col_ids]] = l # entire annotation row
			else:
				self.headers = l # set column headers
			linenum += 1 # increment line number
		print(self, len(self.annotations), 'annotations')
		
	# The current object is represented as the file basename
	def __str__(self):
		return 'Annotation [' + os.path.basename(self.fname) + ']\t'

# Computes RPKM given 'c': number of reads mapping to a transcript, 'n': the
# number of reads mapping to the entire build, and 'l': length of the transcript.
def rpkm(c, n, l):
	return (1e9 * float(c)) / (n * l)	

# Helpful function to write a string to an output file-handle
def write(outstr, handle):
	handle.write(outstr)
	handle.flush()

# Performs the bulk functionality of this script
def initialize(args):
	# count_frequency the user-provided annotations
	a = GenomeAnnotations(fname=args['a'])
	a.set_args(args={'start': args['start'], 'end': args['end'], 'ids': args['ids']})
	a.count_frequency() # count_frequency the user-provided CSV annotation
	
	# encapsulate user-provided RNA-Seq counts
	build_control = CountBuild(fname=args['baseline'], is_baseline=True)
	build_query = CountBuild(fname=args['query'], is_baseline=False)
	
	# count_frequency both count-sets
	build_control.count_frequency()
	build_query.count_frequency()
	
	# merge both CountBuilds into one so transcript abundance can be deduced
	countset = CountSetFactory(query=build_query, baseline=build_control)
	countset.merge()
	countset.add_annotations(annot=a, count_min=args['count'], rpkm_min=args['rpkm'])
	countset.write_countset() # last stage
	
def main():
	desc = 'Script to help map RNA-Seq read counts to functional annotations.'
	p = argparse.ArgumentParser(description=desc, add_help=False) # parser object
	
	# Groups to encapsulate user-provided arguments
	p_reqd = p.add_argument_group('Required Arguments')
	p_opts = p.add_argument_group('Optional Arguments')
	
	# Required command-line arguments
	p_reqd.add_argument('-a', help='CSV annotations [na]', metavar='FILE', 
				required=True)
	p_reqd.add_argument('-query', help='Query counts [na]', metavar='FILE',
				required=True)
	p_reqd.add_argument('-baseline', help='Control counts [na]', metavar='FILE',
				required=True)
	
	# Optional command-line arguments
	p_opts.add_argument('-rpkm', help='Minimum-acceptable log2 RPKM cutoff [0.5]', type=float, 
				default=0.5, metavar='FLOAT')
	p_opts.add_argument('-count', help='Minimum-acceptable read count [5]', 
				default=5, type=int, metavar='INT')
	p_opts.add_argument('-ids', help='Column representing transcript IDs [1]', 
				default=1, type=int, metavar='INT')
	p_opts.add_argument('-start', help='Column representing start indices [4]', 
				default=4, type=int, metavar='INT')
	p_opts.add_argument('-end', help='Column representing end indices [5]', 
				default=5, type=int, metavar='INT')
	p_opts.add_argument('-h','--help', action='help',
					help='Show this help screen and exit')
	initialize(vars(p.parse_args())) # begin analysis given arguments

if __name__ == '__main__':
	main()