# A trivial MAF+ alignment parser; a file-format produced using the 
# LASTZ alignment tool (Miller lab, PSU).

import argparse, csv, re

		
# Models a user-provided functional annotation
class GenomeAnnotations():
	def __init__(self, fname):
		self.fname = fname # annotation filename; provided CSV file
		self.annotations = {} # key => transcript ID, value => full annotation
		self.col_chrom = 0 # columns representing chromosome of annotation 
		self.col_start = 0 # columns representing start index of annotation
		self.col_end = 0 # columns representing end index of annotation
	
	# Sets arguments pertaining to the annotation such as start and end indices
	def set_args(self, args):
		self.col_chrom = args['chrom']
		self.col_start = args['start']
		self.col_end = args['end']
		
	# Parse the user-provided CSV file; use in-built module for parsing.
	def parse(self):
		linenum = 0
		for l in csv.reader(open(self.fname)):
			if linenum != 0: # the first line is the GFF3 comment
				# Represents CSV columns; zero-indexing
				chrom = l[self.col_chrom]
				start, end = int(l[self.col_start]), int(l[self.col_end])
				if chrom not in self.annotations:
					self.annotations[chrom] = []
				entry = {'chrom': chrom, 'start': start, 'end': end, 'annot': l}
				self.annotations[chrom].append(entry)
			linenum += 1 # increment line number
		print('Success:',len(self.annotations), 'chromosomes added')

# Represents an abstraction given a MAF+ entry  
class MAFEntry():
	def __init__(self):
		self.identity = 0
		self.coverage = 0
		self.continuity = 0
		self.cigar = ''
		self.score = 0
		self.query = None
		self.target = None

# Represents a specific chromosomal region that the alignment pertains to. 
class ChromosomeCoordinate():
	def __init__(self, data):
		self._data = data # references the actual data to be parsed
		self.align_len = 0 # length of alignment
		self.chrom_len = 0 # length of the actual chromosome
		self.coordinate = 0 # region where alignment hits the chromosome
		self.chrom = '' # chromosome
		self.strand = '' # references sequence strand; + or -
		self.seq = '' # references actual sequence masked or not
		self._parse() # given the raw data, parse and assign to attributes
	
	# Each entry is space-delimited; split and assign to each attribute
	def _parse(self):
		self._data = re.sub(' +', ' ', self._data).split(' ')
		self.chrom, self.coordinate = self._data[1], int(self._data[2])
		self.align_len, self.strand = int(self._data[3]), self._data[4]
		self.chrom_len, self.seq = int(self._data[5]), self._data[6]
		self._data = None # raw data not needed since it is already processed 

# Performs the bulk of the parsing.
class MAFParser():
	def __init__(self, fname):
		self.fname = fname # references MAF+ filename
		self.entries = [] # references all MAF entries
	
	# parse the provided MAF+ file
	def parse(self):
		data = open(self.fname)
		maf_accn = None
		for entry in data:
			entry = entry.strip()
			if entry.startswith('# identity'):
				maf_accn = MAFEntry()
				maf_accn.identity = entry
				self.entries.append(maf_accn)
			elif entry.startswith('# coverage'):
				maf_accn.coverage = entry
			elif entry.startswith('# continuity'):
				maf_accn.continuity = entry
			elif entry.startswith('# cigar'):
				maf_accn.cigar = entry
			elif entry.startswith('a score'):
				maf_accn.score = entry
			elif entry.startswith('s '):
				if not maf_accn.query: # set alignment coordinates to query
					maf_accn.query = ChromosomeCoordinate(entry)
				else: # set alignment coordinates to target
					maf_accn.target = ChromosomeCoordinate(entry)
		if len(self.entries) == 0:
			raise IOError('Error: No MAF+ entries identified.')
		else:
			print('Success:', len(self.entries), 'parsed')
			self.remove_metadata() # remove metadata from all MAF+ states
	
	# Removes metadata from the MAF+ identity attribute
	def _cleanse_identity(self, maf):
		# A typical MAF+ entry looks like this: # identity=64/73 (87.7%)
		# Only keep the percentage.
		perc = maf.identity[maf.identity.find('(')+1: maf.identity.find(')')]
		maf.identity = perc
	
	# Removes metadata from the MAF+ coverage attribute
	def _cleanse_coverage(self, maf):
		# A typical MAF+ entry looks like this: # coverage=62/366924 (0.0%)
		# Only keep the percentage.
		cov = maf.coverage[maf.coverage.find('(')+1: maf.coverage.find(')')]
		maf.coverage = cov
	
	# Removes metadata from the MAF+ continuity attribute
	def _cleanse_continuity(self, maf):
		# A typical MAF+ entry looks like this: # coverage=62/366924 (0.0%)
		# Only keep the percentage.
		c = maf.continuity[maf.continuity.find('(')+1: maf.continuity.find(')')]
		maf.continuity = c
	
	# Removes metadata from the MAF+ cigar attribute
	def _cleanse_cigar(self, maf):
		maf.cigar = maf.cigar.replace('# cigar=', '')
		
	# Removes metadata from the MAF+ score attribute
	def _cleanse_score(self, maf):
		maf.score = maf.score.replace('a score=', '')	
		
	# Each accession 
	def remove_metadata(self):
		for maf_obj in self.entries:
			self._cleanse_identity(maf_obj)
			self._cleanse_coverage(maf_obj)
			self._cleanse_continuity(maf_obj)
			self._cleanse_cigar(maf_obj)
			self._cleanse_score(maf_obj)

# Run the factory
def init(args):
	parser = MAFParser(fname=args['i'])
	parser.parse()
	query_annot = GenomeAnnotations(fname=args['annot_query'])
	query_annot.set_args(args)
	query_annot.parse()

if __name__ == '__main__':
	p = argparse.ArgumentParser(add_help=False) # parser object
	
	# Groups to encapsulate user-provided arguments
	p_reqd = p.add_argument_group('Required Arguments')
	p_opts = p.add_argument_group('Optional Arguments')
	
	# Required command-line arguments
	p_reqd.add_argument('--annot_query', help='Query CSV annotations [na]', 
					metavar='FILE', required=True)
	p_reqd.add_argument('-i', metavar='FILE', required=True, # provide input file
				help='MAF+ input file [na]', default=None)
	
	# Optional command-line arguments
	p_opts.add_argument('-chrom', help='Column representing chromosome IDs [3]', 
				default=3, type=int, metavar='INT')
	p_opts.add_argument('-start', help='Column representing start indices [4]', 
				default=4, type=int, metavar='INT')
	p_opts.add_argument('-end', help='Column representing end indices [5]', 
				default=5, type=int, metavar='INT')
	p_opts.add_argument('-h','--help', action='help',
					help='Show this help screen and exit')
	
	init(vars(p.parse_args())) # run the parser