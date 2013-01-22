# A trivial MAF+ alignment parser; a file-format produced using the 
# LASTZ alignment tool (Miller lab, PSU).

import argparse

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
	def __init__(self):
		self.align_len = 0
		self.coordinate = 0
		self.chrom = ''
		self.strand = ''
		self.seq = ''

# Performs the bulk of the parsing.
class MAFParser():
	def __init__(self, fname):
		self.fname = fname # references MAF+ filename
		self.entries = [] # references all MAF entries
		self._parse()
	
	# parse the provided MAF+ file
	def _parse(self):
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
					maf_accn.query = ChromosomeCoordinate()
				else: # set alignment coordinates to target
					maf_accn.target = ChromosomeCoordinate()
		if len(self.entries) == 0:
			raise IOError('Error: No MAF+ entries identified.')
		else:
			print('Success:', len(self.entries), 'parsed')
					
	# Each accession 
	def _cleanse(self):
		pass
	
if __name__ == '__main__':
	p = argparse.ArgumentParser()
	p.add_argument('-i', metavar='FILE', required=True, # provide input file
				help='MAF+ input file [na]', default=None)
	args = vars(p.parse_args())
	MAFParser(fname=args['i'])