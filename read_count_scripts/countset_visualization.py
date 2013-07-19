""" 
Useful script for the visualization of a CountSet file
"""

import argparse
from collections import OrderedDict

no_hit = '---' # symbol to represent no annotation is provided

class CountSet():
	"""
	A CountSet is a file which contains read-counts from a query
	and baseline CountBuild pair of files, hence its name. Each
	CountSet is then mapped onto user-provided annotations to aide
	deducing homology-based annotations of the differentially
	expressed genes 
	"""
	def __init__(self, fname):
		self.fname = fname # references raw filename
		self.col_rpkm = 0 # initially, no column yet reference RPKMs
		self.col_annot = 0 # references annotations to aide graphing
		self.counts = {} # K => annotation, V => induced, suppressed count
	
	# Set column-number referencing that of annotations to graph
	def set_col_annot(self, c):
		self.col_annot = c
	
	# Set column-number referencing that of RPKMs
	def set_col_rpkm(self, c):
		self.col_rpkm = c
	
	# Parse the user-provided CountSet file
	def count_frequency(self):
		for counter, line in enumerate(open(self.fname)):
			line = line.strip().split('\t')
			if counter > 0: # the first line is headers; skip this
				rpkm, annot = float(line[self.col_rpkm]), line[self.col_annot]
				if annot != no_hit: # ignore unannotated entries 
					for a in annot.split(no_hit):
						a = a.strip() # annotations may have multiple entries
						if a not in self.counts:
							# add number of induced or suppressed annotations
							self.counts[a] = {'ind': 0, 'supp': 0}
						if rpkm < 0:
							self.counts[a]['supp'] += 1
						else:
							self.counts[a]['ind'] += 1				
	
	# Filters counts given counts from both a query and baseline
	def filter_counts(self, diff):
		dict_filtered = {} # hash of filtered counts
		for i in self.counts:
			d = abs(self.counts[i]['ind'] - self.counts[i]['supp'])
			if d >= diff:
				dict_filtered[i] = self.counts[i]
		
		# once filtered, sort by name
		self.counts = OrderedDict()
		sorted_keys = sorted(dict_filtered.keys(), key=str.lower)
		for i in sorted_keys:
			self.counts[i] = dict_filtered[i]
		

		
# A high-level class which creates plots given CountSet objects.		
class VisualizationFactory():
	
	# Create a stacked bar-chart
	@staticmethod
	def stacked_barchart(cs):
		writer = open('plot_'+cs.fname, 'w')
		writer.write('Position\tLabel\tInduced\tSuppressed\n') # write header
		writer.flush()
		for pos, i in enumerate(cs.counts):
			# write output string
			s = str(pos+1) + '\t{' + i + '}\t' + str(cs.counts[i]['ind']) +\
				'\t' + str(cs.counts[i]['supp'])
			writer.write(s + '\n')
			writer.flush()
		writer.close()
		print('[ Analysis complete ]')
	
	@staticmethod
	def heatmap(*countset):
		pass

# Run analyses to create a stacked barchart given a CountSet
def run_barchart_factory(args):
	cset = CountSet(fname=args['file'])
	cset.set_col_annot(args['col_annot'])
	cset.set_col_rpkm(args['col_rpkm'])
	cset.count_frequency()
	cset.filter_counts(diff=args['diff'])
	VisualizationFactory.stacked_barchart(cs=cset)

def run_heatmap_factory(args):
	fnames_csets = args['files']
	csets = [] # list of CountSet objects
	for i in fnames_csets:
		cset = CountSet(i)
		cset.set_col_annot(args['col_annot'])
		cset.set_col_rpkm(args['col_rpkm'])
		cset.count_frequency()
		csets.append(cset)

if __name__ == '__main__':
	desc = 'Creates plots to help visualize annotations given RPKMs'
	p = argparse.ArgumentParser(description = desc)
	
	# Group of parameters pertaining to specific functionality
	p_bar = p.add_argument_group('Stacked stacked_barchart CountSet annotations')
	p_heatmap = p.add_argument_group('Heatmap of several CountSet')
	
	p_bar.add_argument('-file', metavar='FILE',
			help='Input CountSet file [na]')
	p_bar.add_argument('-diff', metavar='INT', default=6, type=int,
			help='Threshold difference given induced & suppressed counts [6]')
	p_bar.add_argument('-col_rpkm', metavar='INT', default=5, type=int,
			help='Column number for log2(query/baseline) RPKMs [5]')
	p_bar.add_argument('-col_annot', metavar='INT', default=22, type=int,
			help='Column number to visualize annotations [22]')
	
	# Parameters for creating a heatmap of multiple CountSets
	p_heatmap.add_argument('-files', metavar='', nargs='+',
			help='List of CountSet files [na]')
	
	args = vars(p.parse_args())
	if args['file']:
		run_barchart_factory(args)
	if args['files']:
		run_heatmap_factory(args)