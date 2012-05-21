
'''
A simple script which takes two files:
1) A line-separated file of only accessions, A
2) A tab-delimited file such that column 1 maps to items to file A. Each
consecutive column thereof represents BLAST homology-based functional 
annotations for the respective annotation
'''

from optparse import OptionParser
import os

def parse_accessions(filename):
	# traverse through each accession and simply add to a list
	accessions = []
	for line in open(filename, 'r'):
		line = line.strip()
		accessions.append(line)
	return accessions

def parse_annotations(filename):
	annotations = {}
	# traverse through each annotation; key => header, val => its annotation
	for line in open(filename, 'r'):
		line = line.strip().split('\t')
		header, an_annot = line[0], '\t'.join(line[1:])
		annotations[header] = an_annot
	return annotations

def merge(accessions, annotations):
	print len(accessions), 'accessions'
	print len(annotations), 'annotations'
	
	for counter, each_accn in enumerate(accessions):
		if each_accn in annotations:
			print str(counter)+'\t'+each_accn+'\t' + annotations[each_accn]
		else:
			print each_accn, 'not in'

if __name__ == '__main__':	
	# provide command-line arguments for passing-in parameters
	parser = OptionParser()
	parser.add_option("-l", "--accessions", help="Line-separated accessions")
	parser.add_option("-a", "--annotations", help="BLAST annotations; from blast scripts in repository")
	(options, args) = parser.parse_args()
	
	# if no parameters are provided, alert the user
	if not options.accessions or not options.annotations:
		print 'All parameters must have a respective value'
	else:
		# otherwise, check to see if filenames are valid
		if os.path.isfile(options.accessions) and os.path.isfile(options.annotations):
			accessions = parse_accessions(options.accessions)
			annotations = parse_annotations(options.annotations)
			merge(accessions, annotations)
		
		# ... or else, alert the user that the provided filenames are invalid
		else:
			print 'Please make sure annotations and accessions are valid files'