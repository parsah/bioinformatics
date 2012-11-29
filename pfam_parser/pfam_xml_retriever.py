
''' 
Functions to parse a list of PFAM IDs and retrieve their respective PFAM
annotation. Such activity is dependent on the PFAM REST web-service.
'''

import urllib2
import os, sys
from optparse import OptionParser

#import easy to use xml parser called minidom:
from xml.dom.minidom import parseString
#all these imports are standard on most modern python implementations

def retrieve_annotation(pfam_id):
	annotation = ''
	try:
		# retrieve the xml file for a specific PFAM ID
		handle = urllib2.urlopen('http://pfam.sanger.ac.uk/family/'+pfam_id+'?output=xml')
		data = handle.read()
		handle.close()
		dom = parseString(data)
		
		# each 'description' tag has a CDATA. Traverse through these nodes, trim
		# whitespace and what is left is the actual annotation
		name = dom.getElementsByTagName('description')
		annotation = " ".join(t.nodeValue for t in name[0].childNodes).strip()
		return annotation
	except Exception:
		print 'error with', pfam_id,'; ID not existent or references another ID'
	finally:
		return annotation

def parse_inputfile(filename):
	for counter, line in enumerate(open(filename)):
		line = line.strip().split('\t')
		header, pfam_entry = line[0], line[2]
		# some additional information accompanies the PFAM ID; trim this
		index_pfam = pfam_entry.find('PF')
		pfam_entry = pfam_entry[index_pfam:]
		# next, remove the PFAM ID revision number denoted by a period
		index_period = pfam_entry.find('.')
		pfam_entry = pfam_entry[0: index_period]
		
		# some accessions may not have PFAM IDs, therefore filter such
		annot = ''
		if len(pfam_entry) > 0:
			annot = retrieve_annotation(pfam_entry)
		else:
			# both dashes represent no PFAM accn + PFAM annot respectively
			pfam_entry = '-'
			annot = '-\t-'
		
		# once complete, stdout results
		sys.stdout.write(str(counter) + '\t' + header + '\t' + pfam_entry + '\t' + annot+'\n')
		sys.stdout.flush()
	

if __name__ == '__main__':
	# create command-line parser
	parser = OptionParser()
	parser.add_option("-i", "--infile", help="Tab-delim file w/PFAM IDs in first column")
	(options, args) = parser.parse_args()
	
	filename = options.infile
	# if an invalid file is passed, cease execution
	if not os.path.isfile(filename):
		print 'Invalid user file provided. Please try again.'
	else:
		parse_inputfile(filename)
