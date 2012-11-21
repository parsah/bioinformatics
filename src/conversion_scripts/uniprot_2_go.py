""" 
Given a set of uniprot accessions, fetch GO annotations for each accession.
Since a uniprot accession can have multiple GO annotations, metadata can be
supplied so as to aide in mapping GO results.
"""

import argparse, urllib2
from xml.etree import ElementTree 

# Parses user-provided uniprot file	
def load_input(filename):
	handle = open(filename)
	accns = []
	for entry in handle:
		entry = entry.strip().split('\t')
		if len(entry) == 1: # if no metadata is provided, provide empty string
			entry.append('')
		if len(entry) > 0:
			accns.append(entry)
	return accns

# Cleanses output so that GO component, function and process are aligned
def to_out(uniprot, metadata, component, function, process):
	max_go = max([len(component), len(function), len(process)])
	comp = ['-'] * max_go
	func = ['-'] * max_go
	proc = ['-'] * max_go
	comp[0:len(component)] = component
	func[0:len(function)] = function
	proc[0:len(process)] = process
	
	for i in range(max_go): # print out all the Go annotations
		print uniprot+'\t'+comp[i]+'\t'+ func[i]+'\t'+ proc[i]+'\t'+metadata

def run_analysis(accn):
	try:
		uniprot, metadata = accn # only 2x columns allowed in input file
		xml_fname = 'http://www.uniprot.org/uniprot/' + uniprot + '.xml'
		data = urllib2.urlopen(xml_fname).read()
		tree = ElementTree.XML(data)
		func, proc, comp = set(), set(), set() # stores Go ontologies
		for node in tree.getiterator():
			str_curr_tag = node.tag.replace('{http://uniprot.org/uniprot}', '')
			# Each node = {'type': 'GO', 'id': 'GO:0005351', 'key': '17'}
			# Each node has children nodes which point to a set of ontologies, eg:	
			if str_curr_tag == 'dbReference' and node.attrib['type'] == 'GO':
				# {'type': 'term', 'value': 'F:sugar:hydrogen symporter activity'}
				# {'type': 'evidence', 'value': 'IEA:InterPro'}
				# The above node only points to a function, therefore process and
				# component are not yet annotated.
				
				# Extracting out the 'value' will help us map it an ontology.
				children = node.getchildren()
				for child in children:
					value = child.attrib['value']
					if 'F:' in value:
						func.add(value[2:]) # skip 'F:'
					if 'P:' in value:
						proc.add(value[2:]) # skip 'P:'
					if 'C:' in value:
						comp.add(value[2:]) # skip 'C:'
		# there may be multiple GO components, functions or processes. In
		# such a case, print each on a different line
		to_out(uniprot, metadata, list(comp), list(func), list(proc))
	except urllib2.HTTPError: # if an invalid URL is provided, print blanks
		print uniprot + '\t' + metadata + '\t---\t---\t---'

if __name__ == '__main__':
	p = argparse.ArgumentParser(description='Uniprot to GO-annotations script')
	p.add_argument('-in', metavar='FILE', required=True, 
				help='File of uniprot accessions; line separated')
	args = vars(p.parse_args())
	accns = load_input(filename=args['in'])
	for accn in accns: # get all GOs per uniprot accession
		run_analysis(accn)