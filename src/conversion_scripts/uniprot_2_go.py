""" 
Given a set of uniprot accessions, fetch GO annotations for each accession.
This logic is implemented using 'multiprocessing' and 'Elementree' 
"""

import argparse, urllib2
from xml.etree import ElementTree 

schema = '{http://uniprot.org/uniprot}'
no_annot = '---' # represents uniprot accessions lacking a specific GO function 

# Parses user-provided uniprot file	
def load_in_file(filename):
	handle = open(filename)
	list_accns = []
	for accn in handle:
		accn = accn.strip().split('\t')
		if len(accn) == 1: # if no metadata is provided, provide empty string
			accn.append('')
		if len(accn) > 0:
			list_accns.append(accn)
	return list_accns

# Cleanses output so that GO component, function and process are aligned
def to_out(uniprot, metadata, component, function, process):
	max_go = max([len(component), len(function), len(process)])
	comp = ['-'] * max_go
	func = ['-'] * max_go
	proc = ['-'] * max_go
	comp[0:len(component)] = component
	func[0:len(function)] = function
	proc[0:len(process)] = process
	
	for i in range(max_go):
		print uniprot+'\t'+comp[i]+'\t'+ func[i]+'\t'+ proc[i]+'\t'+metadata

def run_analysis(i, accn):
	try:
		uniprot, metadata = accn
		xml_fname = 'http://www.uniprot.org/uniprot/' + uniprot + '.xml'
		data = urllib2.urlopen(xml_fname).read()
		tree = ElementTree.XML(data)
		func, proc, comp = set(), set(), set()
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
						#print '\t', value
						func.add(value[2:]) # skip 'F:'
					if 'P:' in value:
						#print '\t', value
						proc.add(value[2:]) # skip 'P:'
					if 'C:' in value:
						#print '\t', value
						comp.add(value[2:]) # skip 'C:'
		# there may be multiple GO components, functions or processes. In
		# such a case, print each on a different line
		to_out(uniprot, metadata, list(comp), list(func), list(proc))
	except ElementTree.ParseError:
		print accn+ '\t*** NOT FOUND ***'

if __name__ == '__main__':
	p = argparse.ArgumentParser(description='Uniprot to GO-annotations script')
	p.add_argument('-in', metavar='FILE', required=True, 
				help='File of uniprot accessions; line separated')
	args = vars(p.parse_args())
	uniprot_accns = load_in_file(filename=args['in'])
	for i, accn in enumerate(uniprot_accns):
		run_analysis(i, accn)
			