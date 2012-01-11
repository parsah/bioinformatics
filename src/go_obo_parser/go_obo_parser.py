
"""
Given a Gene Ontology (GO) OBO-XML file, parse this xml file and map its
respective GO annotations with blast txt output which has GO IDs.
As a result, you get GO annotations given a set of GO IDs.
"""


import argparse, re
from lxml import etree

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('-x','-xml', help='OBO-XML file from Gene Ontology')
	parser.add_argument('-b','-blast_out', help='BLAST against GO txt output')
	args = parser.parse_args()
	return args.x, args.b

def parse_obo(go_xml):
	tree = etree.parse(go_xml)
	all_terms = tree.xpath('term') # fetching all GOs which start w/'term'
	print '#/entries found:', format(len(all_terms), ',d')
	print 'parsing entries...',
	
	dict_all_ontologies = {} # all GOs per line stored here
	for each_term in all_terms: # for each ontology...
		dict_ontologies_line = {'id': None, 'name': None, 'namespace': None}
		for child in each_term: # ... and for each id, name and GO, add -> dict
			if child.tag == 'id':
				dict_ontologies_line['id'] = child.text
			if child.tag == 'name':
				dict_ontologies_line['name'] = child.text
			if child.tag == 'namespace':
				dict_ontologies_line['namespace'] = child.text
		# each line of GO is indexed based on the GO ID
		dict_all_ontologies[dict_ontologies_line['id']] = dict_ontologies_line
		dict_ontologies_line = None
	print '[OK]', format(len(dict_all_ontologies),',d'), 'entries parsed.'
	return dict_all_ontologies


def execute(file_obo, file_blast):
	xml_results = parse_obo(file_obo)

if __name__ == '__main__':
	file_obo, file_blast = parse_arguments()
	try:
		execute(file_obo, file_blast)
	except IOError:
		print 'Error: file \'', file_obo,'\' is not a valid file.'		