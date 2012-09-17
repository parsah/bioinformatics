
""" 
Given a set of uniprot accessions, fetch GO annotations for each accession.
This logic is implemented using 'multiprocessing' and 'Elementree' 
"""

import urllib2, multiprocessing
from xml.etree import ElementTree 

schema = '{http://uniprot.org/uniprot}'
non_accn = '---'

LIST_ALL_RESULTS = []
	
def load_in_file(filename):
	handle = open(filename, 'r')
	list_accns = []
	for accn in handle:
		accn = accn.strip()
		if accn != non_accn and len(accn) > 0:
			list_accns.append(accn)
	print '#/entries:', format(len(list_accns), ',d')
	return list_accns

def run_analysis(counter, accn):
	result = ''
	try:
		xml_fname = 'http://www.uniprot.org/uniprot/' + accn + '.xml'
		data = urllib2.urlopen(xml_fname).read()
		tree = ElementTree.XML(data)
		set_function, set_process, set_component = set(), set(), set()
		for node in tree.getiterator():
			str_curr_tag = node.tag.replace('{http://uniprot.org/uniprot}', '')
			# Each node = {'type': 'GO', 'id': 'GO:0005351', 'key': '17'}
			# Each node has children nodes which point to a set of ontologies, eg:	
			if str_curr_tag == 'dbReference' and node.attrib['type'] == 'GO':
	#			print node.attrib
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
						set_function.add(value[2:]) # skip 'F:'
					if 'P:' in value:
						#print '\t', value
						set_process.add(value[2:]) # skip 'P:'
					if 'C:' in value:
						#print '\t', value
						set_component.add(value[2:]) # skip 'C:'
		component = ' ### '.join(set_component)
		function = ' ### '.join(set_function)
		process = ' ### '.join(set_process)
		result = accn +'\t'+component+'\t'+function+'\t'+process
	except ElementTree.ParseError:
		result = accn+ '\t*** NOT FOUND ***'
	return str(counter)+'\t'+result

def cb(returned_value):
	LIST_ALL_RESULTS.append(returned_value)
	print returned_value

if __name__ == '__main__':
	uniprot_accns = load_in_file('uniprot_accns.txt')
	out_writer = open('go_ontologies_from_uniprot.txt', 'w')
	
	pool = multiprocessing.Pool(processes=1)
	for counter, accn in enumerate(uniprot_accns):
		pool.apply_async(run_analysis, (counter, accn,), callback=cb)
		
	pool.close()
	pool.join()
	
	for i in LIST_ALL_RESULTS:
		out_writer.write(i+'\n')
		out_writer.flush()
	out_writer.close()
	print 'done'
	