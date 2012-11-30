""" 
A script which fetches functional annotations for a uniprot ID.
Annotations include GO ontologies, EC accession and uniprot description.
"""

import argparse, urllib2
from xml.etree import ElementTree 

schema = '{http://uniprot.org/uniprot}'
no_hit = '---'
	
# Print headers so columns can be easily identified
def print_headers():
	h = ["Uniprot","EC","Description","Component","Function","Process"]
	print '\t'.join(h)
	
# Parses user-provided uniprot file	
def load_input(fname):
	handle = open(fname)
	accns = []
	for entry in handle:
		entry = entry.strip()
		if len(entry) == 0:
			break
		else:
			accns.append(entry.strip())
	return accns

# A UniprotTree is an encapsulation of the XML file given its web-service
class UniprotTree():
	def __init__(self, accn):
		self.accn = accn # raw uniprot accession
		self.tree = None # represents a uniprot XML encapsulation
		self._parse()
		
	# parse the XML referencing the uniprot accession
	def _parse(self):
		xml_fname = 'http://www.uniprot.org/uniprot/' + self.accn + '.xml'
		data = urllib2.urlopen(xml_fname).read()
		self.tree = ElementTree.XML(data)

	# Get the EC description
	def get_ec_desc(self):
		return [n.text for n in self.tree.getiterator(schema+'fullName')][0]
	
	# Get all EC accessions
	def get_ecs(self):
		ec_tree = self.tree.getiterator(schema+'dbReference')
		ec = '|'.join([node.attrib['id'] 
						for node in ec_tree if 'EC' in node.attrib['type']])
		return no_hit if ec == '' else ec
	
	# Return the uniprot accession
	def get_accn(self):
		return self.accn
	
	def string_go(self, c, f, p):
		# if no GO ontology has an annotation, give it a default blank
		if len(c) == 0:
			c.append(no_hit)
		if len(f) == 0:
			f.append(no_hit)
		if len(p) == 0:
			p.append(no_hit)
		return ' --- '.join(c).strip() + '\t'+\
			 ' --- '.join(f).strip()+'\t' +' --- '.join(p).strip()
	
	# Get all GO ontologies
	def get_go(self):
		comp, func, proc = [], [], [] # lists of GO ontologies
		for node in self.tree.getiterator():
			str_curr_tag = node.tag.replace('{http://uniprot.org/uniprot}', '')
			if str_curr_tag == 'dbReference' and node.attrib['type'] == 'GO':
				children = node.getchildren()
				for child in children:
					value = child.attrib['value']
					if 'F:' in value:
						func.append(value[2:]) # add 'F:'
					if 'P:' in value:
						proc.append(value[2:]) # add 'P:'
					if 'C:' in value:
						comp.append(value[2:]) # add 'C:'
		return self.string_go(c=comp, f=func, p=proc)
		
# Run analysis				
def run_analysis(accn):
	try:
		up = UniprotTree(accn)
		print up.get_accn()+'\t'+up.get_ecs()+'\t'+up.get_ec_desc()+'\t'+up.get_go()
	except urllib2.HTTPError:
		# print empty spaces to represent no-hits; 3x at the end for each GO
		print accn +'\t'+ no_hit +'\t'+ no_hit +('\t' + no_hit) * 3

if __name__ == '__main__':
	try:
		p = argparse.ArgumentParser(description='Uniprot to GO-annotations script')
		p.add_argument('-in', metavar='FILE', required=True, # input file
					help='File of uniprot accessions; line separated')
		args = vars(p.parse_args())
		accns = load_input(fname=args['in']) # parse input
		print_headers() # first line of get_output
		
		for accn in accns: # run analysis and ultimately save to file
			run_analysis(accn)
		
	except KeyboardInterrupt:
		pass	