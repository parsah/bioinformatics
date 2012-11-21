""" 
A script which fetches functional annotations for a uniprot ID.
Annotations include GO ontologies, EC accession and uniprot description.
"""

import argparse, urllib2
from xml.etree import ElementTree 

schema = '{http://uniprot.org/uniprot}'

# Wraps uniprot annotation results
class OutputResult():
	def __init__(self, accn, desc, metadata):
		self.uniprot_accn = accn # is set at runtime
		self.metadata = metadata # references input metadata
		self.uniprot_desc = desc # uniprot description
		self.ec = None # references EC accession(s)
		self.go_proc = set() # an accession may have multiple GO entries
		self.go_comp = set() # references GO components
		self.go_func = set() # references GO functions
		
	# Add GO component
	def add_component(self, item):
		self.go_comp.add(item)
	
	# Add GO function
	def add_function(self, item):
		self.go_func.add(item)
	
	# Add GO process
	def add_process(self, item):
		self.go_proc.add(item)
	
	# Set EC accessions mapping to the uniprot accession
	def set_ec_accns(self, ecs):
		self.ec = ecs
		
	# Cleanses output so that GO component, function and process are aligned
	def to_out(self):
		max_go = max([len(self.go_comp), len(self.go_func), len(self.go_proc)])
		comp = ['-'] * max_go
		func = ['-'] * max_go
		proc = ['-'] * max_go
		comp[0:len(self.go_comp)] = self.go_comp
		func[0:len(self.go_func)] = self.go_func
		proc[0:len(self.go_proc)] = self.go_proc
		
		for i in range(max_go): # print out all the Go annotations
			print self.uniprot_accn+'\t'+self.ec+'\t'+\
				self.uniprot_desc+'\t'+comp[i]+'\t'+ func[i]+'\t'+\
				proc[i]+'\t'+self.metadata

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

def run_analysis(accn):
	try:
		uniprot, metadata = accn # only 2x columns allowed in input file
		xml_fname = 'http://www.uniprot.org/uniprot/' + uniprot + '.xml'
		data = urllib2.urlopen(xml_fname).read()
		tree = ElementTree.XML(data)
		ec_desc = [n.text for n in tree.getiterator(schema+'fullName')][0]
		ec_tree = tree.getiterator(schema+'dbReference')
		# get the number of EC maps corresponding to this GO ID
		ecs = '|'.join([node.attrib['id'] 
					for node in ec_tree if 'EC' in node.attrib['type']])
		outresult = OutputResult(accn=uniprot, desc=ec_desc, metadata=metadata)
		outresult.set_ec_accns(ecs) # sets EC accessions
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
						outresult.add_function(value[2:]) # skip 'F:'
					if 'P:' in value:
						outresult.add_process(value[2:]) # skip 'P:'
					if 'C:' in value:
						outresult.add_component(value[2:]) # skip 'C:'
		# there may be multiple GO components, functions or processes. In
		# such a case, print each on a different line
		outresult.to_out()
		
	except urllib2.HTTPError: # if an invalid URL is provided, print blanks
		print uniprot + '\t' + metadata

if __name__ == '__main__':
	try:
		p = argparse.ArgumentParser(description='Uniprot to GO-annotations script')
		p.add_argument('-in', metavar='FILE', required=True, 
					help='File of uniprot accessions; line separated')
		args = vars(p.parse_args())
		accns = load_input(filename=args['in'])
		for accn in accns: # get all GOs per uniprot accession
			run_analysis(accn)
	except KeyboardInterrupt:
		pass	