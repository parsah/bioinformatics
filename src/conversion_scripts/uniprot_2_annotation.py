""" 
A script which fetches functional annotations for a uniprot ID.
Annotations include GO ontologies, EC accession and uniprot description.
"""

import argparse, multiprocessing, urllib2
from xml.etree import ElementTree 

schema = '{http://uniprot.org/uniprot}'

# Wraps uniprot annotation results
class OutputResult():
	def __init__(self, accn, desc, metadata):
		self.uniprot = accn # is set at runtime
		self.metadata = metadata # references input metadata
		self.desc = desc # uniprot description
		self.ec = None # references EC accession(s)
		self.go_proc = [] # an accession may have multiple GO entries
		self.go_comp = [] # references GO components
		self.go_func = [] # references GO functions
		
	# Add GO component
	def add_component(self, item):
		self.go_comp.append(item)
	
	# Add GO function
	def add_function(self, item):
		self.go_func.append(item)
	
	# Add GO process
	def add_process(self, item):
		self.go_proc.append(item)
	
	# Set EC accessions mapping to the uniprot accession
	def set_ec_accns(self, ecs):
		self.ec = ecs
	
	# Print headers so columns can be easily identified
	@staticmethod
	def print_headers():
		h = ["Uniprot","EC","Description","Component","Function","Process","Metadata"]
		print '\t'.join(h)
	
	# Function to convert all GO ontologies to printable strings
	def stringify_annotations(self):
		# if no GO ontology has an annotation, give it a default blank
		if len(self.go_comp) == 0:
			self.go_comp.append('---')
		if len(self.go_func) == 0:
			self.go_func.append('---')
		if len(self.go_proc) == 0:
			self.go_proc.append('---')
			
		return ' --- '.join(self.go_comp).strip() + '\t'+\
			 ' --- '.join(self.go_func).strip()+'\t' +' --- '.join(self.go_proc).strip()
	
	# Helper-function to get_output findings to standard-out	
	def get_output(self):		
		return self.uniprot+'\t'+self.ec+'\t'+ self.desc+'\t'+\
			self.stringify_annotations() + '\t'+self.metadata
			
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
		out = OutputResult(accn=uniprot, desc=ec_desc, metadata=metadata)
		if len(ecs) == 0:
			out.set_ec_accns('---') # sets EC accessions
		else:
			out.set_ec_accns(ecs) # sets EC accessions
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
						out.add_function(value[2:]) # skip 'F:'
					if 'P:' in value:
						out.add_process(value[2:]) # skip 'P:'
					if 'C:' in value:
						out.add_component(value[2:]) # skip 'C:'
		# Output results to standard-out
		return out.get_output()
		
	except urllib2.HTTPError: # if an invalid URL is provided, print blanks
		return uniprot + '\t---\t---\t---\t---\t---\t' + metadata

# Callback function given a 
def cb(retval):
	print(retval)

if __name__ == '__main__':
	try:
		p = argparse.ArgumentParser(description='Uniprot to GO-annotations script')
		p.add_argument('-in', metavar='FILE', required=True, # input file
					help='File of uniprot accessions; line separated')
		p.add_argument('-n', metavar='INT', default=2, type=int, # workers
					help='Worker processes [2]')
		args = vars(p.parse_args())
		accns = load_input(filename=args['in']) # parse input
		OutputResult.print_headers() # first line of get_output
		pool = multiprocessing.Pool(processes=args['n'])
		for accn in accns: # run analysis and ultimately save to file
			pool.apply_async(func=run_analysis, args=(accn,), callback=cb)
		pool.close() # when complete, close pool
		pool.join()
		pool.terminate()
		
	except KeyboardInterrupt:
		pass	