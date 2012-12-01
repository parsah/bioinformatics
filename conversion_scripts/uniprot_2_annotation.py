""" 
A useful script for retrieving annotations given a list of such accessions.
Several 'modes' exist so that annotations only of this type are obtained;
useful in cases when you seek to only update annotations of a particular type.
"""

import argparse, urllib2
from xml.etree import ElementTree 

modes = ['ec', 'go', 'desc', 'pfam'] # modes of analysis
no_hit = '---' # displayed when no annotation is possible
div = ' '+str(no_hit)+' ' # divides multiple accessions

# This class performs the bulk of all annotation efforts
class AnnotationFactory():
	def __init__(self, fname, mode):
		self.schema = '{http://uniprot.org/uniprot}' # for xml parsing
		self.fname = fname # references input uniprot accessions
		self.mode = mode # whether to annotate ECs, GO, or others
		self.uniprot_accns = [] # contains uniprot accessions once parsed
		self._parse()
		
	# parse user-provided uniprot accessions
	def _parse(self):
		handle = open(self.fname)
		for entry in handle:
			entry = entry.strip()
			if len(entry) == 0:
				break
			else:
				self.uniprot_accns.append(entry.strip())
	
	# Query accessions against uniprot and get accessions based on 'mode'
	def query(self):
		# modes are [ec, go, desc], therefore map to their respective function
		route = dict(zip(modes, 
				[self.get_ecs, self.get_gos, self.get_desc, self.get_pfam]))
		route[self.mode]() # execute the desired function
	
	# Get PFAM IDs given the uniprot accession-list
	def get_pfam(self):
		print 'UNIPROT\tPFAM' # print header
		for accn in self.uniprot_accns:
			try:
				u = UniprotXML(accn)
				pf_tree = u.tree.getiterator(self.schema+'dbReference')
				pfam = div.join([node.attrib['id'] for node in pf_tree 
							if 'Pfam' in node.attrib['type']])
				pfam = no_hit if pfam == '' else pfam
				print accn + '\t' + pfam
			except Exception:
				print accn + '\t' + no_hit
		
	# get all ECs given the list of uniprot accessions
	def get_ecs(self):
		print 'UNIPROT\tEC' # print header
		for accn in self.uniprot_accns:
			try:
				u = UniprotXML(accn)
				ec_tree = u.tree.getiterator(self.schema+'dbReference')
				ec = div.join([node.attrib['id']
							for node in ec_tree if 'EC' in node.attrib['type']])
				ec = no_hit if ec == '' else ec
				print accn + '\t' + ec
			except Exception:
				print accn + '\t' + no_hit
	
	# get all GO ontologies given the list of uniprot accessions
	def get_gos(self):
		print 'UNIPROT\tGO_COMP\tGO_FUNC\tGO_PROC'
		for acc in self.uniprot_accns:
			try:
				u = UniprotXML(acc)
				comp, func, proc = [], [], [] # lists of GO ontologies
				for node in u.tree.getiterator(self.schema + 'dbReference'):
					if node.attrib['type'] == 'GO':
						children = node.getchildren()
						for child in children:
							value = child.attrib['value']
							if 'F:' in value:
								func.append(value[2:]) # add 'F:'
							if 'P:' in value:
								proc.append(value[2:]) # add 'P:'
							if 'C:' in value:
								comp.append(value[2:]) # add 'C:'
				c, f, p = self._string_go(c=comp, f=func, p=proc)
				print acc + '\t' + c + '\t' + f + '\t' + p
			except Exception:
				print acc + '\t' + no_hit + '\t' + no_hit + '\t' + no_hit
	
	# Helper-function to convert GO ontologies to string objects
	def _string_go(self, c, f, p):
		# if no GO ontology has an annotation, give it a default blank
		if len(c) == 0:
			c.append(no_hit)
		if len(f) == 0:
			f.append(no_hit)
		if len(p) == 0:
			p.append(no_hit)
		return div.join(c).strip(), div.join(f).strip(), div.join(p).strip()
		
	# get the uniprot description for the accession 
	def get_desc(self):
		print 'UNIPROT\tUNIPROT_DESC'
		for accn in self.uniprot_accns:
			try:
				u = UniprotXML(accn) # uniprot xml encapsulation
				desc = [n.text for n in u.tree.getiterator(self.schema+'fullName')]
				print accn + '\t' + desc[0] # get the description
			except Exception:
				print accn + '\t' + no_hit
			
# A UniprotXML is an encapsulation of the XML file given its web-service
class UniprotXML():
	def __init__(self, accn):
		self.accn = accn # raw uniprot accession
		self.tree = None # represents a uniprot XML encapsulation
		self._parse()
		
	# _parse the XML referencing the uniprot accession
	def _parse(self):
		xml_fname = 'http://www.uniprot.org/uniprot/' + self.accn + '.xml'
		data = urllib2.urlopen(xml_fname).read()
		self.tree = ElementTree.XML(data)

if __name__ == '__main__':
	try:
		desc = """Script to retrieve uniprot annotations. Several modes are 
			used to facilitate retrieval of annotations given uniprot 
			accessions. Accepted modes include: [ec, go, desc, pfam].""" 
		p = argparse.ArgumentParser(description=desc)
		p.add_argument('-in', metavar='FILE', required=True, # input file
					help='File of uniprot accessions; line separated')
		p.add_argument('-mode', metavar='MODE', # annotation type
					choices=modes, required=True,
					help='Annotation retrieval mode')
		args = vars(p.parse_args())
		
		factory = AnnotationFactory(fname=args['in'], mode=args['mode'])
		factory.query() # get uniprot accessions given the desired mode
		
	except KeyboardInterrupt:
		pass	