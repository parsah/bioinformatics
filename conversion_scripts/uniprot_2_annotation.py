""" 
A useful script for retrieving annotations given a list of such accessions.
Several 'modes' exist so that annotations only of this type are obtained;
useful in cases when you seek to only update annotations of a particular type.
"""

import argparse, urllib2
from xml.etree import ElementTree 

modes = ['ec', 'go', 'desc'] # modes of analysis
	
# This class performs the bulk of all annotation efforts
class AnnotationFactory():
	def __init__(self, fname, mode):
		self.schema = '{http://uniprot.org/uniprot}' # for xml parsing
		self.no_hit = '---' # displayed when no annotation is possible
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
		route = dict(zip(modes, [self.get_ecs, self.get_gos, self.get_desc]))
		route[self.mode]() # execute the desired function
	
	# get all ECs given the list of uniprot accessions
	def get_ecs(self):
		print 'UNIPROT\tEC' # print header
		for accn in self.uniprot_accns:
			try:
				u = UniprotXML(accn)
				ec_tree = u.tree.getiterator(self.schema+'dbReference')
				ec = '|'.join([node.attrib['id']
							for node in ec_tree if 'EC' in node.attrib['type']])
				ec = self.no_hit if ec == '' else ec
				print accn + '\t' + ec
			except urllib2.HTTPError:
				print accn + '\t' + self.no_hit
	
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
			except urllib2.HTTPError:
				print acc + '\t' + self.no_hit + '\t' + self.no_hit + '\t' + self.no_hit
	
	# Helper-function to convert GO ontologies to string objects
	def _string_go(self, c, f, p):
		# if no GO ontology has an annotation, give it a default blank
		div = ' '+str(self.no_hit)+' ' # separates like-GOs from one another
		if len(c) == 0:
			c.append(self.no_hit)
		if len(f) == 0:
			f.append(self.no_hit)
		if len(p) == 0:
			p.append(self.no_hit)
		return div.join(c).strip(), div.join(f).strip(), div.join(p).strip()
		
	# get the uniprot description for the accession 
	def get_desc(self):
		print 'UNIPROT\tUNIPROT_DESC'
		for accn in self.uniprot_accns:
			try:
				u = UniprotXML(accn) # uniprot xml encapsulation
				desc = [n.text for n in u.tree.getiterator(self.schema+'fullName')]
				print accn + '\t' + desc[0] # get the description
			except urllib2.HTTPError:
				print accn + '\t' + self.no_hit
			
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
		p = argparse.ArgumentParser(description='Uniprot to GO-annotations script')
		p.add_argument('-in', metavar='FILE', required=True, # input file
					help='File of uniprot accessions; line separated')
		p.add_argument('-mode', metavar='LIST', # annotation type
					choices=modes, required=True,
					help='Annotation retrieval; choices [na]')
		args = vars(p.parse_args())
		
		factory = AnnotationFactory(fname=args['in'], mode=args['mode'])
		factory.query() # get uniprot accessions given the desired mode
		
	except KeyboardInterrupt:
		pass	