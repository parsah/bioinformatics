
import argparse, os, urllib2
from xml.etree import ElementTree # library to parse xml file

base_url = 'http://www.uniprot.org/uniprot/' # base url for accessing uniprot
schema = '{http://uniprot.org/uniprot}'
temp_results = 'tempfile.xml' # xml containing information for the accession
invalid_ids = set(['','-']) # set of invalid uniprot accessions 

# parse user-provided input file
def _parse_uniprot_file(fname):
	handle, accns = open(fname, 'r'), [] # input filename and list of input IDs
	for line in handle:
		# column 1 => uniprot accession, column 2 => metadata, eg. Phytozome ID
		line = line.strip().split('\t')
		if len(line) == 1: # if uniprot accessions are provided, metadata is null
			line.append('')
		if len(line) == 0: # if line with no accession is encountered, break
			break
		else: # ... otherwise, continue
			accns.append(line)
	return accns

# Helper-function which downloads the XML file for a respective uniprot ID
def uniprot_to_xml(uniprot_id):
	uniprot_url = base_url + uniprot_id + '.xml' # obide-by uniprot web-service
	url_contents = urllib2.urlopen(url=uniprot_url).read()
	handle = open(temp_results, 'w')
	handle.write(url_contents+'\n') # save downloaded XML locally
	handle.flush()
	return temp_results # the output XML filename

# iteratively retrieve EC accessions for each uniprot ID (if existent)
def get_ecs(uniprot_ids):
	for i, entry in enumerate(uniprot_ids):
		try:
			uniprot, metadata = entry
			out_fname = uniprot_to_xml(uniprot) # get filename XML is stored in
			# get the name of the EC also
			tree = ElementTree.parse(out_fname)
			ec_name = [n.text for n in tree.getiterator(schema+'fullName')][0]
			ec_refs = tree.getiterator(schema+'dbReference')
			# only get 'dbReference' nodes: nodes which have annotations
			# get the number of EC maps corresponding to this GO ID
			num_ecs = '|'.join([node.attrib['id'] 
					for node in ec_refs if 'EC' in node.attrib['type']])
			print str(i) +'\t'+ uniprot+'\t'+num_ecs+'\t'+ ec_name + '\t' + metadata
		except urllib2.HTTPError:
			print str(i) +'\t' + uniprot+'\t-' # null-output
	os.remove(temp_results) # when all is complete, delete temp file

if __name__ == '__main__':
	# Create a parser so file of uniprot IDs can have their EC IDs retrieved
	parser = argparse.ArgumentParser(description='Uniprot -> EC ID script')
	parser.add_argument('-in', metavar='FILE', required=True,
					help='Line-delimited file of uniprot IDs')
	args = vars(parser.parse_args())
	try:
		uniprot_ids = _parse_uniprot_file(args['in']) # parse user-provided file
		get_ecs(uniprot_ids) # run the script
	except KeyboardInterrupt:
		print 'Keyboard-interrupt [Analysis cancelled]'