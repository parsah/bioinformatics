
import argparse, urllib2
from xml.etree import ElementTree # library to parse xml file

base_url = 'http://www.uniprot.org/uniprot/' # base url for accessing uniprot
schema = '{http://uniprot.org/uniprot}'
invalid_ids = set(['','-']) # set of invalid uniprot accessions 

# parse user-provided input file
def _parse_uniprot_file(fname):
	handle, accns = open(fname, 'r'), [] # input filename and list of input IDs
	for line in handle:
		line = line.strip()
		if len(line) == 0: # if line with no accession is encountered, break
			break
		else: # ... otherwise, continue
			accns.append(line)
	return accns

# Helper-function which downloads the XML file for a respective uniprot ID
def uniprot_to_xml(uniprot_id):
	output_fname = 'tempfile.xml'
	uniprot_url = base_url + uniprot_id + '.xml' # obide-by uniprot web-service
	url_contents = urllib2.urlopen(url=uniprot_url).read()
	handle = open(output_fname, 'w')
	handle.write(url_contents+'\n') # save downloaded XML locally
	handle.flush()
	return output_fname # the output XML filename

# iteratively retrieve EC accessions for each uniprot ID (if existent)
def get_ecs(uniprot_ids):
	for counter, uniprot_id in enumerate(uniprot_ids):
		try:
			out_fname = uniprot_to_xml(uniprot_id) # get filename XML is stored in
			# only get 'dbReference' nodes: nodes which have annotations
			tree = ElementTree.parse(out_fname).getiterator(schema+'dbReference')
			# get the number of EC maps corresponding to this GO ID
			num_ecs = '|'.join([node.attrib['id'] 
					for node in tree if 'EC' in node.attrib['type']])
			print str(counter) +'\t' + uniprot_id+'\t'+num_ecs 
		except urllib2.HTTPError:
			print str(counter) +'\t' + uniprot_id+'\t-' # null-output
			
	print 'done'

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
	
	