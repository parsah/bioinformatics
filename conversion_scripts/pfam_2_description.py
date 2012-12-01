import argparse
from uniprot_2_annotation import no_hit, div
import urllib2
from xml.etree import ElementTree

baseurl = 'http://pfam.sanger.ac.uk/family/' # REST url for accessing IDs
schema = '{http://pfam.sanger.ac.uk/}' # required to extract description

# Many PFAMs can be mapped to an accession, therefore they must
# be each appreciated by separating the string.
def parse_pfam_ids(fname):
	ids = []
	for line in open(fname):
		ids.append(line.strip())
	return ids

# There may be multiple PFAMs delimited by the same symbol used to delimit
# uniprot annotations. Annotate each PFAM given this delimiter produce output.  
def get_description(pfam_id):
	ids = pfam_id.split(no_hit) # split PFAMs given the delimiter
	desc = [] # contains descriptions for each accession
	try:
		for acc in ids: # iterate over each split-PFAM
			acc = acc.strip()
			data = urllib2.urlopen(url=baseurl+acc+'?output=xml').read()
			tree = ElementTree.XML(data)
			tree = tree.getiterator(schema + 'description')
			# append the description per PFAM to the list referencing ids
			desc.append([node.text for node in tree][0].strip())
		print pfam_id + '\t' +div.join(desc)
	except ElementTree.ParseError:
		print pfam_id + '\t' + no_hit # if PFAM is invalid, produce no-hit 

def main():
	p = argparse.ArgumentParser(description='PFAM -> description script')
	p.add_argument('-in', required=True, help='PFAM IDs; can be delimited', 
				metavar='FILE')
	args = vars(p.parse_args())
	entries = parse_pfam_ids(fname= args['in'])
	for pfam_id in entries: # iterate over each PFAM ID and get description
		get_description(pfam_id)

if __name__ == '__main__':
	main()