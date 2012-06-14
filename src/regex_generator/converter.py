
from xml.etree import ElementTree
from collections import defaultdict
from regex_factory import generate_regex

''' 
Takes an XML with the following organization:

<family name="TF family name">
	<description>
	<elements>
		<binding_site seq="sequence" source="database">
			<transcription_factor orgn="org">TF_name</transcription_factor>
		<binding_site>
		...
	</elements>
	</description>
</family>

The above schema is then parsed and organized in such a way that TFBSs with
the same length and TF gene get converted into a regular expression.
'''

# Function to parse an XML file and represent it as a hash of TFBS sizes
def __parse_xml(family_name):
	print family_name,'\n'
	results = defaultdict(list) # for storing results from XML parsing
	tree = ElementTree.parse('../binding_elements.xml')
	for family in tree.getiterator('family'): # only process a specific family
		if family.attrib['name'] == family_name:
			for bind_site in family.getiterator('binding_site'):
				sequence, source = bind_site.attrib['seq'], bind_site.attrib['source']
				tf_gene = bind_site.getchildren()[0] # get TF gene element
				orgn, tf_name = tf_gene.attrib['orgn'], tf_gene.text
				
				# represent each xml child-parent relationship as a hash
				node_info = {'orgn': orgn, 'tf_name': tf_name, 'seq': sequence,
							'source': source}
				# map each result by TFBS length and TF gene name
				results[(len(sequence), tf_name)].append(node_info)
	return results # return all the size-TF gene name mappings	

# Iterates through all mappings and produces unified representations of TFBSs
# should they exist in multiple forms
def __consolidate(mappings):
	for result in mappings:
		set_source, set_orgn = set(), set() # stores copies of source & orgn
		sequences = []
		for each_map in mappings[result]:
			set_source.add(each_map['source'])
			set_orgn.add(each_map['orgn'])
			sequences.append(each_map['seq']) # get seq for each variant
		
		# next, create xml tags for each element
		if len(sequences) < 2: # print-out as-is as no regex is needed
			create_xml_tags(each_map['seq'], each_map['source'], 
						each_map['orgn'], each_map['tf_name'])
			
		else: # feed into regex converter
			sequence = generate_regex(sequences)
			create_xml_tags(sequence, each_map['source'], 
						each_map['orgn'], each_map['tf_name'])
	
# Trivial functions to print the XML header and footer
def print_header(family_name):
	print '<family name=\"' + family_name + '\">'
	print '\t<description/>'
	print '\t<elements>'
	
def print_footer():
	print '\t</elements>'
	print '<family>'

# Create XML element tags to model TFBSs, TF genes, TFBS organism and source
def create_xml_tags(tfbs_sequence, tfbs_source, tf_orgn, tf_name, ):
	print '\t\t<binding_site seq=\"' + tfbs_sequence +'\" source=\"' + tfbs_source+'\">'
	print '\t\t\t<transcription_factor orgn=\"' + tf_orgn + '\">' + tf_name + "</transcription_factor>"
	print '\t\t</binding_site>'

def main():
	family_name = 'bHLH' # the TF family to search for
	mappings = __parse_xml(family_name)
	print_header(family_name)
	__consolidate(mappings)
	print_footer()

if __name__ == '__main__':
	main()