import argparse, urllib2

url_base = 'http://enzyme.expasy.org/EC/' # base url for retrieving info.
header = 'DE' # all descriptions begin with the DE key

# Trivial parsing of file populated with ECs; line-separated
def parse_ecs(infile):
	entries = []
	for line in open(infile):
		entries.append(line.strip())
	return entries

# Traverse the record and only get the line which has a description (DE)
def get_description(ec):
	try:
		url = url_base + ec + '.txt'
		str_url = str(urllib2.urlopen(url).read()).splitlines()
		description = None # no description is yet extracted
		for line in str_url: # iterate over each line in the record
			if line.startswith(header):
				description = line.replace(header, '').strip() # cleanse string
		print ec + '\t' + description
	except urllib2.HTTPError: # if invalid EC (i.e. URL) print null-character
		print ec + '\t-'
		
def main():
	p = argparse.ArgumentParser(description='EC -> description script')
	p.add_argument('-in', required=True, help='EC IDs', metavar='FILE')
	args = vars(p.parse_args())
	entries = parse_ecs(infile= args['in'])
	for ec_accn in entries: # iterate over each EC and get its description
		get_description(ec=ec_accn)

if __name__ == '__main__':
	main()