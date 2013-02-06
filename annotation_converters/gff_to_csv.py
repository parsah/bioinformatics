
import argparse

# Parse GFF3 file
def parse_gff3(fname):
	counter = 0
	for line in open(fname):
		if not line.startswith('##'):
			line = line.strip().split('\t')
			chrom, src, start = line[0], line[2], line[3] # extract regions
			end, strand, annot = line[4], line[6], line[8]
			if src == 'mRNA': # only focus on mRNA entries
				annot = annot.split(';')
				gene = annot[1].replace('Name=', '') # phytozome delimiter
				transcript = annot[-1].replace('Parent=', '') # phytozome delimiter
				if strand == '+':
					strand = '1'
				elif strand == '-':
					strand = '-1'
				else:
					print('Invalid strand character:', strand)
					break
				counter += 1 # unique index
				print(str(counter)+','+transcript +',' + gene +',' +\
					chrom +','+start+','+ end +','+strand)
		
if __name__ == '__main__':
	p = argparse.ArgumentParser()
	p.add_argument('-in', metavar='FILE', help='GFF3 file [na]', required=True)
	args = vars(p.parse_args())
	parse_gff3(fname=args['in'])		