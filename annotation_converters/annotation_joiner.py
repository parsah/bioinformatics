'''
A useful script which serves much like RDBMS inner-joins. This script merges
a baseline set of entries with a more informative set of entries that have
detailed functional annotations. Each entry is the baseline is line-delimited
whereas each entry in the more informative annotation is tab-delimited.
'''

import argparse

def parse_baseline(fname):
    ''' 
    Parse a user-provided baseline file whereby each line is line delimited
    and contains one or more occurances of each accession.
    @param fname: Input filename.
    '''
    contents = []
    for line in open(fname):
        line = line.strip().split('\t')
        if len(line) != 1:
            raise IOError('Entries in the baseline file must be 1x column.')
        contents.append(line[0]) # add the only entry in that row
    return contents

def parse_annotations(fname):
    ''' 
    Parse a user-provided annotations file whereby each column is tab-delimited
    and contains an arbitrary number of columns. The very first column
    however must equal that in the baseline dataset.
    @param fname: Input filename.
    '''
    contents = {}
    for line in open(fname):
        line = line.strip().split('\t')
        key, value = line[0], line[1:]
        contents[key] = value
    return contents

def merge(baseline, annotations):
    ''' 
    Merges baseline accessions with that in the annotations dataset.
    @param baseline: Baseline entries.
    @param annotations: Annotated entries referencing the baseline.
    '''
    for b in baseline:
        if b in annotations:
            print(b + '\t' + '\t'.join(annotations[b]))
        else:
            print(b + '\t-')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-baseline', metavar='FILE', required=True,
                        help='Line-delimited collection of entries [na]')
    parser.add_argument('-annotations', metavar='FILE', required=True,
                        help='Tab-delimited annotation(s) per entry [na]')
    args = vars(parser.parse_args())
    baseline = parse_baseline(fname=args['baseline'])
    annotations = parse_annotations(fname=args['annotations'])
    merge(baseline, annotations) # prints-out streaming data to the console.
    