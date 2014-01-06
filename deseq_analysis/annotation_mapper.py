'''
Takes a DESeq output file and maps its differentially-expressed features onto
user-provided functional annotations.
'''
import argparse
import csv

def parse_deseq_output(fname, fc_cutoff):
    ''' 
    Parses a DESeq output file by only parsing-out the feature name, log-2
    fold-change, p-value as well as adjusted p-value.
    '''
    contents = {}
    for linenum, line in enumerate(csv.reader(open(fname))):
        if linenum != 0:
            feature, log2FC, pval, adjPVal = line[1], line[6], line[7], line[8]
            if abs(float(log2FC)) >= fc_cutoff:
                contents[feature] = log2FC+'\t'+pval+'\t'+adjPVal
    return contents
    
def parse_annotations(fname):
    '''
    Parses user-provided annotations by indexing annotations solely by their
    feature name.
    '''
    contents = {}
    for linenum, line in enumerate(csv.reader(open(fname))):
        if linenum != 0:
            feature = line[1]
            contents[feature] = '\t'.join(line)
    return contents

def merge(degs, annotations):
    ''' 
    Merges differentially--expressed genes with corresponding annotations.
    '''
    for deg in degs:
        if deg in annotations:
            print(deg + "\t" + degs[deg] + '\t' + annotations[deg])
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', required=True, metavar='FILE',
                        help='DESeq output file [na]')
    parser.add_argument('-fc', metavar='FLOAT', type=float, default=2.0,
                        help='Fold-change cutoff (abs value) [2]')
    parser.add_argument('-a', required=True, metavar='FILE',
                        help='Function annotations CSV [na]')
    args = vars(parser.parse_args())
    contents = parse_deseq_output(fname=args['in'],fc_cutoff=args['fc'])
    annotations = parse_annotations(fname=args['a'])
    merge(degs=contents, annotations=annotations)