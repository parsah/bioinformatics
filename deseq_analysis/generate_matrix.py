'''
Given counts mapping to a BAM / SAM file which are derived using
sort | uniq -c, this script parses such counts and creates a
matrix-like output so it can be parsed by DESeq.
'''

import argparse

def _merger(query, baseline):
    ''' 
    Parse query and baseline count-files and only get a union-set of
    accessions shared across all counts regardless of experiment.
    @param query: Query counts input file.
    @param baseline: Baseline (control) counts input file.
    @return: dictionary of counts per accession across all files.
    '''
    union = {} # key => accession in counts file, value => read count.
    all_files = [query, baseline] # references input count-files.
    for infile  in all_files:
        for line in open(infile):
            line = line.strip().split()
            accn = line[1]
            union[accn] = ['0']*len(all_files) # default 0-count per accession.
    return union # return map containing accession and its respective counts.

def to_matrix(fnames, baseline, shared_keys):
    ''' 
    Builds a counts-table which can be fed into DESeq. 
    '''
    headers = []
    all_files = [fnames, baseline]
    # process all filename-datasets.
    for filenum, infile in enumerate(all_files):
        headers.append(infile)
        for line in open(infile):
            line = line.strip().split()
            num, accn = line[0], line[1] # count and accession, respectively.
            shared_keys[accn][filenum] = num # save accession-specific count.
    headers.insert(0, '') # insert blank entry for R data-frame parsing
    headers_string = ','.join(headers)
    
    print(headers_string) # print headers
    for i in shared_keys: # for each key (accessions), print its counts
        out_str = i + ',' + ','.join(shared_keys[i])
        print(out_str)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-query', metavar="", nargs='+', required=True,
                        help='List of query counts [na]')
    parser.add_argument('-baseline', metavar="", nargs='+', 
                        required=True, help='List of baseline counts [na]')
    args = vars(parser.parse_args())
    union = _merger(query=args['query'], baseline=args['baseline'])
    to_matrix(shared_keys=union, fnames=args['query'], 
                       baseline=args['baseline'])
    