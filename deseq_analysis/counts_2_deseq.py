import argparse

'''
Given counts mapping to a BAM / SAM file which are derived using
sort | uniq -c, this script parses such counts and creates a
table-like output so it can be parsed by DESeq.
'''

def get_union(query, baseline):
    ''' 
    Parse query and baseline count-files and only get a union-set of
    accessions shared across all counts regardless of experiment.
    '''
    union = {}
    all_files = []
    all_files.extend(query)
    all_files.extend(baseline)
    for infile  in all_files:
        for line in open(infile):
            line = line.strip().split()
            accn = line[1]
            union[accn] = ['1']*len(all_files) # default 1-count per accession.
    return union

def to_timecourse_counts(fnames, baseline, shared_keys):
    ''' 
    Builds a counts-table which can be fed into DESeq. 
    '''
    headers = []
    all_files = []
    all_files.extend(fnames)
    all_files.extend(baseline)
    # process fnames datasets.
    for filenum, infile in enumerate(all_files):
        headers.append(infile)
        for line in open(infile):
            line = line.strip().split()
            num, accn = line[0], line[1]
            shared_keys[accn][filenum] = num

    headers.insert(0, '') # insert blank entry for R data-frame parsing
    
    headers_string = ','.join(headers)
    print(headers_string)
    for i in shared_keys:
        out_str = i + ',' + ','.join(shared_keys[i])
        print(out_str)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-query', metavar="", nargs='+', required=True,
                        help='List of query counts [na]')
    parser.add_argument('-baseline', metavar="", nargs='+', 
                        required=True, help='List of baseline counts [na]')
    args = vars(parser.parse_args())
    union = get_union(query=args['query'], baseline=args['baseline'])
    to_timecourse_counts(shared_keys=union, fnames=args['query'], 
                       baseline=args['baseline'])
    