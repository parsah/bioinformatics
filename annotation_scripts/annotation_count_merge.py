import argparse
import csv

def parse_csv(fname, colnum):
    """ 
    Parse the user-provided CSV file and index each entry by transcript ID.
    """
    handle = open(fname)
    entries = {}
    for line in csv.reader(handle):
        transcript = line[colnum]
        if transcript not in entries:
            entries[transcript] = '\t'.join(line)
    return entries

def parse_counts(fname):
    """ 
    Parse user-provided transcript counts which are produced using
    unix sort | uniq -c. From such a function, each transcript will
    have a corresponding read-count that will ultimately be mapped
    onto the user-provided annotations.
    """
    counts = {}
    for i in open(fname):
        i = i.strip().split(' ')
        a_count, transcript = int(i[0]), i[-1]
        counts[transcript] = a_count
    return counts

def merge(annotations, counts):
    """ 
    Merges parsed counts with user-provided annotations. As a result, there
    will be produced a single output file with each annotation entry
    referencing a single read-count.
    """
    for transcript_id in annotations:
        if transcript_id not in counts:
            # if the transcript does not have any reads, set count to 0.
            print(transcript_id+"\t" + str(0)+"\t"+\
                 annotations[transcript_id])
        else:
            # if the transcript has reads, get its count.
            print(transcript_id+"\t" + str(counts[transcript_id])+"\t"+\
                 annotations[transcript_id]) 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-csv', required=True, metavar='',
                        help='CSV annotation file [na]')
    parser.add_argument('-transcript', metavar='', default=1, type=int,
                        help='Feature column # used in counts [1]')
    parser.add_argument('-in', metavar='', required=True,
                        help='Transcript counts from unix uniq -c [na]')
    args = vars(parser.parse_args())
    annotations = parse_csv(fname=args['csv'], colnum=args['transcript'])
    counts = parse_counts(fname=args['in'])
    merge(annotations, counts)
    