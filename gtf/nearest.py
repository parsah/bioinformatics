'''
Given both a BED file and a corresponding GTF file, this trivial script maps
each BED entry to its nearest GTF. Such mapping is performed in a
strand--independent manner. This script also enables filtering of keywords
which may not be necessarily warranted to be part of the nearest derivation,
e.g. pseudogene.
'''

import argparse
import pandas
import sys


def parse_gtf(f, kw):
    '''
    Parse a user-provided GTF file.
    @param f: GTF filename.
    @param kw: keywords to omit if found in the source column.
    @return: hash referencing chromosome names and its entries.
    '''

    data = {}  # key => chromosome name, values => chromosome entries.
    for i in open(f):
        i = i.strip().split('\t')
        chrom, source = i[0: 2]

        # if chromosome doesn't begin with chr
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom

        # test whether the source keyword is in the omitted keyword list
        if source in kw:
            if chrom not in data:  # if chromosome is not in the dataset
                data[chrom] = []
            data[chrom].append(i)
    return data


def parse_bed(f):
    '''
    Parses a user-provided BED file and only holds the first three critical
    columns (chromosome, start, and end).

    @param f: BED filename.
    @return: pandas data-frame referencing the BED file.
    '''

    df = pandas.read_table(f, header=None, sep='\t')
    idx = df[0] + '.' + df[1].map(str) + '.' + df[2].map(str)
    df = df.set_index(idx)  # set index as chrom.start.end
    df = df[[0, 1, 2]]
    df.columns = ['chrom', 'start', 'end']
    return df  # only focus on the first three important columnns


def compute_nearest(bed, gtf, is_verbose=False):
    '''
    @param bed: parsed BED file from parse_bed.
    @param gtf: parsed GTF file from parse_gtf.
    '''
    for row in bed.index:
        row = bed.loc[row]
        bed_chrom, bed_start, bed_end = row['chrom'], row['start'], row['end']
        top_dist = sys.maxsize  # save the nearest distance
        top_feat = None  # keep-track of the best feature

        gtfs = gtf[bed_chrom]  # get GTF entries for the chromosome
        for a_gtf in gtfs:
            gtf_start, gtf_end = int(a_gtf[3]), int(a_gtf[4])

            # compare distance between BED end and GTF start.
            if abs(bed_end - gtf_start) <= top_dist:
                top_dist = abs(bed_end - gtf_start)
                top_feat = a_gtf

            if abs(bed_end - gtf_end) <= top_dist:
                top_dist = abs(bed_end - gtf_end)
                top_feat = a_gtf

            if abs(bed_start - gtf_start) <= top_dist:
                top_dist = abs(bed_start - gtf_start)
                top_feat = a_gtf

            if abs(bed_start - gtf_end) <= top_dist:
                top_dist = abs(bed_start - gtf_end)
                top_feat = a_gtf
        if is_verbose:
            line = [bed_chrom, bed_start, bed_end, top_dist] + top_feat
            line = [str(i) for i in line]
        else:
            line = [str(i) for i in [bed_chrom, bed_start, bed_end, top_dist]]
        print(','.join(line))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', metavar='BED', required=True,
                        help='BED file [req]')
    parser.add_argument('-gtf', metavar='GTF', required=True,
                        help='GTF annotations file [req]')
    parser.add_argument('-w', metavar='LIST', nargs='+', default=[],
                        help='Words to keep from GTF annotations [none]')
    parser.add_argument('--verbose', default=False, action='store_true',
                        help='Verbose output [false]')
    args = vars(parser.parse_args())
    try:
        gtf = parse_gtf(f=args['gtf'], kw=args['w'])
        bed = parse_bed(f=args['bed'])
        compute_nearest(bed, gtf, is_verbose=args['verbose'])
    except KeyboardInterrupt:
        print()
