"""
A trivial but helpful script to derive the average distance between
all user-provided BED entries.
"""

import argparse
import pandas

INVALID='_'  # chromosomes must not having this character in their name.

def parse_bed(f):
    df = pandas.read_table(f, sep='\t', header=None)
    df = df[[0, 3, 4]]
    df.columns = ['Chromosome', 'Start', 'End']
    return df

def compute_distance(df):
    chromosomes = df['Chromosome'].unique()
    distances = []
    for c in chromosomes:
        if INVALID not in c:
            chromosome_data = df[df['Chromosome'] == c].sort(['Start'])
            chromosome_data.index = range(0, chromosome_data.shape[0])

            for i in chromosome_data.index:
                try:
                    if i != 0:
                        distance = abs(chromosome_data.loc[i+1, 'Start'] - \
                                   chromosome_data.loc[i, 'End'])
                        distances.append(distance)
                except KeyError:
                    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parameters necessary for combinatorial BED analysis
    parser.add_argument('-bed', metavar='FILE', required=True,
                        help='BED file [reqd].')
    args = vars(parser.parse_args())
    compute_distance(parse_bed(f=args['bed']))
