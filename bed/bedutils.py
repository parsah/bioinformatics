import pandas


def parse_bed(f):
    '''
    Parses a user-provided BED file and only holds the first three critical
    columns (chrom, start, and end).

    @param f: BED filename.
    @return: pandas data-frame referencing the BED file.
    '''

    df = pandas.read_table(f, header=None, sep='\s')
    idx = df[0] + '.' + df[1].map(str) + '.' + df[2].map(str)
    df = df.set_index(idx)  # set index as chrom.start.end
    return df


def parse_multibed(f):
    '''
    Parses a multiIntersectBed file; similar to a traditional BED file,
    however a multiIntersectBed file is not only generated from bedtools but
    also can contain header information of potential value when looking
    at BED entry intersections.

    @param f: BED filename.
    @return: pandas data-frame referencing the BED file.
    '''

    df = pandas.read_table(f, sep='\s')
    idx = df['chrom'] + '.' + df['start'].map(str) + '.' + df['end'].map(str)
    df = df.set_index(idx)  # set index as chrom.start.end
    return df
