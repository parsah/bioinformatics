"""
Given a user-provided BED file and a set of MAF files, all MAF alignments for each BED entry are identified
and corresponding organisms are extracted.
"""

import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor


CUTOFF = 0.75  # organisms must be in at least this much alignment blocks.


def __maf_to_features(line):
    """
    Extracts valuable information given a user-provided
    MAF string.

    :param line: String representing a MAF file line.
    :return: organism, chromosome, and position MAF values.
    """

    line = list(filter(None, line.split(' ')))
    position = int(line[2])
    organism, chromosome = line[1].split('.')[0: 2]
    return organism, chromosome, position


def parse_maf(maf):
    """
    Parses a user-provided MAF file. Such a file is parsed such that
    the chromosome ID is a dictionary key. Values per key are themselves
    dictionaries whereby each index (bp) referencing all concordant
    organisms.

    :param maf: MAF file.
    :return: Dictionary referencing the parsed MAF file.
    """

    hg_chromosome = ''
    hg_position = 0
    contents = {}  # key => chromosome, value => {index: [list of organisms]}
    for line in open(maf):
        line = line.strip()
        if line.startswith('s'):
            if line.startswith('s hg'):
                hg19, hg_chromosome, hg_position = __maf_to_features(line)
                if hg_chromosome not in contents:
                    contents[hg_chromosome] = {}
                contents[hg_chromosome][hg_position] = set([hg19])  # each element is found in hg19 by default
            else:
                organism, org_chromosome, org_position = __maf_to_features(line)
                contents[hg_chromosome][hg_position].add(organism)
    return contents


def map_intervals(chromosome, start, end, maf):
    """
    Given a BED chromosome, start, and end, this triplet is
    queried to see whether a MAF entry maps between this
    BED start and END interval.

    :param chromosome: BED chromosome string.
    :param start: BED start position.
    :param end: BED end position.
    :param maf: Parsed MAF object; see parse_maf(...) function.
    :return:  BED chromosome, start, end, and all mapped organisms.
    """

    diff = end - start
    organisms = []
    for position in maf:
        if (start - diff) < position < (end + diff):
            organisms.extend(maf[position])
    valid_orgns = dict(Counter(organisms))  # enumerate organisms.
    if len(valid_orgns) > 0:
        max_count = max(valid_orgns.values())  # find most frequenct organism
        valid_orgns = {k: float(valid_orgns[k])/max_count for k in valid_orgns}
        valid_orgns = [org for org in valid_orgns if valid_orgns[org] >= CUTOFF]
    return chromosome, start, end, valid_orgns


def cb(obj):
    """
    Callback function for use during concurrent operations.
    :param obj: Object returned from map_intervals(...)
    """

    chromosome, start, end, organisms = obj.result()
    if len(organisms) > 0:
        for organism in organisms:
            print(chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + str(organism))


if __name__ == '__main__':
    try:
        if len(sys.argv) < 3:
            print('Arg1 => BED file')
            print('Arg2 => MAF file')
        else:
            executor = ThreadPoolExecutor(2)
            maf_data = parse_maf(sys.argv[2])
            for bed_entry in open(sys.argv[1]):
                bed_entry = bed_entry.strip().split('\t')
                bed_chromosome, bed_start, bed_end = bed_entry[0], int(bed_entry[1]), int(bed_entry[2])
                if bed_chromosome in maf_data:
                    future = executor.submit(map_intervals, bed_chromosome, bed_start, bed_end, maf_data[bed_chromosome])
                    future.add_done_callback(cb)
            executor.shutdown()
    except KeyboardInterrupt:
        print()
    except KeyError as e:
        print(e)
