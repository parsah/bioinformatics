"""
Given a user-provided BED file and a set of MAF files, all MAF alignments for each BED entry are identified
and corresponding organisms are extracted.
"""

import sys
from collections import Counter
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor


CUTOFF = 0.7  # organisms must be in at least this much alignment blocks.
HG19 = 'hg19'

def __maf_to_features(line):
    """
    Extracts valuable information given a user-provided
    MAF string.

    :param line: String representing a MAF file line.
    :return: organism, chromosome, and position MAF values.
    """

    line = list(filter(None, line.split(' ')))
    position = int(line[2])
    sequence = line[6]
    organism, chromosome = line[1].split('.')[0: 2]
    return sequence, chromosome, position, organism


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
            value = __maf_to_features(line)
            if line.startswith('s hg'):
                global HG19
                hg19_seq, hg_chromosome, hg_position, HG19 = value
                if hg_chromosome not in contents:
                    contents[hg_chromosome] = {}
                contents[hg_chromosome][hg_position] = {HG19: hg19_seq}
            else:
                organism_seq, org_chromosome, org_position, org = value
                contents[hg_chromosome][hg_position][org] = organism_seq
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
    list_mapped_orgs = []  # list of successfully-mapped organisms.
    for position in maf:
        if (start - diff) < position < (end + diff):
            organism_sequences = list(maf[position].keys())  # get all concordant sequences
            for org_name in organism_sequences:
                org_seq = maf[position][org_name]
                hg_seq = maf[position][HG19]
                perc_identity = sum([1.0 for i in range(len(hg_seq)) if hg_seq[i] ==  org_seq[i]]) / len(hg_seq)
                if perc_identity >= CUTOFF:
                    list_mapped_orgs.append(org_name)
    valid_orgns = dict(Counter(list_mapped_orgs))  # enumerate organisms.
    if len(valid_orgns) > 0:
        max_count = max(valid_orgns.values())  # find most frequenct organism
        valid_orgns = {k: float(valid_orgns[k])/max_count for k in valid_orgns}
        valid_orgns = [org for org in valid_orgns if valid_orgns[org] >= CUTOFF]
    return chromosome, start, end, valid_orgns


def prettify(future):
    """
    Callback function for use during concurrent operations.
    :param future: Futures object returned from map_intervals(...)
    """

    chromosome, start, end, organisms = future.result()
    if len(organisms) > 0:
        for organism in organisms:
            sys.stdout.write(chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + str(organism) + '\n')
            sys.stdout.flush()


if __name__ == '__main__':
    try:
        if len(sys.argv) < 3:
            print('Arg1 => BED file')
            print('Arg2 => MAF file')
        else:
            futures = []
            executor = ThreadPoolExecutor(1)
            maf_data = parse_maf(sys.argv[2])
            for bed_entry in open(sys.argv[1]):
                bed_entry = bed_entry.strip().split('\t')
                bed_chromosome, bed_start, bed_end = bed_entry[0], int(bed_entry[1]), int(bed_entry[2])

                if bed_chromosome in maf_data:
                    #print(map_intervals(bed_chromosome, bed_start, bed_end, maf_data[bed_chromosome]))
                    futures.append(executor.submit(map_intervals, bed_chromosome, bed_start, bed_end, maf_data[bed_chromosome]))
            for future in concurrent.futures.as_completed(futures):
                prettify(future)

    except KeyboardInterrupt:
        print()
    except KeyError as e:
        print(e)
