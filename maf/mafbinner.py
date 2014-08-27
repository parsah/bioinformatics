"""
Given a user-provided BED file and a set of MAF files, all MAF alignments for each BED entry are identified
and corresponding organisms are extracted.
"""

import sys
from concurrent.futures import ThreadPoolExecutor

def __maf_to_features(line):
    line = list(filter(None, line.split(' ')))
    position = int(line[2])
    organism, chromosome = line[1].split('.')[0: 2]
    return organism, chromosome, position


def parse_maf(maf):
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
                # if 'scaffold' not in line:
                organism, org_chromosome, org_position = __maf_to_features(line)
                contents[hg_chromosome][hg_position].add(organism)
    return contents


def map_intervals(chromosome, start, end, maf):
        organisms = set()
        diff = end - start
        for position in maf:
            if (start - diff) < position < (end + diff):
                organisms.update(maf[position])
        return chromosome, start, end, organisms


def cb(obj):
    chromosome, start, end, organisms = obj.result()
    if len(organisms) > 0:
        for organism in organisms:
            print(chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + organism)

if __name__ == '__main__':
    try:
        if len(sys.argv) < 3:
            print('Arg1 => BED file')
            print('Arg2 => MAF file')
        else:
            executor = ThreadPoolExecutor(8)
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