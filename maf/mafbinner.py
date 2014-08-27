"""
Given a user-provided BED file and a set of MAF files, all MAF alignments for each BED entry are identified
and corresponding organisms are extracted.
"""

import sys

def __maf_to_features(line):
    line = filter(None, line.split(' '))
    position = int(line[2])
    organism, chromosome = line[1].split('.')[0: 2]
    return (organism, chromosome, position)

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
                if 'scaffold' not in line:
                    organism, org_chromosome, org_position = __maf_to_features(line)
                    contents[hg_chromosome][hg_position].add(organism)
    return contents

def map_intervals(bed, maf):
    for bed_entry in open(bed):
        bed_entry = bed_entry.strip().split('\t')
        bed_chromosome, bed_start, bed_end = bed_entry[0], int(bed_entry[1]), int(bed_entry[2])
        if bed_chromosome in maf:
            positions = maf[bed_chromosome]
            organisms = set()
            for position in positions:
                if bed_start-(bed_end-bed_start) < position < bed_end+(bed_end-bed_start):
                    organisms.update(maf[bed_chromosome][position])
            if organisms:
                for organism in organisms:
                    print str(bed_chromosome) + '\t' + str(bed_start) + '\t' + str(bed_end) + '\t' + organism


if __name__ == '__main__':
    try:
        if len(sys.argv) < 3:
            print('Arg1 => BED file')
            print('Arg2 => MAF file')
        else:
            maf_data = parse_maf(sys.argv[2])
            map_intervals(sys.argv[1], maf_data)
    except KeyboardInterrupt:
        print()
    except KeyError as e:
        print(e)