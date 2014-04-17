'''
Given a BED file and a set of genomic FASTA file, each genomic region in the
BED file is paired with an equal-length genomic segment, matching its GC and
repeat content. Selecting such sequences therefore helps extract sequences
representative of a BED file which match its intrinsic properties; serving as
a baseline or control-set.
'''

import random
import sys
from Bio import SeqIO
from Bio.SeqUtils import GC

# allowable difference between candidate and query GC and N percentage.
DIFF = 0.2


def parse_genomes(d):
    dict_genome = dict()
    for entry in SeqIO.parse(open(d), 'fasta'):
        dict_genome[entry.description] = entry.seq
    return dict_genome


def get_n_perc(seq):
    return (str(seq).count('N') / len(seq)) * 100


def parse_bed(f, g):
    bed = {}
    for bed_entry in open(f):
        key = '.'.join(bed_entry.strip().split('\t')[0:3])
        line = bed_entry.strip().split('\t')
        chrom, start, end = line[0], int(line[1]), int(line[2])
        seq = g[chrom][start: end]
        bed[key] = {'chrom': chrom, 'start': start,
                         'end': end, 'seq': seq, 'gc_perc': round(GC(seq), 2),
                         'n_perc': round(get_n_perc(seq), 2)}
    return bed


def is_match_gc(query_gc, control_gc):
    gc_diff = abs(query_gc - control_gc)
    return gc_diff <= DIFF


def is_match_n(query_n, control_n):
    n_diff = abs(query_n - control_n)
    return n_diff <= DIFF


def run(bed, genome):
    for bed_key in bed:
        bed_value = bed[bed_key]

        # control cannot be from same chromosome as peak
        genome_chrom_names = list(genome.keys())
        genome_chrom_names.remove(bed_value['chrom'])

        # next, pick a random chromosome and use this to get a random sequence
        has_match = False   # flag to help deduce whether sequence is found.
        suitable_seq = None  # represents matching-GC, N, and length control.
        while not has_match:
            cont_chrom_name = random.choice(genome_chrom_names)
            cont_chrom_seq = genome[cont_chrom_name]
            for idx in range(len(cont_chrom_seq)):
                cont_seq = str(cont_chrom_seq[idx:idx + len(bed_value['seq'])])
                cont_gc = GC(cont_seq)
                cont_n = get_n_perc(cont_seq)
                is_match = is_match_gc(bed_value['gc_perc'], cont_gc) and \
                    is_match_n(bed_value['n_perc'], cont_n)
                if is_match:
                    has_match = True
                    suitable_seq = cont_seq  # set the control sequence
                    print('>matched.' + bed_key + '\n' + suitable_seq)
                    break

if __name__ == '__main__':
    try:
        if len(sys.argv) < 2:
            raise OSError('Provide BED file and genome file, respectively.')
        else:
            genome = parse_genomes(d=sys.argv[2])
            bed_file = parse_bed(f=sys.argv[1], g=genome)
            run(bed_file, genome)
    except IOError as e:
        print(e)
    except KeyboardInterrupt:
        print()
