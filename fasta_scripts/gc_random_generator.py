'''
A very simple but helpful script for generating variable-length DNA sequences
with a predefined GC percentage.
'''

import argparse
import random

def validate_args(args):
    num_seqs, seq_len, gc_perc = args['n'], args['l'], args['gc']
    if seq_len <= 0:
        raise IOError('Positive sequence length (-l) required [error].')
    if num_seqs <= 0:
        raise IOError('Positive number of sequence (-n) required [error].')
    if gc_perc > 100 or gc_perc < 0:
        raise IOError('GC percentage (-gc) must be between 0 .. 100 [error].')

def generator(args):
    # extract-out command line arguments
    num_seqs, seq_len, gc_perc = args['n'], args['l'], args['gc']/100.0
    for i in range(num_seqs):
        # begin by making an AT repeat-sequence of the user-desired length
        seq_list = list('AT'*seq_len)[:seq_len]
        num_gc_reqd = int(len(seq_list) * gc_perc) # number of GCs required
        # create list of unique indices
        gc_positions = list(range(0, len(seq_list)))
        random.shuffle(gc_positions) # jumble their positions and add G or C
        gc_positions = gc_positions[: num_gc_reqd]
        for position in gc_positions:
            g_or_c = random.choice(['G', 'C'])
            seq_list[position] = g_or_c # insert either a G or C
        seq_str = ''.join(seq_list)
        print('>sequence_'+str(i+1) + '\n' + seq_str) # display as FASTA format

if __name__ == '__main__':
    desc = 'Sequence generator given GC percentage'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-gc', type=int, required=False, default=50,
                        metavar='INT', help='GC percentage [50]')
    parser.add_argument('-l', type=int, required=False, default=10,
                        metavar='INT', help='Sequence length [10]')
    parser.add_argument('-n', type=int, required=False, default=1,
                        metavar='INT', help='# / sequences to generate [1]')
    args = vars(parser.parse_args())
    try:
        validate_args(args)
        generator(args) # begin the sequence generator
    except IOError as e:
        print(e)
        