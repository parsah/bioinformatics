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

def output(seqs):
    ''' 
    Prints all the random sequences be-they generated from adhoc or template
    arguments. The resultant sequences are sent to standard-out.
    '''
    for i, seq in enumerate(seqs):
        # shuffle bases so that if a sequence with 0 GC% are not only AT dimers
        seq = list(seq)
        random.shuffle(seq)
        print('>sequence_' + str(i + 1) + '\n' + ''.join(seq))

def parse_fasta(fname):
    ''' 
    Trivial function to parse a user-provided FASTA file
    '''
    try:
        seqs, header = {}, ''
        for line in open(fname):
            line = line.strip()
            if line.startswith('>'):
                header = line[1:] # remove the fasta header
                seqs[header] = ''
            else:
                seqs[header]+=line
        return seqs
    except IOError:
        print('input file:', fname, 'not found [error]')

def template_generator(fname):
    ''' 
    Generates random sequences given the inherent properties of a user-defined
    FASTA file. 
    '''
    entries = parse_fasta(fname)
    seqs = [] # a list of all generated sequences given a template
    for entry in entries:
        seq = str(entries[entry]).upper()
        seq_len = len(seq)
        seq_gc = round(((seq.count('G') + seq.count('C')) / seq_len) * 100)
        seqs.extend(adhoc_generator(n=1, l=seq_len, gc=seq_gc))
    return seqs

def adhoc_generator(n, l, gc):
    ''' 
    Builds random sequences given GC content and length. This trivial function
    is designed for instances when analyzing promoter sequences of a query
    sequence set against a random baseline. If random sequences are desired to
    match GC content of a predefined set of sequences, use the 
    template_generator(...) function instead. 
    '''
    num_seqs, seq_len, gc_perc = n, l, gc/100.0
    seqs = [] # a list of all the generated random sequences
    for _ in range(num_seqs):
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
        seqs.append(seq_str) # save as FASTA
    return seqs

if __name__ == '__main__':
    desc = 'Sequence generator given GC percentage'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-in', required=False,
                        metavar='FILE', help='FASTA file to build sequences around [na]')
    parser.add_argument('-gc', type=int, required=False, default=50,
                        metavar='INT', help='GC percentage [50]')
    parser.add_argument('-l', type=int, required=False, default=10,
                        metavar='INT', help='Sequence length [10]')
    parser.add_argument('-n', type=int, required=False, default=1,
                        metavar='INT', help='# / sequences to generate [1]')
    args = vars(parser.parse_args())
    try:
        validate_args(args)
        if not args['in']: # if no FASTA file is provided, generate sequences
            # generate adhoc sequences and send them to standard-out
            output(adhoc_generator(n=args['n'], l=args['l'], gc=args['gc']))            
        else:
            output(template_generator(fname=args['in'])) # use template
    except IOError as e:
        print(e)
        