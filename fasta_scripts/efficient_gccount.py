''' 
A helpful script for counting GC percentage of a user-provided sequence.
Such functionality is very useful in instances where you wish to select 
sequences sharing similar GC or repeat percentages.
'''
import argparse
from Bio import SeqIO

def count_gc(sequence, win):
    prior_gc = 0 # define initial GC count
    prior_seq = '' # define empty prior sequence
    for i in range(len(sequence)):
        if i == 0: # for the first sliding window, compute GC content
            prior_seq = sequence[i: win]
            prior_gc = prior_seq.count('G') + prior_seq.count('C')
            print(prior_seq+'\t'+str(round((prior_gc / win) * 100, 2)))
        else:
            curr_seq = sequence[i: i+win]
            if len(curr_seq) != win: # do not pull pre-mature windows
                break
            
            # if the prior sequence's head is G or C, decrease count
            # since this base is not part of current window.
            if prior_seq[0] == 'G' or prior_seq[0] == 'C':
                prior_gc -= 1
            
            # if the current sequence's tail (new addition) contains
            # a G or C, increase count since this base is within window.
            if curr_seq[-1] == 'G' or curr_seq[-1] == 'C':
                prior_gc += 1
            
            # update reference
            prior_seq = curr_seq
            
            print(curr_seq+'\t'+str(round((prior_gc / win) * 100, 2)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', help='Input FASTA file [na]', 
                        required=True, metavar='FILE')
    parser.add_argument('-window', help='Sliding window size [5]',
                        type=int, metavar='INT', default=5)
    args = vars(parser.parse_args())
    records = SeqIO.parse(args['in'], 'fasta')
    for record in records:
        print('sequence\tGC-perc\tN-perc')
        count_gc(sequence = str(record.seq), win = args['window'])
        print()
