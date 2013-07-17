'''
This script enables the ability to extract sequences from a candidate-set which 
match a user-defined sequence set. Unlike generation of random sequences, this
method enables selection of sequences which have similar GC content. Assuming
such candidates are biological sequences, this script facilitates selection
of true sequences so as to match the query.
@author: Parsa Hosseini
'''

import argparse
import os

def validate_params(args):
    arg_d, arg_query, arg_candidate = args['d'], args['i'], args['c']
    if arg_d < 0:
        raise IOError('argument \'d\' must be positive [error]')
    if not os.path.exists(arg_query):
        raise IOError('arg \'i\' must be valid query FASTA file [error]')
    if not os.path.exists(arg_candidate):
        raise IOError('arg \'c\' must be valid candidate FASTA file [error]')

def parse_fasta(fname):
    ''' 
    Trivial function to parse a user-provided FASTA file.
    @param fname: Input FASTA filename.
    '''
   
    seqs, header = {}, ''
    for line in open(fname):
        line = line.strip()
        if line.startswith('>'):
            header = line[1:] # remove the fasta header
            seqs[header] = {'seq': '', 'used': False}
        else:
            seqs[header]['seq'] += line
    return seqs

def run_script(args):
    ''' 
    Begins the execution of the script, assuming all user-provided arguments
    are valid and sound.
    @param args: Arguments as provided by the user at runtime.
    '''
    arg_d, arg_query, arg_candidate = args['d'], args['i'], args['c']
    print('Parsing both query and candidate sequences ...')
    seqs_query = parse_fasta(fname = arg_query)
    seqs_candidate = parse_fasta(fname = arg_candidate)
    out_handle = open(args['o'], 'w') # save to output filename
    
    for query_header in seqs_query: # iterate over each query
        query_seq = seqs_query[query_header]
        query_len = len(query_seq)
        query_gc_pct = gc_percentage(query_seq)
        is_matched = False # flag whether the query has a corresponding baseline
        
        for candidate_header in seqs_candidate: # ... and over each candidate
            candidate_seq = seqs_candidate[candidate_header]['seq']
            candidate_len = len(candidate_seq)
            candidate_gc_pct = gc_percentage(candidate_seq)
            
            # test whether the candidate sequence is current paired or not
            if not seqs_candidate[candidate_header]['used']:
                if abs(query_gc_pct - candidate_gc_pct) <= arg_d:
                    seqs_candidate[candidate_header]['used'] = True
                    is_matched = True # whether the query has a match
                    out_handle.write('>' + candidate_header + '\n' + candidate_seq + '\n')
                    out_handle.flush()
                    break
        if not is_matched: # if a query lacks a pair, throw error.
            raise IOError(query_header, 'lacks baseline. Try increasing \'d\'')
        
    out_handle.close()
    print('Analysis complete [ok]')

def gc_percentage(s):
    ''' 
    Trivial function to derive sequence GC content.
    @param s: DNA sequence.
    '''
    seq = str(s)
    gc_count = seq.count('G') + seq.count('C') # count number of Gs and Cs
    return int((gc_count / float(len(seq)))*100) # cast to integer

if __name__ == '__main__':
    desc = 'Baseline sequence selector driven by GC content'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', required=True,
                        metavar='FILE', help='Query FASTA file [na]')
    parser.add_argument('-c', required=True,
                        metavar='FILE', help='Candidate FASTA file [na]')
    parser.add_argument('-d', type=int, required=False, default=5,
                        metavar='INT', help='Query-candidate GC content delta [5]')
    parser.add_argument('-o', required=False, default='baseline.fasta',
                        metavar='FILE', help='Output baseline FASTA [./baseline.fasta]')
    args = vars(parser.parse_args())
    try:
        validate_params(args) # test if arguments are valid
        run_script(args) # now, begin running the script
    except KeyboardInterrupt:
        print('Analysis terminated by user.')
    except IOError as e:
        print(e)
        