'''
Takes both a query and control output file generated from running tfSearch
and molds such data into a count-matrix that can be used for machine-learning
and classification purposes.
'''

import argparse
import pandas as pd
from collections import OrderedDict

ENTRY_QUALIFIER = '::' # valid tfSearch outputs contain the '::' string

def normalize(m):
    ''' 
    Length-normalize each entry in the count matrix.
    @param m: Count matrix.
    '''
    pass

def get_count(s):
    ''' 
    Extracts the numerical count of a PWM in any given tfSearch output string.
    @param s: tfSearch output string.
    @return: count representing PWM-count in a user-provided tfSearch string.
    '''
    str_found = s[s.rfind(ENTRY_QUALIFIER): ]
    str_found = str_found.replace(ENTRY_QUALIFIER, '').strip()
    num = str_found[0: str_found.find('(')] # right-brace ends the count
    return int(num)

def parse(f):
    ''' 
    Extract all PWMs shared across both query and control input files.
    @param files: List of both control and query input files.
    '''
    d = {} # key => accession, value => PWM and its respective count
    for i in open(f):
        i = i.strip()
        if ENTRY_QUALIFIER in i: # only print line if valid entry
            num = get_count(i)
            i = i.split(' ')
            # idx 2 => PWM ID, idx 3 => PWM name; remove trailing comma
            accn, pwm = i[0], i[2] + ' ' + i[3][:-1]
            if accn not in d:
                d[accn] = {pwm: 0}
            d[accn][pwm] = num
    return d

def unique_pwms(control, query):
    ''' 
    Extracts all PWMs which are shared across control and query parsed files.
    The end result is a collection representing all PWMs which serve as the
    columns of a count-matrix.
    @param control: Parsed control input file.
    @param query: Parsed query input file.
    @return: set of PWMs shared across both control and query input files.
    '''
    union = set()
    for d in [control, query]:
        for accn in d:
            pwms = set(list(d[accn].keys()))
            union.update(pwms)
    return union

def write_matrix(m, f):
    df = pd.DataFrame(m)
    df.to_csv(f)

def build_matrix(control, query, bool_array):
    ''' 
    The matrix is such that you have n rows. Each row is a sequence from both
    the control and query datasets. Each PWM column, j, references a list of 
    counts of length n. Thus, the index [ i , j ] can be used to retrieve the
    PWM-count in a given sequence.
    @param control: Parsed control input file.
    @param query: Parsed query input file.
    @param bool_array: Array referencing which dataset is control or not
    '''
    name_seq, name_target = 'Sequence', 'Target'
    num_rows = len(control) + len(query) # references the total number of sequences    
    m = OrderedDict({ name_seq: [''] * num_rows }) # references strings  
    m.update({ pwm: [0] * num_rows for pwm in  unique_pwms(control, query )}) # PWMs are our matrix columns
    m.update({ name_target: ['None'] * num_rows })

    row_num = 0 # begin row counter
    for i, dataset in enumerate([control, query]):
        for accn in dataset:
            #print(accn, row_num)
            m[name_seq][row_num] = accn # set accession name
            m[name_target][row_num] = int(bool_array[i]) # set target variable
            pwms = dataset[accn] # get PWMs mapping to the accession 
            for pwm in pwms:
                count = pwms[pwm] # get the PWM count as well
                m[pwm][row_num] = count
                #print('\t', accn, pwm, count)
            row_num += 1
    return m

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-control', metavar='FILE', required=True,
                        help='tfSearch output; control file [req]')
    parser.add_argument('-query', metavar='FILE', required=True,
                        help='tfSearch output; query file [req]')
    parser.add_argument('-out', metavar='FILE', required=False,
                        default='./out.csv',
                        help='Output matrix file [out.csv]')
    parser.add_argument('--norm', action='store_true', 
                        help='Length-normalize PWM counts [true]')
    args = vars(parser.parse_args())
    
    try:
        control = parse(f = args['control'])
        query = parse(f = args['query'])
        m = build_matrix(control, query, [False, True])
        if args['norm']:
            pass
        write_matrix(m, args['out']) 
    except KeyboardInterrupt:
        print()
