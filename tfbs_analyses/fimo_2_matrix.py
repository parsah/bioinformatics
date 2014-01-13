'''
Creates a count-matrix given FIMO output using its --text format.
'''

import argparse
import pandas
from collections import OrderedDict

def build_skeleton(cont, query):
    pwms = set() # ordering is not important
    all_seqs = []
    
    # loop through both files to derive row-count
    for dataset in [cont, query]:
        seqs = set()
        for linenum, line in enumerate(open(dataset)):
            if linenum != 0:
                line = line.strip().split('\t')
                pwm, seq = line[0: 2]
                pwms.add(pwm)
                seqs.add(seq)
        all_seqs.extend(seqs)
    
    # wrap counts in a DataFrame
    m = OrderedDict({'Sequence': all_seqs})
    m.update({pwm: [0] * len(all_seqs) for pwm in pwms})
    m.update({'Target': ['None'] * len(all_seqs)})
    df = pandas.DataFrame(m, index=all_seqs)
    return df # return dataframe capturing the matrix

def populate(df, cont, query):
    for i, dataset in enumerate([cont, query]):
        for linenum, line in enumerate(open(dataset)):
            if linenum != 0:
                line = line.strip().split('\t')
                pwm, seq = line[0: 2]
                df[pwm][seq] += 1 # increment sequence-PWM count
                df['Target'][seq] = i
    return df # contains actual counts

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-control', required=True, metavar='FILE',
                        help='Control output from FIMO --text mode [req]')
    parser.add_argument('-query', required=True, metavar='FILE',
                        help='Query output from FIMO --text mode [req]')
    parser.add_argument('-csv', metavar='FILE', default='./out.csv',
                        help='Output file [./out.csv]')
    args = vars(parser.parse_args())
    cont_fname, query_fname = args['control'], args['query']
    df = build_skeleton(cont_fname, query_fname)
    df = populate(df, cont_fname, query_fname)
    df.to_csv(args['csv'])
