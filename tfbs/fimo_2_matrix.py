'''
Creates a count-matrix given FIMO output using its --text format.
'''

import argparse
import pandas
from pwm_info_content import MEMEPWMParser
from collections import OrderedDict

COLNAME_SEQ = 'Sequence'
COLNAME_TARGET = 'Target'

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
    m = OrderedDict({COLNAME_SEQ: all_seqs})
    m.update({pwm: [0] * len(all_seqs) for pwm in pwms})
    m.update({COLNAME_TARGET: ['None'] * len(all_seqs)})
    df = pandas.DataFrame(m, index=all_seqs)
    return df # return dataframe capturing the matrix

def write(df, csv):
    for i, seq in enumerate(df[COLNAME_SEQ]):
        seq = df[COLNAME_SEQ][i] 
        # replace unique delimiters so both control and query sequences
        # share possible sub-sequences., i.e. match.chr1.111.553
        # is related to chr1.111.553, and so on.
        df[COLNAME_SEQ][i] = seq.replace(':', '.').replace('-', '.')
    df.to_csv(csv)

def populate(df, cont, query):
    for i, dataset in enumerate([cont, query]):
        for linenum, line in enumerate(open(dataset)):
            if linenum != 0:
                line = line.strip().split('\t')
                pwm, seq = line[0: 2]
                df[pwm][seq] += 1 # increment sequence-PWM count
                df[COLNAME_TARGET][seq] = i
    return df # contains actual counts

# class FIMOMatrixBuilder():
#     # builds the matrix and outputs a CSV
#     def __init__(self):
#         pass

class FIMOFilter():
    def __init__(self):
        self.pwm_info = None # references PWM information content
        self.resultset = None # FIMO result-set
        
    def set_resultset(self, rs):
        assert isinstance(rs, FIMOResultSet)
        self.resultset = rs
    
    def set_meme_pwms(self, pwms):
        self.pwm_info = pwms    
    
    def filter(self):
        dataset = self.resultset.to_sequence_counts() # PWMs for respective FIMO file.
        for seq_a in dataset: # iterate over each sequence
            for seq_b in dataset: # iteration is pairwise
                motifs_a = dataset[seq_a]
                motifs_b = dataset[seq_b]
                set_motifs = set(motifs_a.keys()) # set of PWMs A and B
                set_motifs.update(motifs_b.keys())
                print(seq_a, motifs_a, '\t=>' , seq_b, motifs_b, '->', set_motifs)
                
                for pwm_a in set_motifs:
                    print(pwm_a, motifs_a[pwm_a])
                    for pwm_b in set_motifs:
                        if pwm_a != pwm_b:
                            # get intersection of hit-indices 
                            intersect = set(motifs_a[pwm_a]).intersection(motifs_b[pwm_b]) # overlaps
                            set_pwm_a = motifs_a[pwm_a] # indices of A in sequence
                            set_pwm_b = motifs_b[pwm_b] # indices of B in sequence
                            frac_a, frac_b = 0.0, 0.0 # derive PWM A and B ratios
                            if len(intersect) > 0: # to avoid zero-division
                                frac_a = len(intersect) / len(set_pwm_a)
                                frac_b = len(intersect) / len(set_pwm_b)
                            print(seq_a, seq_b, '\tA:', pwm_a, motifs_a[pwm_a],'B:', pwm_b, motifs_b[pwm_b],
                                  '-> intersect:', intersect, 'idx A:', set_pwm_a, 'idx B:', set_pwm_b, 'fA:', frac_a, 'fB', frac_b)

class FIMOResultSet():
    def __init__(self):
        self.mapping = {}
        self.pwms = set()
    
    def parse(self, f):
        for num, line in enumerate(open(f)):
            line = line.strip().split()
            if num != 0:
                motif_id, seq, hit = line[0: 3]
                if seq not in self.mapping:
                    self.mapping[seq] ={}
                if motif_id not in self.mapping[seq]:
                    self.mapping[seq][motif_id] = set()
                self.mapping[seq][motif_id].add(int(hit))
                self.pwms.add(motif_id)
    
    @staticmethod
    def fill_empty(query, control):
        assert isinstance(query, FIMOResultSet)
        assert isinstance(control, FIMOResultSet)
        shared_pwms = set(control.get_pwms())
        shared_pwms.update(query.get_pwms())
        
        for dataset in [query, control]:
            for seq in dataset.mapping:
                for pwm in shared_pwms:
                    if pwm not in dataset.mapping[seq]:
                        dataset.mapping[seq][pwm] = set()
    
    def to_pwm_counts(self):
        counts = {}
        for seq in self.mapping:
            for motif in self.mapping[seq]:
                if motif not in counts:
                    counts[motif] = set()
                # update map with indices of where PWM hits.
                counts[motif].update(self.mapping[seq][motif])
        return counts
        
    def to_sequence_counts(self):
        return self.mapping # by-default, sequences reference their PWM
    
    def get_pwms(self):
        return list(self.pwms) # return only the PWMs mapped to this file.
     
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-control', required=True, metavar='FILE',
                        help='FIMO --text mode output; control [req]')
    parser.add_argument('-query', required=True, metavar='FILE',
                        help='FIMO --text mode output; query [req]')
    parser.add_argument('-m', required=False, metavar='FILE',
                        help='MEME PWM file; removes low-information PWMs [na]')
    parser.add_argument('-csv', metavar='FILE', default='./out.csv',
                        help='Output file [./out.csv]')
    args = vars(parser.parse_args())
    cont_fname, query_fname = args['control'], args['query']
    
    control = FIMOResultSet() # parse the query and control file.
    query = FIMOResultSet() # parse the query and control file.
    control.parse(f = args['control'])
    query.parse(f = args['query'])
    FIMOResultSet.fill_empty(query, control)
        
    if args['m']:
        meme_pwms = MEMEPWMParser(args['m']).get_records() # PWM info. content
        filter_obj = FIMOFilter()
        filter_obj.set_resultset(query)
        filter_obj.set_meme_pwms(meme_pwms)
        
        filter_obj.filter()
        
        #meme_parser = meme_summary.MEMEPWMParser(args['m']) # parse MEME PWMs
        #meme_parser.parse()
#    df = build_skeleton(cont_fname, query_fname)
#    df = populate(df, cont_fname, query_fname)
#    write(df, args['csv'])
