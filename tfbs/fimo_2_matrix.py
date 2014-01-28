'''
Creates a count-matrix given FIMO output using its --text format.
'''

import argparse
import pandas
import tempfile
from pwm_info_content import MEMEPWMParser
from collections import OrderedDict

COLNAME_SEQ = 'Sequence'
COLNAME_TARGET = 'Target'

class FIMOMatrixBuilder():
    def __init__(self, frs_a, frs_b):
        self.control_frs = frs_a
        self.query_frs = frs_b
        self.df = None
        self.build_skeleton()
        
    def build_skeleton(self):
        pwms = set() # ordering is not important
        all_seqs = []
        
        # loop through both files to derive row-count
        for dataset in [self.control_frs.f, self.query_frs.f]:
            seqs = set()
            for linenum, line in enumerate(open(dataset)):
                if linenum != 0:
                    line = line.strip().split()
                    pwm, seq = line[0: 2]
                    pwms.add(pwm)
                    seqs.add(seq)
            all_seqs.extend(seqs)
        
        # wrap counts in a DataFrame
        m = OrderedDict({COLNAME_SEQ: all_seqs})
        m.update({pwm: [0] * len(all_seqs) for pwm in pwms})
        m.update({COLNAME_TARGET: ['None'] * len(all_seqs)})
        self.df = pandas.DataFrame(m, index=all_seqs) # captures count-matrix
    
    def write(self, csv):
        for i, seq in enumerate(self.df[COLNAME_SEQ]):
            seq = self.df[COLNAME_SEQ][i] 
            # replace unique delimiters so both control and query sequences
            # share possible sub-sequences., i.e. match.chr1.111.553
            # is related to chr1.111.553, and so on.
            self.df[COLNAME_SEQ][i] = seq.replace(':', '.').replace('-', '.')
        self.df.to_csv(csv)
    
    def populate(self):
        assert self.df # data-frame must be valid
        for i, dataset in enumerate([self.control_frs.f, self.query_frs.f]):
            for linenum, line in enumerate(open(dataset)):
                if linenum != 0:
                    line = line.strip().split()
                    pwm, seq = line[0: 2]
                    self.df[pwm][seq] += 1 # increment sequence-PWM count
                    self.df[COLNAME_TARGET][seq] = i

class FIMOFilter():
    def __init__(self):
        self.pwm_info = None # references PWM information content
        self.resultset = None # FIMO result-set
        self.pwms_remove = set()
        
    def set_resultset(self, rs):
        assert isinstance(rs, FIMOResultSet)
        self.resultset = rs
    
    def set_meme_pwms(self, pwms):
        self.pwm_info = pwms    
    
    def as_fimo_handle(self):
        fname = self.resultset.f+'.filtered'
        out_handle = open(fname, 'w')
        out_handle.write('#pattern sequence start\n')
        out_handle.flush()
        for seq in self.resultset.mapping:
            seq_hits = self.resultset.mapping[seq]
            for motif in seq_hits:
                for hit in seq_hits[motif]:
                    out_handle.write(motif + ' ' + seq + ' ' + str(hit) + '\n')
                    out_handle.flush()
        out_handle.close()
        return fname
    
    def filter(self):
        ''' 
        Filtering a FIMO file entails the process of determining if there 
        are possibly redundant PWMs. In other words, if there are PWMs that
        each map to the same index on the same sequence, such redundant
        entries need to get removed as they add redundant features to the
        resultant count-matrix.
        '''
        dataset = self.resultset.to_sequence_counts() # PWMs of FIMO file.
        pwm_counts = self.resultset.to_pwm_counts() # number of total PWM hits 
        
        for seq in dataset: # iterate over each sequence
            # print(seq)
            for pwm_a in dataset[seq]:
                for pwm_b in dataset[seq]:
                    if pwm_a != pwm_b:
                        all_hits = dataset[seq] # get all sequence mappings
                        idx_intersect = set(all_hits[pwm_a]).intersection(all_hits[pwm_b])
                        frac_pwm_a, frac_pwm_b = 0.0, 0.0 # PWM A and B ratios
                        if len(idx_intersect) > 0: # e.g. nAB / nA and nAB / nB
                            frac_pwm_a = round(len(idx_intersect) / len(pwm_counts[pwm_a]), 6)
                            frac_pwm_b = round(len(idx_intersect) / len(pwm_counts[pwm_b]), 6)
                        if max(frac_pwm_a, frac_pwm_b) >= 0.25:
                            pwm_a_info = self.pwm_info[pwm_a] # get information
                            pwm_b_info = self.pwm_info[pwm_b]
                            if pwm_a_info >= pwm_b_info: # select worst info
                                self.pwms_remove.add(pwm_a)
                            else:
                                self.pwms_remove.add(pwm_b)

class FIMOResultSet():
    def __init__(self, f):
        self.f = f
        self.mapping = {}
        self.pwms = set()
    
    def parse(self):
        for num, line in enumerate(open(self.f)):
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
        shared_pwms = set(control.get_pwms()).union(query.get_pwms())
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
    
    def hit_count(self):
        num = 0
        for seq in self.mapping:
            for motif in self.mapping[seq]:
                num += len(self.mapping[seq][motif])
        return num
    
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
    
    try:
        control = FIMOResultSet(f = args['control']) # parse the query and control file.
        query = FIMOResultSet(f = args['query']) # parse the query and control file.
        control.parse()
        query.parse()
        FIMOResultSet.fill_empty(query, control)
            
        if args['m']:
            meme_pwms = MEMEPWMParser(args['m']).get_records() # PWM info. content
            filter_cont = FIMOFilter()
            filter_cont.set_resultset(control)
            filter_cont.set_meme_pwms(meme_pwms)
            filter_cont.filter()

            filter_query = FIMOFilter()
            filter_query.set_resultset(query)
            filter_query.set_meme_pwms(meme_pwms)
            filter_query.filter()

            control = FIMOResultSet(filter_cont.as_fimo_handle())
            query = FIMOResultSet(filter_query.as_fimo_handle())
            control.parse()
            query.parse()
            FIMOResultSet.fill_empty(query, control)
        fmb = FIMOMatrixBuilder(control, query)
        fmb.populate()
        fmb.write(csv = args['csv'])
            
    except KeyboardInterrupt:
        print()
