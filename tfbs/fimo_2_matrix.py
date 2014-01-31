'''
A collection of classes which collectively derive a count-matrix given
a baseline and query FIMO output file.
'''

import argparse
import pandas
from pwm_info_content import MEMEPWMParser
from collections import OrderedDict

PROGRESS_CHAR = '*'
PROGRESS_STEP = 5e4

class FIMOMatrixBuilder():
    ''' 
    Constructs a count-matrix given a baseline and query FIMO
    output file. Such a file must be the result of running FIMO in
    --text mode.
    '''
    COLNAME_SEQ = 'Sequence'
    COLNAME_TARGET = 'Target'
    
    def __init__(self, frs_a, frs_b):
        ''' 
        Constructs a bare-bones count-matrix.
        @param control_frs: Control (baseline) FIMOResultSet object.
        @param query_frs: Query FIMOResultSet object.
        '''
        self.control_frs = frs_a
        self.query_frs = frs_b
        self.df = None # references the count-matrix (aka. data-frame)
        self.filter_cont = None # references control FIMOFilter object
        self.filter_query = None # references query FIMOFilter object
        self._build_skeleton()
        
    def _build_skeleton(self):
        ''' 
        Builds the row and column layout for the current matrix object.
        This resultant matrix is empty and is yet to be populated with
        actual count-data.
        '''
        print('Building skeleton-matrix ', end='')
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
                if linenum % PROGRESS_STEP == 0:
                    print(PROGRESS_CHAR, end='')
            all_seqs.extend(seqs)
        print()
        
        # wrap counts in a DataFrame
        m = OrderedDict({FIMOMatrixBuilder.COLNAME_SEQ: all_seqs})
        m.update({pwm: [0] * len(all_seqs) for pwm in pwms})
        m.update({FIMOMatrixBuilder.COLNAME_TARGET: ['None'] * len(all_seqs)})
        self.df = pandas.DataFrame(m, index=all_seqs) # captures count-matrix
    
    def set_filters(self, filter_cont, filter_query):
        ''' 
        Sets FIMOFilter objects which will be used to remove 
        redundant PWMs from the count-matrix.
        @param filter_cont: Control FIMOFilter object.
        @param filter_query: Query FIMOFilter object. 
        '''
        if filter_cont and filter_query:
            assert isinstance(filter_cont, FIMOFilter)
            assert isinstance(filter_query, FIMOFilter)
            self.filter_cont = filter_cont
            self.filter_query = filter_query
    
    def write(self, csv):
        ''' 
        Writes the current count-matrix to a user-provided CSV file.
        @param csv: Output filename.
        '''
        print('Saving output')
        for i, seq in enumerate(self.df[FIMOMatrixBuilder.COLNAME_SEQ]):
            seq = self.df[FIMOMatrixBuilder.COLNAME_SEQ][i] 
            # replace unique delimiters so both control and query sequences
            # share possible sub-sequences., i.e. match.chr1.111.553
            # is related to chr1.111.553, and so on.
            self.df[FIMOMatrixBuilder.COLNAME_SEQ][i] = seq.replace(':', '.').replace('-', '.')
        self.df.to_csv(csv)
    
    def populate(self):
        ''' 
        Adds counts to current matrix object, a process necessary to
        ultimately build a functional object.
        '''
        assert self.df # data-frame must be valid
        print('Populating matrix ', end='')
        for i, dataset in enumerate([self.control_frs.f, self.query_frs.f]):
            for linenum, line in enumerate(open(dataset)):
                if linenum != 0:
                    line = line.strip().split()
                    pwm, seq = line[0: 2]
                    self.df[pwm][seq] += 1 # increment sequence-PWM count
                    self.df[FIMOMatrixBuilder.COLNAME_TARGET][seq] = i
                if linenum % PROGRESS_STEP == 0:
                    print(PROGRESS_CHAR, end='')
        print()
                
        if self.filter_cont and self.filter_query: # if valid filters provided
            pwms_del = set(self.filter_cont.pwms_remove).union(self.filter_query.pwms_remove)
            self.df.drop(list(pwms_del), axis=1) # delete PWMs flagged for removal

class FIMOFilter():
    ''' 
    Oftentimes, a PWM may map to the same sequence indices as another PWM.
    In such instances, there is no need to keep both PWMs since the 
    mapping dynamics of the former can be captured by the latter. By
    filtering for such redundancy, the feature space is reduced; allowing
    for easier feature discrimination and user interpretation.
    '''
    
    def __init__(self):
        ''' 
        Creates a filter object that captures redundant PWMs. 
        '''
        self.pwm_info = None # references PWM information content
        self.resultset = None # FIMO result-set
        self.pwms_remove = set() # references PWMs to be removed.
        
    def set_resultset(self, rs):
        ''' 
        Set the FIMOResultSet object utilized by this object.
        @param rs: FIMOResultSet object.
        '''
        assert isinstance(rs, FIMOResultSet)
        self.resultset = rs
    
    def set_meme_pwms(self, pwms):
        ''' 
        Sets the MEME-format PWMs used to deduce whether one PWM is
        redundant given another. The information content per PWM is
        used as a factor in such determination.
        @param pwms: List of MEMEPWMRecord objects.
        '''
        self.pwm_info = pwms    
    
    def as_fimo_handle(self):
        ''' 
        Following filtering, the filtered-file is to be saved so that
        it can be fed back into the analysis pipeline so that a count-matrix
        can be built around this high-quality data-set. 
        '''
        fname = self.resultset.f+'.filtered'
        out_handle = open(fname, 'w')
        out_handle.write('#pattern sequence start\n') # headers.
        out_handle.flush()
        for seq in self.resultset.mapping:
            seq_hits = self.resultset.mapping[seq]
            for motif in seq_hits:
                for hit in seq_hits[motif]: # write each entry to file.
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
                                print('[INFO]',pwm_a, 'is flagged for removal')
                            else:
                                self.pwms_remove.add(pwm_b)
                                print('[INFO]',pwm_b, 'is flagged for removal')

class FIMOResultSet():
    ''' 
    Models contents of a FIMO output file produced from running fimo in 
    --text mode. There may be instances in such a file whereby a PWM can
    hit the same sequence index multiple times. In such an instance, there
    may be a differing p-value. This class only models one index and not
    all such indices.
    '''
    def __init__(self, f):
        ''' 
        Builds an object of type FIMOResultSet strictly based on a 
        user-provided FIMO --text mode input file.
        @param f: Input FIMO file.
        @param mapping: Hash of sequences and their mapped PWMs.
        @param pwms: Collection of PWMs mapping to file.
        '''
        self.f = f
        self.mapping = {}
        self.pwms = set()
    
    def parse(self):
        ''' 
        Parses a FIMO file so as to keep track of what PWMs mapped
        to what sequences, and how many unique PWMs there are in the
        input file.
        '''
        print('Parsing', self.f, '' ,end='')
        for num, line in enumerate(open(self.f)):
            line = line.strip().split()
            if num != 0:
                motif_id, seq, hit = line[0: 3]
                if seq not in self.mapping:
                    self.mapping[seq] ={}
                if motif_id not in self.mapping[seq]:
                    self.mapping[seq][motif_id] = set()
                self.mapping[seq][motif_id].add(int(hit))
                if motif_id not in self.pwms:
                    self.pwms.add(motif_id)
            if num % PROGRESS_STEP == 0:
                print(PROGRESS_CHAR, end='')
        print()
    
    def to_pwm_counts(self):
        ''' 
        Unlike the default mapping, this behavior enables derivation of
        a hash that references where each PWM maps across the entire 
        sequence. This function is useful in instances in-which you wish
        to investigate similarity between where PWM x maps with PWM y.
        @return: dictionary; key: PWM, value: set of mapped indices.
        '''
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
    
    try:
        # === parse the user-provided FIMO files === 
        control = FIMOResultSet(f = args['control']) # parse the query and control file.
        query = FIMOResultSet(f = args['query']) # parse the query and control file.
        control.parse()
        query.parse()
        filter_query, filter_cont = None, None # reference FIMO filteration
        
        # === perform FIMO filtering on query, control if requested ===
        if args['m']:
            meme_pwms = MEMEPWMParser(args['m']).get_records() # PWM info. content
            filter_cont = FIMOFilter()
            filter_cont.set_resultset(control)
            filter_cont.set_meme_pwms(meme_pwms)
            filter_cont.filter()

            filter_query = FIMOFilter() # filter query file.
            filter_query.set_resultset(query)
            filter_query.set_meme_pwms(meme_pwms)
            filter_query.filter()

            # === parse the FIMO-filtered file ===
            control = FIMOResultSet(filter_cont.as_fimo_handle())
            query = FIMOResultSet(filter_query.as_fimo_handle())
            control.parse()
            query.parse()
        
        # === build a count-matrix and save results to file === 
        fmb = FIMOMatrixBuilder(control, query)
        fmb.set_filters(filter_cont, filter_query) # set the filters
        fmb.populate()
        fmb.write(csv = args['csv'])
            
    except KeyboardInterrupt:
        print()
