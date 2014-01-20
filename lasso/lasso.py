import argparse
import sys

COLNAME_TARGET = 'Target'
COLNAME_SEQUENCE = 'Sequence'
VERBOSE = False

try:
    import pandas
    from sklearn.linear_model import LassoCV
    from sklearn.metrics import auc
except ImportError as e:
    print(format(e))
    sys.exit(1) # if no errors found, exit program

class LASSOClassifier():
    def __init__(self):
        pass
    
    def classify(self):
        pass
    
    def get_auc(self):
        pass

class ClassifiableCountDataset():
    def __init__(self):
        self.df = None
        self.is_checked = False
    
    def parse(self, f):
        if VERBOSE:
            print('Parsing input ...')
        self.df = pandas.read_csv(f)
        self.sanity_check()
    
    def get_data(self):
        return self.df
    
    def split_counts(self, nfold):
        pass
    
    def get_control(self):
        pass
    
    def get_query(self):
        pass
    
    def sanity_check(self):
        ''' 
        Performs several tests to determine if the matrix is suitable
        for analytical purposes. Such tests include whether a 'Sequence'
        column is present and references 
        '''
        if VERBOSE:
            print('Running matrix sanity checks ...')
        if COLNAME_TARGET in self.df:
            target_col = set(self.df[COLNAME_TARGET])
            if target_col != set([0, 1]) and COLNAME_SEQUENCE not in self.df:
                raise IOError('Target must be binary [error]')
        else:
            raise IOError('a \'Target\' vector is required [error]')
        
        if COLNAME_SEQUENCE in self.df:
            if self.df[COLNAME_SEQUENCE].dtype != object:
                raise IOError('Sequence vector must be non-numeric [error]')
        else:
            raise IOError('a \'Sequence\' column is required [error]')
        self.df = self.df.set_index(COLNAME_SEQUENCE, drop=False) # index to sequences
        del self.df['Unnamed: 0']
        self.is_checked = True
        
    def is_sanitized(self):
        return self.is_checked

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', metavar='FILE', required=True,
                        help='Input CSV matrix [req].')
    parser.add_argument('-nfolds', metavar='INT', type=int,
                        default=5, help='Number of folds [5].')
    parser.add_argument('--verbose', default=False, action='store_true',
                        help='Verbose output [false]')
    parser.add_argument('-outdir', metavar='DIR', default='./output/',
                        help='Output project folder [./output].')
    try:
        args = vars(parser.parse_args())
        VERBOSE = args['verbose']
        dataset = ClassifiableCountDataset() # encapsulates input CSV
        dataset.parse(f = args['file']) # parse input file
        
    except IOError as e:
        print(e)