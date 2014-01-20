import argparse
import sys

COLNAME_TARGET = 'Target'
COLNAME_SEQUENCE = 'Sequence'
VERBOSE = False # whether output should be displayed or not. 

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
    
    def set_data(self, df):
        assert isinstance(df, pandas.DataFrame)
        self.df = df
    
    def data(self):
        return self.df
    
    def get_control(self):
        ''' 
        Extract data-frame regions that reference control rows.
        @return: New DataFrame object.
        '''
        cdf = ClassifiableCountDataset()
        cdf.set_data(self.data()[self.data()[COLNAME_TARGET] == 0])
        return cdf
        
    def get_query(self):
        ''' 
        Extract data-frame regions that reference query rows.
        @return: New DataFrame object.
        '''
        cdf = ClassifiableCountDataset()
        cdf.set_data(self.data()[self.data()[COLNAME_TARGET] == 1])
        return cdf

class BinaryClassificationFile():
    def __init__(self, f):
        self.f = f
        self.df = None
        self.is_checked = False
    
    def parse(self):
        if VERBOSE:
            print('Parsing input ...')
        self.df = pandas.read_csv(self.f)
        self.sanity_check()

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
    
    def get_data(self):
        return self.df
    
#    def split_counts(self, nfold):
#         pass
    
#     def get_shape(self):
#         ''' 
#         Return number of rows and columns within data-frame
#         @return: list referencing number of rows and column.
#         '''
#         return self.df.shape
        


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
        file = BinaryClassificationFile(f = args['file']) # encapsulates input CSV
        file.parse() # parse input file
        cdf = ClassifiableCountDataset()
        cdf.set_data(file.get_data())
        print(cdf.get_control())
        
    except IOError as e:
        print(e)