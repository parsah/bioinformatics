'''
A LASSO classifier and predictor to aide in classification.
'''

import argparse
import pandas

# variables referencing important columns.
SEQUENCE = 'Sequence'
TARGET = 'Target'

class ClassifiableCountMatrix():
    def __init__(self, m):
        self.matrix = m

class CountMatrixParser():
    ''' 
    A CountMatrixParser parses a user-provided count-matrix and renders
    it as an object fit for manipulation in classification and predictive
    initiatives.
    '''
    
    def __init__(self, f):
        self.f = f
        self.df = None
        
    def parse(self):
        ''' 
        Parse a user-provided count-matrix file. Such parsing renders
        the count-matrix fit for usage in classification and
        predictive endeavors.
        '''
        self.df = pandas.read_csv(self.f)
        # the first column is usually an incremental integer    
        self.df = self.df.drop(labels=self.df.columns[0], axis=1)
        self.df = self.df.set_index(SEQUENCE)
        if not self.is_binary_target():
            raise IOError('Only binary values allowed in target vector.')
    
    def is_binary_target(self):
        ''' 
        Tests whether the count-matrix target is binary. Such testing
        therefore ensures that each observation maps to either a signal
        or background scalar.
        @return: boolean indicative of the target binary state.
        '''
        targets = set(self.df[TARGET])
        return len(targets) == 2
    
    def get_count_matrix(self):
        ''' 
        Return the matrix as an object of type ClassifiableCountMatrix.
        @return: object of type ClassifiableCountMatrix.
        '''
        return ClassifiableCountMatrix(self.df)
    
if __name__ == '__main__':
    desc = 'Classification and prediction using the LASSO algorithm.'
    parser = argparse.ArgumentParser(add_help = False, description=desc)

    # An argument group to encapsulate classification functions
    parser_classify = parser.add_argument_group('Classification')
    parser_classify.add_argument('-classify', metavar='CSV', 
                        help='Input CSV file to classify [na]')
    parser_classify.add_argument('-nfolds', metavar='INT', type=int, 
                        default=5, help='Number of cross-validations [5]')
    parser_classify.add_argument('-iter', metavar='INT', type=int, 
                        default=5, help='Number of iterations [5]')
    
    # An argument group to encapsulate predictive functions
    parser_pred = parser.add_argument_group('Model prediction')
    parser_pred.add_argument('-model', metavar='MODEL',
                        help='Fitted LASSO model [na]')
    parser_pred.add_argument('-test', metavar='CSV',
                        help='Test data to predict [na]')

    # An argument group to encapsulate optional functions    
    parser_opt = parser.add_argument_group('Optional arguments')
    parser_opt.add_argument("-h", "--help", action="help", 
                            help="Show this help message and exit")
    args = vars(parser.parse_args()) # return parsed arguments
    
    try:
        if args['classify']:
            cmp = CountMatrixParser(f = args['classify'])
            cmp.parse()
            matrix = cmp.get_count_matrix()
    except IOError as e:
        print(e)
        