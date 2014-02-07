'''
A LASSO classifier and predictor to aide in classification.
'''

import argparse
import pandas
from sklearn import linear_model

# variables referencing important columns.
SEQUENCE = 'Sequence'
TARGET = 'Target'

class LASSOClassifier():
    ''' 
    A LASSOClassifier object performs classification so as to 
    ultimately discriminate between features and their propensity
    to be found in either query or baseline data-sets.
    '''
    
    def __init__(self, m, n, i):
        ''' 
        Creates a new LASSOClassifier object.
        @param m: ClassifiableCountMatrix object.
        @param n: number of cross-validation folds.
        @param i: number of iterations used for predictive purposes.
        '''
        assert isinstance(m, ClassifiableCountMatrix)
        self.matrix = m
        self.nfolds = n
        self.iters = i
        self.model = None
    
    def classify(self):
        assert isinstance(self.matrix, ClassifiableCountMatrix)
        print('Building classifier ...')
        # pull-out matrix counts exclusing targets (x) and targets (y)
        x = self.matrix.df[self.matrix.df.columns[0: -2]].as_matrix()
        y = self.matrix.df[TARGET].as_matrix()
        self.model = linear_model.LassoCV(tol = 0.05, cv=self.nfolds)
        self.model.fit(x, y)
        
    def get_coefficients(self):
        coef_vector = self.model.coef_
        features = self.matrix.df.columns[0: -2]
        df = pandas.DataFrame({'PWM': features, 'Weights': coef_vector})
        #df = df.set_index('PWM')
        return df
    
    def get_ratios(self):
        pass
    
    def homogenize(self):
        pass
    
    def prediction_vector(self):
        pass

class ClassifiableCountMatrix():
    ''' 
    A ClassifiableCountMatrix encapsulates a count-matrix which is
    to be used for classification purposes. The value of this matrix
    comes from its ability to have a target vector which maps each
    feature to either a signal or background data-space.
    '''
    
    def __init__(self, df):
        ''' 
        Creates a new ClassifiableCountMatrix.
        @param df: References a pandas data-frame object.
        '''
        assert isinstance(df, pandas.DataFrame)
        self.df = df
        
    def get_query(self):
        ''' 
        Returns on the features mapped as a query.
        @return: new matrix referencing only query features.
        '''
        return ClassifiableCountMatrix(self.df[self.df[TARGET] == 1])
    
    def get_baseline(self):
        '''
        Returns on the features mapped as a baseline (control).
        @return: new matrix referencing only baseline features.
        '''
        return ClassifiableCountMatrix(self.df[self.df[TARGET] == 0])
    
    def sample(self, p):
        ''' 
        Samples the matrix based on a predefined fraction.
        '''
    
    def __str__(self):
        return str(self.df)
    

class CountMatrixParser():
    ''' 
    A CountMatrixParser parses a user-provided count-matrix and renders
    it as an object fit for manipulation in classification and predictive
    initiatives.
    '''
    
    def __init__(self, f):
        ''' 
        Given an input filename, create an object of type 
        CountMatrixParser.
        @param f: Input CSV file.
        '''
        self.f = f
        self.df = None
        
    def parse(self):
        ''' 
        Parse a user-provided count-matrix file. Such parsing renders
        the count-matrix fit for usage in classification and
        predictive endeavors.
        @return: object of type ClassifiableCountMatrix.
        '''
        print('Parsing file ...')
        self.df = pandas.read_csv(self.f)
        # the first column is usually an incremental integer; delete this.
        self.df = self.df.drop(labels=self.df.columns[0], axis=1)
        self.df = self.df.set_index(SEQUENCE)
        if not self.is_binary_target():
            raise IOError('Only binary values allowed in target vector.')
        return ClassifiableCountMatrix(self.df)
    
    def is_binary_target(self):
        ''' 
        Tests whether the count-matrix target is binary. Such testing
        therefore ensures that each observation maps to either a signal
        or background scalar.
        @return: boolean indicative of the target binary state.
        '''
        targets = set(self.df[TARGET])
        return targets == set([0, 1]) # only 0 and 1 allowed.
    
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
    parser_opt.add_argument("-o", metavar='DIR', default='./output',
                            help="Output folder [./output]")
    
    args = vars(parser.parse_args()) # return parsed arguments
    
    try:
        if args['classify']:
            matrix = CountMatrixParser(f = args['classify']).parse()
            lasso = LASSOClassifier(matrix, args['nfolds'], args['iter'])
#             lasso.classify()
            lasso.get_ratios()
    except IOError as e:
        print(e)
    except KeyboardInterrupt as e:
        print()
        