'''
'''

import argparse
import pandas
import time
from sklearn import linear_model
from sklearn import cross_validation
from sklearn import metrics


class MatrixParser():
    '''
    The MatrixParser class provides behaviors and states for parsing a matrix
    CSV file. Such files must have their first column represent desired
    observations. The very last column must explicitly represent a binary
    target vector (0, 1).
    '''
    COLNAME_PWM = 'Unnamed: 0'
    COLNAME_TARGET = 'Target'

    def __init__(self, f, is_classifiable=True):
        '''
        Creates an object of type MatrixParser. Each instance of this object
        enables the ability to parse a matrix saved as a CSV file and mold its
        contents to ultimately be classified.
        @param f: Input CSV file.
        @param is_classifiable: Whether has file has a target vector or not.
        '''
        self.mat = pandas.read_csv(f)
        self.is_classifiable = is_classifiable
        self.sanity_check()
        print('Matrix:', f)

    def sanity_check(self):
        '''
        Perform sanity checks to ensure that the first and last columns of a
        user-provided CSV file reference PWMs and Targets (if applicable).
        The former check is only applicable if a classifiable matrix is
        provided. Ideally, such a file should have an unnamed first column
        which references PWMs, and an explicitly-named target column that is
        the very last column in the file.
        '''
        cols = list(self.mat.columns)
        has_pwm = cols[0] == MatrixParser.COLNAME_PWM  # PWMs must be found.
        has_target = cols[-1] == MatrixParser.COLNAME_TARGET
        if self.is_classifiable:  # targets, PWMs needed in classifiable files.
            if not all([has_pwm, has_target]):
                raise IOError('PWMs, Targets must be first and last columns.')
            if set(self.mat[MatrixParser.COLNAME_TARGET]) != set([0, 1]):
                raise IOError('Target vector must be binary (0, 1) only.')
        else:
            if not has_pwm:
                raise IOError('PWMs must be the first column.')

    def parse(self):
        '''
        Parse the CSV file. Following such parsing, the matrix is accessible
        through its column-names and row-names.
        '''
        print('Parsing matrix ...')
        self.mat = self.mat.set_index(MatrixParser.COLNAME_PWM, drop=True)
        self.mat.index.names = ['PWM']
        if 'Sequence' in self.mat.columns:  # drop column already indexed.
            self.mat = self.mat.drop(['Sequence'], axis=1)

    def to_classifiable_matrix(self):
        '''
        Wraps the results following matrix parsing as an object of type
        ClassifiableMatrix. An object of this class is therefore fit for
        matrix manipulation as it is now fit to perform classification on.
        @return: object of type ClassifiableMatrix.
        '''
        cm = ClassifiableMatrix(self.mat)
        if cm.nrows() != cm.get_query().nrows() + cm.get_control().nrows():
            raise IOError('Error in matrix; ensure target is binary.')
        return cm


class ClassifiableMatrix():
    '''
    A ClassifiableMatrix enables the manipulation of a matrix data-structure
    with intention of classification. Such an object must have observations
    as designated by row-names, and features as designated by column names.
    A corresponding binary target vector need also be present; found as the
    last column in the matrix.
    '''
    def __init__(self, m):
        '''
        Creates an object of type ClassifiableMatrix.
        @param m: pandas Dataframe object generated from MatrixParser.
        '''
        assert isinstance(m, pandas.DataFrame)
        self.mat = m

    def debug(self):
        n_query, n_cont = self.get_query().nrows(), self.get_control().nrows()
        print('[ {0} rows x {1} columns ]'.format(self.nrows(), self.ncols()))
        print('[ => {0} query, {1} control ]'.format(n_query, n_cont))

    def normalize(self, n):
        '''
        In cases whereby a query-set is small, a corresponding control-set
        cannot be of equal size since this control-set may be capture the the
        properties of the global system. Thus, rather than generating 1 control
        for each query, you may have to generate N controls. Doing so could
        help capture the systematic properties of the query by giving it a
        large control-set to be contrasted against. This constant, N, is what
        is used to 'amplify' counts and must be normalized so as to not skew
        classifications.
        '''
        control = self.mat[MatrixParser.COLNAME_TARGET] == 0  # get controls
        self.mat.loc[control] = self.mat.loc[control] / n  # normalize

    def get_control(self):
        '''
        Retrieve the background observations which are observations referencing
        a target value of 0.
        @return: ClassifiableMatrix object referencing control observations.
        '''
        control = self.mat[self.mat[MatrixParser.COLNAME_TARGET] == 0]
        return ClassifiableMatrix(control)

    def get_query(self):
        '''
        Retrieve the signal (query) observations which are observations
        referencing a target value of 1.
        @return: ClassifiableMatrix object referencing query observations.
        '''
        query = self.mat[self.mat[MatrixParser.COLNAME_TARGET] == 1]
        return ClassifiableMatrix(query)

    def nrows(self):
        '''
        Retrieve the total number of rows.
        @return: integer row count.
        '''
        return self.mat.shape[0]

    def ncols(self):
        '''
        Retrieve the total number of columns.
        @return: integer column count.
        '''
        return self.mat.shape[1]


class ClassificationFactory():
    '''
    The ClassificationFactory class enables the ability to classify a
    ClassifiableMatrix object. Successful classification is indeed a much
    desired objective and can be obtainable with the help of numerous
    runtime parameters such as cross-validation.
    '''
    def __init__(self, cm, iters, cv, thr, homogen):
        '''
        Creates a ClassificationFactory object, an object designed solely for
        classifying features based on the observation target vector.
        @param m: ClassificationMatrix object.
        @param iters: Number of predictive iterations.
        @param cv: Number of cross-validations.
        @param thr: Predictive threshold cutoff.
        @param homogen: Whether to homogenize counts (see Taher, et. al, 2013).
        '''
        assert isinstance(cm, ClassifiableMatrix)
        self.cm = cm
        self.iters = iters
        self.cv = cv
        self.thr = thr
        self.homogenize = homogen

    def classify(self):
        '''
        Perform classification using the object's ClassifiableMatrix object
        and all classification-specific arguments. Based on how many
        cross-validations are desired, classifiers are built and their
        respective receiver operating characteristic (ROC) curve is
        subsequently produced.
        '''
        print('Building cross-validated classifier ...')
        x = self.cm.mat.drop(MatrixParser.COLNAME_TARGET, 1).values  # counts
        y = self.cm.mat[MatrixParser.COLNAME_TARGET].values  # target vector
        kf = cross_validation.StratifiedKFold(y, self.cv)

        for i, (train_idx, test_idx) in enumerate(kf):
            x_train, x_test = x[train_idx], x[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            classifier = linear_model.LassoCV(normalize=True)
            classifier.fit(x_train, y_train)
            y_preds = classifier.predict(x_test)
            fpr, tpr = metrics.roc_curve(y_test, y_preds)[: -1]
            auc = metrics.auc(fpr, tpr)
            print('[', time.strftime('%I:%M:%S'), '] # ', i + 1, 'of', len(kf), classifier.alpha_)
            kf = cross_validation.StratifiedKFold(y, self.cv)
        #model = linear_model.LassoCV(normalize=True)
        #model.fit(x, y)
        #print(model.alphas_, model.alpha_)

    def get_weights(self):
        df = pandas.DataFrame()
        df['Weights'] = self.classifier.coef_
        df = df.set_index(self.cm.mat.columns[0: -1])
        df.to_csv('weights.csv')


def main(args):
    if args['in']:
        m = MatrixParser(f=args['in'])
        m.parse()
        cm = m.to_classifiable_matrix()
        cm.normalize(n=args['amplify'])  # normalize controls if amplified.
        cf = ClassificationFactory(cm=cm, iters=args['iters'], cv=args['cv'],
                                   thr=args['thr'], homogen=args['homogenize'])
        cf.classify()
        #cf.get_weights()

if __name__ == '__main__':
    desc = 'LASSO classification and prediction.'
    parser = argparse.ArgumentParser(description=desc, add_help=False)
    parser_class = parser.add_argument_group('Classification')
    parser_pred = parser.add_argument_group('Prediction')
    parser_misc = parser.add_argument_group('Miscellaneous')

    # Define classification parameters
    parser_class.add_argument('-in', metavar='CSV',
                              help='Count-matrix with binary targets [na]')
    parser_class.add_argument('-iters', metavar='INT', default=3, type=int,
                              help='# / homogenization iterations [3]')
    parser_class.add_argument('-cv', metavar='INT', default=5, type=int,
                              help='# / cross-validations [5]')
    parser_class.add_argument('-thr', metavar='FLOAT', default=0.5, type=float,
                              help='Prediction probability threshold [0.5]')
    parser_class.add_argument('-amplify', metavar='INT', default=1, type=int,
                              help='Divide pre-amplified control counts [1]')
    parser_class.add_argument('--homogenize', default=False, metavar='',
                              help='Perform homogenization [false]')

    # Define prediction-specific parameters
    parser_pred.add_argument('-counts', metavar='CSV',
                             help='Count matrix for deriving predictions [na]')
    parser_pred.add_argument('-model', metavar='MODEL',
                             help='Classifier object (.model) [na]')

    # Define miscellaneous parameters
    parser_misc.add_argument('-out', metavar='DIR', default='./proj',
                             help='Output folder [./proj]')
    parser_misc.add_argument('-h', action='help',
                             help='Show this help message and exit.')

    args = vars(parser.parse_args())
    try:
        main(args)
    except IOError as e:
        print(e)
    except OSError as e:
        print(e)
    except KeyboardInterrupt:
        print()
