'''
'''

import argparse
import pandas
import sklearn


class MatrixParser():
    COLNAME_PWM = 'Unnamed: 0'
    COLNAME_TARGET = 'Target'

    def __init__(self, f, is_classifiable=True):
        self.mat = pandas.read_csv(f)
        self.is_classifiable = is_classifiable
        self.sanity_check()

    def sanity_check(self):
        '''
        Perform sanity checks to ensure that the first and last columns of a
        user-provided CSV file reference PWMs and Targets (if applicable).
        The former check is only applicable if a classifiable matrix is
        provided. Ideally, such a file should have an unnamed first column
        which references PWMs, and an explicitly-named target column that is
        the very last column in the file.
        '''
        print('Sanity-checking matrix ...')
        cols = list(self.mat.columns)
        has_pwm = cols[0] == MatrixParser.COLNAME_PWM  # PWMs must be found.
        has_target = cols[-1] == MatrixParser.COLNAME_TARGET
        if self.is_classifiable:  # targets, PWMs needed in classifiable files.
            if not all([has_pwm, has_target]):
                raise IOError('PWMs, Targets must be first and last columns.')
        else:
            if not has_pwm:
                raise IOError('PWMs must be the first column.')

    def parse(self):
        print('Parsing matrix ...')
        self.mat = self.mat.set_index(MatrixParser.COLNAME_PWM, drop=True)
        self.mat.index.names = ['PWM']
        if 'Sequence' in self.mat.columns:  # drop column already indexed.
            self.mat = self.mat.drop(['Sequence'], axis=1)

    def debug(self):
        n_rows, n_cols = self.mat.shape[0], self.mat.shape[1]
        print('  => {0} rows x {1} columns'.format(n_rows, n_cols))

    def to_classifiable_matrix(self):
        return ClassifiableMatrix(self.mat)

#         print(self.mat[self.mat['Target'] == 1])


class ClassifiableMatrix():
    def __init__(self, m):
        assert isinstance(m, pandas.DataFrame)
        self.mat = m


def main(args):
    if args['in']:
        m = MatrixParser(f=args['in'])
        m.parse()
        m.debug()
        cm = m.to_classifiable_matrix()


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
