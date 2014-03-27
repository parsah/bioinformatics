'''
A base-bones BED file contains the chromosome ID and indices dedicated to
its respective features. Oftentimes however, a BED entry could also store
additional information such a matrix-like structure representing whether
a feature is found in one feature versus another. This data-structure
would indeed be useful when contrasting how many features are shared or
exclusive to a feature-set. Tools such as bedtools can generate such a
table. Thus, this script takes such outputs and deduces how many BED
entries are shared between different combinations of BED features.
A root column is required. This root serves as the first node (parent)
used to construct this data-structure. Every column from here-on
must be a valid feature.
'''

import argparse
import bedutils
import itertools
import pandas
import os

FEATURE_START = 5  # column within BED file where entries begin from.
NODE_DELIM = '---'  # delimiter for separating nodes.


def product(seq):

    '''
    Given a sequence (list), combinatorial sets within this collection
    are generated ranging from length 1 to the length of the list.
    Within this list will be a node which will serve as a 'root'. This
    root is solely for serve as a hierarchical aide to understand how
    many data-points branch-off from this node. For instance: if your
    sequence is [1,2,3,4], and your root is 2, then every combination
    must have 2 in it. Thus, possible combinations will be:
    [2, 1,2, 1,2,3, etc].

    @param seq: Collection of objects.
    @return: set of valid combinations built around root combination.
    '''

    combs = set()
    for rep in range(len(seq) - 1):
        prods = list(itertools.product(seq, repeat=rep + 1))
        for i in prods:
            if FEATURE_START in i:  # add combination if contains root feature
                i = set(i)
                if i not in combs:  # for eg. ABC is same as BAC or CAB.
                    combs.add(frozenset(i))  # remove these duplicates.

    if len(combs) == 0:  # if root index is < feature start, yield error.
        raise IOError('0 combinations. Root index must be >= feature start.')
    combs.add(frozenset(list(seq)))
    return combs


class MultiBEDEnumerator():
    '''
    A MultiBEDEnumerator object is responsible for the parsing of
    multiIntersectBed files. Such files are essentially BED files
    whereby each entry references which other BED file a specific
    entry is found in. Presence of a BED entry in another BED file
    is measured by a binary (0, 1) value. Instantiations of this class
    therefore facilitate identification of all possible BED files
    and the set of BED entries within each such combination.
    '''

    def __init__(self, f):
        '''
        Constructs an object of type MultiBEDEnumerator.

        @param f: Input BED file generated from multiIntersectBED.
        '''
        self.df = bedutils.parse_multibed(f)  # BED from multiIntersectBed.
        self.combinations = {}  # K => combination, V => a data-frame.
        self.features = list(self.df.columns[FEATURE_START:])  # BED features.

    def debug(self):
        '''
        Print-out important information about the BED file and which
        features are used for analysis.
        '''
        print('Root feature:', self.df.columns[FEATURE_START])  # root feature
        print('All features:')  # show column numbers of all features
        for num, feature in enumerate(self.features):
            print(FEATURE_START + num, '=>', feature)

    def enumerate(self):
        '''
        Analyzes the parsed data-frame and pulls-out columns that map to
        features within the BED file. Features are therefore enumerated so that
        you can discern how many items are mapped to each combination
        (beginning at the parent-most node which is the root).
        '''

        self.debug()  # print helpful information
        print('Building feature combinations ...')
        seq = list(range(FEATURE_START, len(self.features) + FEATURE_START))
        for comb in product(seq):  # per combination, fetch desired columns.
            rows = self.df['chrom'] + ' ' + self.df['start'].map(str) +\
                 ' ' + self.df['end'].map(str)
            df = pandas.DataFrame(self.df[sorted(list(comb))])
            df = df.set_index(rows)
            df = df[df.sum(axis=1) == (df.shape[1])]  # enhancers in all lines
            self.combinations[comb] = df  # lastly, add to central list

    def get_combinations(self):
        '''
        Return a dictionary of combinations and the data-frames that
        reference such combinations.

        @return: dictionary of combinations and data-frames.
        '''
        return self.combinations

    def get_features(self):
        '''
        Returns all features used for combinatorial analysis.

        @return: list of features.
        '''
        return self.features

    def save(self, outdir):
        '''
        Saves all the BED combinations to an output folder.

        @param outdir: User-provided output folder.
        '''
        if not os.path.exists(outdir):  # create output folder
            os.mkdir(outdir)

        for df in list(self.combinations.values()):
            fname = NODE_DELIM.join(list(df.columns)).replace('.bed', '')
            handle = open(outdir + '/' + fname + '.bed', 'w')
            print('  => Writing BED for', fname)
            for idx in df.index:
                handle.write(idx + '\n')
                handle.flush()
            handle.close()


class BEDNetworkBuilder():
    def __init__(self, enum):
        assert isinstance(enum, MultiBEDEnumerator)
        self.enum = enum

    def to_sif(self, outdir):
        handle = open(outdir + '/nodes.sif', 'w')
        combs = sorted(list(self.enum.get_combinations().keys()), key=len)
        for i in combs:  # iterate over all combinations.
            i = sorted(list(i))  # node-list for node A
            for j in combs:
                j = sorted(list(j))  # node-list for node B
                if i != j:
                    if j == i[0: len(i) - 1]:  # node A must be a child of B
                        # get node-names for both node-list A and B
                        node_a = NODE_DELIM.join([self.enum.df.columns[idx]
                                                  for idx in i])
                        node_b = NODE_DELIM.join([self.enum.df.columns[idx]
                                                  for idx in j])
                        handle.write(node_a + ' pp ' + node_b + '\n')
                        handle.flush()
        handle.close()
        print('SIF file saved [OK]')

    def to_node_attrib(self, outdir):
        '''
        Saves node properties to facilitate additional visualizations.
        This resultant properties file contains the node name, the number
        of enhancers this node references, and how many cell-lines make-up
        this cell-line combination. In addition, there are columns for
        keeping track of the proportion of which lines are in the
        combination. For instance, if the line contains [A, C, D], then
        A, C, and D, would each have a score of 0.33, but B and every
        line not A, C, D, would have a score of 0.0.

        @param outdir: Output folder for saving file to.
        '''
        cols = ['Node', 'Num.Enhancers', 'Depth'] + self.enum.get_features()
        idx = list(range(len(self.enum.get_combinations())))  # rows

        # per combination, derive its attributes for further analysis.
        df = pandas.DataFrame(columns=cols, index=idx)
        df = df.fillna(0.0)  # empty-fill matrix
        for rownum, comb in enumerate(self.enum.get_combinations().values()):
            comb_name = NODE_DELIM.join(comb.columns)
            comb_depth = comb_name.count(NODE_DELIM) + 1
            df.loc[rownum, 'Node'] = comb_name  # set column values
            df.loc[rownum, 'Num.Enhancers'] = comb.shape[0]
            df.loc[rownum, 'Depth'] = comb_depth

            # determine proportion of each feature per combination.
            for feat in self.enum.get_features():
                if feat in comb_name:  # if found, derive proportion
                    df.loc[rownum, feat] = round(1.0 / comb_depth, 2)
        df.to_csv(outdir + '/node-attribs.csv', index=False)
        print('Node attributes saved [OK]')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # parameters necessary for combinatorial BED analysis
    parser.add_argument('-bed', metavar='FILE',
                             help='BED from multiIntersectBED [none].')
    parser.add_argument('-o', metavar='DIR', default='./output',
                             help='Folder to save combinations [./output].')
    args = vars(parser.parse_args())
    worker = MultiBEDEnumerator(args['bed'])
    worker.enumerate()
    worker.save(args['o'])
    net_builder = BEDNetworkBuilder(worker)
    net_builder.to_sif(args['o'])
    net_builder.to_node_attrib(args['o'])
