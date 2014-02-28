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
import itertools
import pandas

FEATURE_START = 5 # column within BED file where entries begin from.

def parse_bed(f):
    ''' 
    Parse a user-provided BED file which was produced using the
    multiIntersectBed application within bedtools.

    @param f: BED file.
    @return: pandas data-frame representative of the BED file.
    '''
    
    df = pandas.read_table(f)
    return df

def build_combinations(seq, root):
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
    @param root: object within seq for serving as the parent combination.
    @return: set of valid combinations built around root combination.
    '''
    
    all_prods = set()
    for rep in range(len(seq)-1):
        prods = list(itertools.product(seq, repeat=rep+1))
        for i in prods:
            if root in i: # only add combination if contains root feature
                i = set(i)
                if i not in all_prods: # for eg. ABC is same as BAC or CAB.
                    all_prods.add(frozenset(i)) # remove these duplicates.

    if len(all_prods) == 0: # if root column index is < feature start, yield error. 
        raise IOError('0 combinations made. Root index must be >= feature start.')
    all_prods.add(frozenset(list(seq)))
    return all_prods

def enumerate_features(df, rootcol):
    ''' 
    Analyzes the parsed data-frame and pulls-out columns that map to
    features within the BED file. Features are therefore enumerated so that
    you can discern how many items are mapped to each combination
    (beginning at the parent-most node which is the root).

    @param df: pandas data-frame object.
    @param rootcol: Column index within data-frame serving as the root. 
    '''
    
    print('Root feature:', df.columns[rootcol]) # display root feature
    print('All features:') # show column numbers of all features
    
    all_cols = list(df.columns[FEATURE_START: ]) # get feature columns
    for num, feature in enumerate(all_cols):
        print(FEATURE_START + num, '=>' ,feature)
    
    # next, generate a sequence for deriving feature combinations given root
    print('Building feature combinations ...')
    seq = list(range(FEATURE_START, len(all_cols) + FEATURE_START))
    combs = build_combinations(seq, rootcol)
    
    # for each feature combination, fetch respective columns and enumerate
    maps = {} # key => feature-set, value => count sequences are found in set.
    for c in combs:
        cols = list(c)
        print(c, cols)
        selected_df = pandas.DataFrame(df[cols])
        selected_df['rowsum'] = selected_df.sum(axis=1)  # row-wise summation
        selected_df = selected_df[selected_df.rowsum == len(cols)]
        print(selected_df.shape)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', required=True, metavar='FILE',
                        help='BED from multiIntersectBed (w / headers) [reqd]')
    parser.add_argument('-root', type=int, default=5, metavar='INT',
                        help='Feature column serving as root; zero idx [5]')
    args = vars(parser.parse_args())
    try:
        df = parse_bed(f = args['bed'])
        enumerate_features(df=df, rootcol=args['root'])
    except IOError as e:
        print(e)
    except KeyboardInterrupt:
        print()