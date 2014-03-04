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
NODE_DELIM = '-' # delimiter for separating nodes.

def parse_bed(f):
    ''' 
    Parse a user-provided BED file which was produced using the
    multiIntersectBed application within bedtools.

    @param f: BED file.
    @return: pandas data-frame representative of the BED file.
    '''
    
    df = pandas.read_table(f)
    return df

def build_combinations(seq):
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
    
    all_prods = set()
    for rep in range(len(seq)-1):
        prods = list(itertools.product(seq, repeat=rep+1))
        for i in prods:
            if FEATURE_START in i: # only add combination if contains root feature
                i = set(i)
                if i not in all_prods: # for eg. ABC is same as BAC or CAB.
                    all_prods.add(frozenset(i)) # remove these duplicates.

    if len(all_prods) == 0: # if root column index is < feature start, yield error. 
        raise IOError('0 combinations made. Root index must be >= feature start.')
    all_prods.add(frozenset(list(seq)))
    return all_prods

def enumerate_nodes(df):
    ''' 
    Analyzes the parsed data-frame and pulls-out columns that map to
    features within the BED file. Features are therefore enumerated so that
    you can discern how many items are mapped to each combination
    (beginning at the parent-most node which is the root).

    @param df: pandas data-frame object.
    @param rootcol: Column index within data-frame serving as the root. 
    '''
    
    print('Root feature:', df.columns[FEATURE_START]) # display root feature
    print('All features:') # show column numbers of all features
    
    all_cols = list(df.columns[FEATURE_START: ]) # get feature columns
    for num, feature in enumerate(all_cols):
        print(FEATURE_START + num, '=>' ,feature)
    
    # next, generate a sequence for deriving feature combinations given root
    print('Building feature combinations ...')
    seq = list(range(FEATURE_START, len(all_cols) + FEATURE_START))
    combs = build_combinations(seq)
    
    # for each feature combination, fetch respective columns and enumerate
    mappings = {}
    for c in combs:
        cols = sorted(list(c)) 
        selected_df = pandas.DataFrame(df[cols])
        selected_df['rowsum'] = selected_df.sum(axis=1)  # row-wise summation
        selected_df = selected_df[selected_df.rowsum == len(cols)]
        num_rows = selected_df.shape[0]
        mappings[tuple(cols)] = num_rows
    return mappings

def to_sif(d, df):
    handle = open('nodes.sif', 'w')
    keys = sorted(list(d.keys()), key=len)
    for i in keys:
        i = list(i)
        for j in keys:
            j = list(j)
            if i != j:
                if j == i[0: len(i)-1]:
                    node_a = NODE_DELIM.join([df.columns[idx] for idx in i])
                    node_b = NODE_DELIM.join([df.columns[idx] for idx in j])
                    handle.write(node_a + ' pp ' +node_b + '\n')
                    handle.flush()
    handle.close()
    print('SIF file saved [OK]')
    
def to_node_attrib(d, df):
    handle = open('nodes.csv', 'w')
    handle.write('Node,Enhancer.Count\n')
    handle.flush()
    
    for i in d:
        node = NODE_DELIM.join([df.columns[idx] for idx in i])
        handle.write(node + ',' + str(d[i]) + '\n')
        handle.flush()
    handle.close()
    print('Node attributes saved [OK]')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', required=True, metavar='FILE',
                        help='BED from multiIntersectBed (w / headers) [reqd]')
    args = vars(parser.parse_args())
    try:
        df = parse_bed(f = args['bed'])
        maps = enumerate_nodes(df=df)
        to_sif(maps, df)
        to_node_attrib(maps, df)
        
    except IOError as e:
        print(e)
    except KeyboardInterrupt:
        print()