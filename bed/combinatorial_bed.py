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
import os

FEATURE_START = 5 # column within BED file where entries begin from.
NODE_DELIM = '---' # delimiter for separating nodes.

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
    for rep in range(len(seq)-1):
        prods = list(itertools.product(seq, repeat=rep+1))
        for i in prods:
            if FEATURE_START in i: # only add combination if contains root feature
                i = set(i)
                if i not in combs: # for eg. ABC is same as BAC or CAB.
                    combs.add(frozenset(i)) # remove these duplicates.

    if len(combs) == 0: # if root column index is < feature start, yield error. 
        raise IOError('0 combinations made. Root index must be >= feature start.')
    combs.add(frozenset(list(seq)))
    return combs

class MultiBEDEnumerator():
    def __init__(self, f):
        self.df = pandas.read_table(f) # input BED generated using multiIntersectBed.
        self.combinations = {} # K => combination, V => corresponding data-frame.
        self.features = [] # list of features comprising the BED file.
    
    def enumerate(self):
        ''' 
        Analyzes the parsed data-frame and pulls-out columns that map to
        features within the BED file. Features are therefore enumerated so that
        you can discern how many items are mapped to each combination
        (beginning at the parent-most node which is the root).
        '''

        print('Root feature:', self.df.columns[FEATURE_START]) # root feature
        print('All features:') # show column numbers of all features
        self.features = list(self.df.columns[FEATURE_START: ]) # feature columns
        for num, feature in enumerate(self.features):
            print(FEATURE_START + num, '=>' ,feature)
        
        print('Building feature combinations ...')
        seq = list(range(FEATURE_START, len(self.features) + FEATURE_START))
        
        for comb in product(seq): # per combination, fetch the desired columns.
            rows = self.df['chrom'] + ' '  + self.df['start'].map(str) + ' ' + self.df['end'].map(str)
            
            df = pandas.DataFrame(self.df[sorted(list(comb))]) # get features
            df = df.set_index(rows)
            df['rowsum'] = df.sum(axis=1)  # row-wise summation
                        
            df = df[df['rowsum'] == (df.shape[1]-1)] # enhancer to be in all lines  
            df = df.drop('rowsum', axis=1) # sums no longer needed
            self.combinations[comb] = df # lastly, add to central list

    def get_combinations(self):
        return self.combinations
    
    def get_features(self):
        return self.features
    
    def save(self, outdir):
        ''' 
        Saves all the BED combinations to an output folder.
        
        @param outdir: User-provided output folder.
        '''
        if not os.path.exists(outdir): # create output folder
            os.mkdir(outdir)
        
        for df in list(self.combinations.values()):
            fname = NODE_DELIM.join(list(df.columns)).replace('.bed', '')
            handle = open(outdir + '/' + fname + '.bed', 'w')
            print('Writing BED file for', fname)
            for idx in df.index:
                handle.write(idx + '\n')
                handle.flush()
            handle.close()
        print('All analyses complete [OK]')

class BEDNetworkBuilder():
    def __init__(self, enum):
        assert isinstance(worker, MultiBEDEnumerator)
        self.enum = enum

    def to_sif(self, outdir):
        handle = open(outdir + '/nodes.sif', 'w')
        keys = sorted(list(self.enum.get_combinations().keys()), key=len)
        for i in keys:
            i = list(i)
            for j in keys:
                j = list(j)
                if i != j:
                    if j == i[0: len(i)-1]:
                        node_a = NODE_DELIM.join([self.enum.df.columns[idx] for idx in i])
                        node_b = NODE_DELIM.join([self.enum.df.columns[idx] for idx in j])
                        handle.write(node_a + ' pp ' +node_b + '\n')
                        handle.flush()
        handle.close()
        print('SIF file saved [OK]')
    
    def to_node_attrib(self, outdir):    
        cols = ['Node', 'Num.Enhancers', 'Depth'] + self.enum.get_features()
        idx = list(range(len(self.enum.get_combinations())))
        
        # per combination, derive its attributes for further analysis.
        df = pandas.DataFrame(columns = cols, index = idx)
        for rownum, df_comb in enumerate(self.enum.get_combinations().values()):
            node_name = NODE_DELIM.join(df_comb.columns)
            node_depth = node_name.count(NODE_DELIM) + 1
            df.loc[rownum, 'Node'] = node_name
            df.loc[rownum, 'Num.Enhancers'] = df_comb.shape[0]
            df.loc[rownum, 'Depth'] = node_depth
            
            # determine proportion of each feature per combination.
            for feat in self.enum.get_features():
                if feat in node_name: # if found, derive proportion
                    df.loc[rownum, feat] = round(1.0 / node_depth, 2)
                else: # if not found, proportion is 0.0.
                    df.loc[rownum, feat] = 0.0
                                
        df.to_csv(outdir + '/node-attribs.csv')
        print('Node attributes saved [OK]')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', required=True, metavar='FILE',
                        help='BED from multiIntersectBed (w / headers) [reqd]')
    parser.add_argument('-o', metavar='DIR', default = './output',
                        help='Output folder for saving BED combinations [./output]')
    args = vars(parser.parse_args())
    try:
        worker = MultiBEDEnumerator(f = args['bed'])
        worker.enumerate()
        worker.save(outdir = args['o'])

        net_builder = BEDNetworkBuilder(worker)
        net_builder.to_sif(outdir = args['o'])
        net_builder.to_node_attrib(outdir = args['o'])
        
    except IOError as e:
        print(e)
    except KeyboardInterrupt:
        print()