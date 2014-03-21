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
NODE_DELIM = '-' # delimiter for separating nodes.

class MultiBEDWorker():
    def __init__(self, f):
        self.df = pandas.read_table(f) # input BED generated using multiIntersectBed.
        self.df_combinations = [] # list of BED combinations

    def enumerate_bed(self):
        ''' 
        Analyzes the parsed data-frame and pulls-out columns that map to
        features within the BED file. Features are therefore enumerated so that
        you can discern how many items are mapped to each combination
        (beginning at the parent-most node which is the root).
        '''
        
        print('Root feature:', self.df.columns[FEATURE_START]) # display root feature
        print('All features:') # show column numbers of all features
        
        all_cols = list(self.df.columns[FEATURE_START: ]) # get feature columns
        for num, feature in enumerate(all_cols):
            print(FEATURE_START + num, '=>' ,feature)
        
        # next, generate a sequence for deriving feature combinations given root
        print('Building feature combinations ...')
        seq = list(range(FEATURE_START, len(all_cols) + FEATURE_START))
        combs = self.get_combinations(seq)
        
        # for each feature combination, fetch respective columns and enumerate
        for c in combs:
            features = sorted(list(c))
            
            # data-frame row-names are its chromosome and start-end indices. 
            rowname = self.df['chrom'] + ' '  +self.df['start'].map(str) + ' '  + self.df['end'].map(str)
        
            # get columns referencing combination
            selected_df = pandas.DataFrame(self.df[features])
            selected_df = selected_df.set_index(rowname)
            selected_df['rowsum'] = selected_df.sum(axis=1)  # row-wise summation
            
            # filter data so only enhancer found in all lines are saved
            selected_df = selected_df[selected_df['rowsum'] == len(features)]
            selected_df = selected_df.drop('rowsum', axis=1) # sums no longer needed
            self.df_combinations.append(selected_df) # lastly, add to central list
            
        assert len(self.df_combinations) > 0 # combinations must be present

    def get_combinations(self, seq):
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
    
    def save_combinations(self, outdir):
        ''' 
        Saves all the BED combinations to an output folder.
        
        @param outdir: User-provided output folder.
        '''
        if not os.path.exists(outdir): # create output folder
            os.mkdir(outdir)
        
        for df in self.df_combinations:
            fname = NODE_DELIM.join(list(df.columns)).replace('.bed', '')
            handle = open(outdir + '/' + fname + '.bed', 'w')
            print('Writing BED file for', fname)
            for idx in df.index:
                handle.write(idx + '\n')
                handle.flush()
            handle.close()
        print('All analyses complete [OK]')

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
    parser.add_argument('-o', metavar='DIR', default = './output',
                        help='Output folder for saving BED combinations [./output]')
    args = vars(parser.parse_args())
    try:
        worker = MultiBEDWorker(f = args['bed'])
        worker.enumerate_bed()
        worker.save_combinations(outdir = args['o'])
        #to_sif(maps, df)
        #to_node_attrib(maps, df)
        
    except IOError as e:
        print(e)
    except KeyboardInterrupt:
        print()