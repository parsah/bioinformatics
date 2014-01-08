'''
Takes a collection of 'Result' output files once Marina has completed analysis
and creates a single table of IPF values for all TFBSs within each 
user-provided Result file. 
'''

import argparse
import os.path

class Merger():
    def __init__(self, fnames):
        self.resultfiles = [ResultFile(fname) for fname in fnames]
        self.tfbs_set = set() # references union-set of TFBSs
        
    def union(self):
        for rf in self.resultfiles:
            rf.consolidate()
            self.tfbs_set.update(set(list(rf.set_tfbs.keys())))
            
    def merge(self):
        header = ['TFBS'] # very first column references TFBSs
        for rf in self.resultfiles:
            header.append(os.path.basename(rf.fname))
        print('\t'.join(header))        
        for tfbs in self.tfbs_set:
            abundances = ['-'] * len(self.resultfiles)
            for num, rf in enumerate(self.resultfiles):
                if tfbs in rf.set_tfbs:
                    abundances[num] = rf.set_tfbs[tfbs] # set its IPF
            print(tfbs + '\t' + '\t'.join(abundances))

class ResultFile():
    def __init__(self, fname):
        self.fname = fname
        self.set_tfbs = {} # key => TFBS, value => IPF score (column 1)
    
    def consolidate(self):
        handle = open(self.fname)
        for linenum, line in enumerate(handle):
            line = line.strip().split('\t')
            if linenum != 0:
                tfbs, ipf = line[0], line[1]
                self.set_tfbs[tfbs] = ipf

if __name__ == '__main__':
    desc = 'Script for merging IPF-normalized Marina output (Result) files.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-in', metavar='', nargs='+', required=True,
                        help='List of Marina Result files [na]')
    args = vars(parser.parse_args())
    merger = Merger(fnames = args['in'])
    merger.union()
    merger.merge()
    
    
