import numpy
import string
import sys
import collections
import argparse

''' 
Takes a tab-delimited text file of real-values and feeds each row into 
the Symbolic Aggregate Approximation (SAX) algorithm.
@author: Parsa Hosseini
'''

class SAXRunner():
    def __init__(self, c, w, a):
        assert w < len(c)
        self.w = w # number of PAA fragments, i.e. length of symbolic string.
        self.n = len(c) # w << i.
        self.a = a # alphabet size.
        self.c = c # data vector.
        self.normed = []
        self.paa_list = []
        
    def paa(self):
        ''' 
        Compute Piecewise Aggregate Approximation (PAA) on a user-provided list
        of Z-normalized values.
        @param c: List of Z-normalized values.
        '''
        for i in range(self.w):
            j = (self.n/self.w) * (i-1)+1
            if j < 0:
                j = 0
            if j <= (self.n/self.w)*i:
                idx_j = int(j)
                self.paa_list.append(round(float(self.c[idx_j]*(self.w/self.n)), 4))
    
    @staticmethod
    def random_normal(i):
        ''' 
        Generates i random numbers from the standard-normal distribution.
        @param i: Number of standard-normal values to be produced.
        @return: sorted values of length i from the standard-normal distribution.
        '''
        # generate a random number where mean is 0 and standard-deviation is 1.0
        vals = list(sorted(numpy.random.normal(loc=0, scale=1.0, size=i)))
        return vals
    
    def normalize(self):
        ''' 
        Z-normalize a vector; i.e. mean of 0 and standard-deviation of 1.
        @param c: List of real values.
        @return: Z-normalized list.
        '''
        mean = numpy.average(self.c)
        sd = numpy.std(self.c)
        v = []
        for i in self.c:
            if sd == 0: # to avoid divide-by-zero, stdev becomes very small
                sd = 1e-5
            normed = (i - mean) / sd
            v.append(normed)
        self.normed = v

    def stratify(self, rand_norm_nums):
        ''' 
        Maps an alphabet character to each segment from the standard-normal 
        function. Mapping is whereby "A" references lowest values and vice-versa.
        '''
        segments = _split_list(rand_norm_nums, parts=self.a)
        segments[0].insert(0, -sys.maxsize) # add upper and lower-bound values.
        segments[-1].append(sys.maxsize)
        chars = list(string.ascii_uppercase)
        strata = collections.OrderedDict() # key => symbol, value => breakpoints.
        for num, segment in enumerate(segments): # characters must be ordered.
            strata[chars[num]] = {'min':min(segment), 'max':max(segment)}
        return strata
    
    def representation(self):
        ''' 
        Maps each PAA onto the standard-normal distribution and identifies
        which strata the PAA maps to. Whichever strata it ends up being,
        the character referencing that strata is then identified.
        @return: symbolic string.
        '''
        sax_out = ''
        for c_i in self.paa_list:
            for s in strata:
                if c_i <= strata[s]['max']:
                    sax_out += s # append the SAX symbol to our string.
                    break
        return sax_out

def _split_list(vals, parts):
    ''' 
    Splits a list into i segments.
    @param vals: List of Z-normalized values
    @param parts: Integer representing number of desired segments.
    @return: collection of lists of maximum size i.
    '''
    frags = []
    frag_len = 1.0/parts*len(vals)
    for i in range(parts):
        frags.append(vals[int(round(i*frag_len)): int(round((i+1)*frag_len))])
    return frags

if __name__ == '__main__':
    rand_norm_nums = SAXRunner.random_normal(i=10000)
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', help='Alphabet size [5]', type=int, 
                        required=False, metavar='INT', default=5)
    parser.add_argument('-w', help='PPA length [6]', type=int, 
                        required=False, default=6)
    parser.add_argument('-in', help='Tab-delimited txt file [na]', 
                        required=True, metavar='FILE')
    args = vars(parser.parse_args())
    
    for line in open(args['in']):
        data = [float(i) for i in line.strip().split('\t')]
        sax = SAXRunner(c=data, w=args['w'], a=args['a'])
        sax.normalize()
        sax.paa()
        strata = sax.stratify(rand_norm_nums)
        sys.stdout.write(sax.representation() + '\n')
