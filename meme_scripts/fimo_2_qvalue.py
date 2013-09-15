import argparse
import os
import rpy2.robjects as robjects

motif_ext = '.tmp' # filename extension of motif files

def get_motifs(fname):
    ''' 
    Get all available motifs within a FIMO output file.
    @param fname: Input FIMO filename.
    '''
    motifs = {}
    for linenum, i in enumerate(open(fname)):
        if linenum > 0:
            i = i.strip().split('\t')
            m = i[0]
            if m not in motifs:
                motifs[m] = open(m+motif_ext, 'w') # file to store p-values
    return motifs

def parse_file(fname):
    ''' 
    Parses the actual FIMO file and partitions into contents into a file for
    each motif.
    @param fname: Input FIMO filename.
    '''
    motifs = get_motifs(fname)
    for linenum, i in enumerate(open(fname)):
        if linenum > 0:
            i = i.strip().split('\t')
            m = i[0]
            motifs[m].write('\t'.join(i) + '\n')
            motifs[m].flush()
    for i in motifs:
        motifs[i].close()
        
def adjust(p_cutoff, m):
    ''' 
    Performs p-value adjustment on a given set of p-values derived for each
    motif.
    @param p_cutoff: Adjusted p-value cutoff.
    @param m: P-value adjustment method.
    '''
    all_files = os.listdir(path=os.curdir)
    motifs = [os.curdir + '/' + i for i in all_files if i.endswith(motif_ext)]
    for each_file in motifs:
        pvals = [] # references p-values
        for each_line in open(each_file):
            each_line = each_line.strip().split('\t')
            pvals.append(float(each_line[6]))
        r_adjust = robjects.r['p.adjust']
        adj_pvals = list(r_adjust(robjects.FloatVector(pvals), method=m))
        
        # once adjusted p-values have been derived, set to their entry
        for linenum, each_line in enumerate(open(each_file)):
            each_line = each_line.strip().split('\t')
            adj_pval = round(adj_pvals[linenum], 4)
            if adj_pval < p_cutoff:
                each_line[7] = str(adj_pval) # set adjusted p-value
                print('\t'.join(each_line))
    clean() # remove all temporary files
        
def clean():
    motifs = [os.curdir + '/' + i for i in os.curdir if i.endswith(motif_ext)]
    for i in motifs:
        os.remove(i)

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', metavar='FILE', required=True,
                            help='Output file with --text mode from FIMO [na]')
        parser.add_argument('-pv', metavar='FLOAT', required=False, default=0.05,
                            help='Adjusted p-value cutoff [0.05]')
        parser.add_argument('-method', metavar='STR', required=False,
                            choices=["holm", "hochberg", "hommel", 
                                     "bonferroni", "BH", "BY","fdr", "none"],
                            default='BH', help='Adjusted p-value cutoff [BH]')
        args = vars(parser.parse_args())
        parse_file(fname = args['i'])
        adjust(p_cutoff = args['pv'], m = args['method'])
        
    except KeyboardInterrupt:
        pass
    finally:
        clean()
