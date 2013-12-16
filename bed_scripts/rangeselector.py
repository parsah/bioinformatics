'''
Given a BED file and a genomic FASTA file, each genomic region in the BED file 
is paired with an equal-length genomic segment, matching its GC and repeat
content. Selecting such sequences therefore helps extract sequences 
representative of a BED file which match its intrinsic properties.
'''

import argparse

class RangeSelector():
    ''' 
    Given genomic ranges defined in a user-provided BED file, this class
    selects corresponding genomic sequences having equivalent GC and
    repetitive properties as the genomic range. These resultant sequences
    are therefore equivalent in intrinsic properties and suitable for 
    contrast-based analyses.
    '''
    def __init__(self):
        pass
    

if __name__ == '__main__':
    desc = 'Script to match BED entries with feature-specific genomic regions.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-bed', metavar='FILE', help='BED file [req]',
                        required=True)
    parser.add_argument('-fasta', metavar='FILE', help='FASTA file [req]',
                        required=True)
    parser.add_argument('-o', metavar='FILE', help='FASTA output file [./sequences.fasta]',
                        default='./sequences.fasta')
    args = vars(parser.parse_args())