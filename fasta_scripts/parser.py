''' 
A high-level module referencing functions used in handling FASTA files.
'''

def parse_fasta(fname):
    ''' 
    Trivial function to parse a user-provided FASTA file
    '''
    try:
        seqs, header = {}, ''
        for line in open(fname):
            line = line.strip()
            if line.startswith('>'):
                header = line[1:] # remove the fasta header
                seqs[header] = ''
            else:
                seqs[header]+=line
        return seqs
    except IOError:
        print('File:', fname, 'not found [error]')