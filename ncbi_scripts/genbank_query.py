'''
A simple script for querying NCBI GenBank given a list of keywords. The 
GenBank record with the keyword in its description is assigned to the keyword.
'''

from Bio import Entrez
from Bio import SeqIO
import argparse
import sys

def query_genbank(fname, db, email):
    ''' 
    Query NCBI GenBank given a filename of user-provided keywords and a
    user-provided email address
    @param fname: User-provided filename of keywords
    @param email: User-provided email address
    '''
    Entrez.email = email
    for i, kw in enumerate(open(fname)):
        kw = kw.strip()
        if len(kw) > 0:
            gi_accession = None
            search_hits = Entrez.esearch(db = db, term = kw) # get search hits
            record = Entrez.read(search_hits)
            if record['Count'] > 0:
                for gi in record['IdList']:
                    handle  = Entrez.efetch(db = db, id=gi, rettype="gb", 
                                            retmode="text").read() # get GI entry
                    out_handle = open('genbank.gb', 'w')
                    out_handle.write(handle) # save genbank contents
                    out_handle.flush()
                    in_handle = SeqIO.read(open('genbank.gb'), 'genbank')
                    if kw in in_handle.description:
                        gi_accession = in_handle # GI accession contains keyword
                        break
            if gi_accession:
                print kw, '->' ,gi_accession.description
            else:
                print kw, '->' ,str(None)
        else:
            print ' ', '->' ,str(None)
                
        
        sys.stdout.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='FILE', required=True,
                            help='File containing keywords [na]')
    parser.add_argument('-db', metavar='db', required=False, 
                        choices = ['nuccore', 'protein'], default='nuccore',
                        help='NCBI database [nuccore]')
    parser.add_argument('-email', metavar='STR', required=True,
                            help='User-provided email [none]')
    args = vars(parser.parse_args())
    query_genbank(fname = args['i'], email = args['email'], db = args['db'])

