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
    for kw in open(fname):
        kw = kw.strip()  # NCBI Entrez ID
        if len(kw) > 0:
            gi_accession = None  # represents best hit
            search_hits = Entrez.esearch(db=db, term=kw)  # get search hits
            record = Entrez.read(search_hits)
            if int(record['Count']) > 0:
                for gi in record['IdList']:
                    handle = Entrez.efetch(db=db, id=gi, rettype="gb",
                                            retmode="text").read()  # GI entry
                    out_handle = open('genbank.gb', 'w')
                    out_handle.write(handle)  # save genbank contents
                    out_handle.flush()
                    in_handle = SeqIO.read(open('genbank.gb'), 'genbank')
                    if kw in in_handle.description:
                        gi_accession = in_handle  # GI accession has keyword
                        break
            if gi_accession:
                print(kw, '\t', gi_accession.id)  # print matching entrez ID
            else:
                print(kw, '\t', str(None))
        else:
            print(' ', '\t', str(None))
        sys.stdout.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='FILE', required=True,
                            help='File containing keywords [na]')
    parser.add_argument('-db', metavar='db', required=False,
                        choices=['nuccore', 'protein'], default='nuccore',
                        help='NCBI database [nuccore]')
    parser.add_argument('-email', metavar='STR', required=True,
                            help='User-provided email [none]')
    args = vars(parser.parse_args())
    try:
        query_genbank(fname=args['i'], email=args['email'], db=args['db'])
    except KeyboardInterrupt:
        sys.stdout.write('Analysis terminated\n')
