'''
Useful script for mining NCBI ESTs given species names of preset nematodes.
'''

import argparse
import sys
from Bio import Entrez

def get_nematodes():
    return ['Ditylenchus africanus', # from nematode.net (from here beyond)
            'Globodera pallida',
            'Globodera rostochiensis',
            'Heterodera glycines',
            'Heterodera schachtii',
            'Meloidogyne arenaria',
            'Meloidogyne artellia',
            'Meloidogyne chitwoodi',
            'Meloidogyne exigua',
            'Meloidogyne hapla',
            'Meloidogyne incognita',
            'Meloidogyne javanica',
            'Meloidogyne paranaensis',
            'Pratylenchus penetrans',
            'Pratylenchus vulnus',
            'Radopholus similis',
            'Xiphinema index',
            ]
def run(args):
    Entrez.email= args["email"]
    outhandle = open(args["out"], "w")
    num_seqs = 0
    for each_org in get_nematodes():
        h = Entrez.esearch(db="nucest",term=each_org+"[Organism]", RetMax=int(1e9))
        record = Entrez.read(h)
        num_counts = int(record["Count"])
        num_seqs += num_counts
        sys.stdout.write(each_org + " => #/ESTs: " + str(num_counts) + "\n")
        entrez_ids = record["IdList"]
        for num, gi_accn in enumerate(entrez_ids):
            sys.stdout.write("\t" + str(num+1)+'/'+str(len(entrez_ids)) +\
                             " => GI: " + gi_accn + "\n")
            sys.stdout.flush()
            res = Entrez.efetch(db="nucest", id=gi_accn, rettype="fasta", 
                                retmode="text")
            outhandle.write(res.read().strip() + "\n")
            outhandle.flush()
    outhandle.close()
    sys.stdout.write("Analysis complete\n")
    
if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-out', help='Output FASTA [na]', metavar='FILE', required=True)
        parser.add_argument('-email', help='Email [na]', metavar='STR', required=True)
        args = vars(parser.parse_args())
        run(args)
    except KeyboardInterrupt:
        sys.stdout.write("Analysis terminated\n")
