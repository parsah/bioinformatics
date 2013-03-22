# Useful script for mining NCBI ESTs given species names of preset
# nematodes.

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
            'Xiphinema index'
            'Pratylenchus penetrans', # from cals.ncsu.edu (from here beyond)
            'Pratylenchus vulnus',
            'Pratylenchus brachyrus',
            'Pratylenchus zea',
            'Pratylenchus coffee',
            'Pratylenchus scribneri'
            'Heterodera avenae',
            'Heterodera trifolii',
            'Globodera tabacum',
            'Ditylenchus dipsaci',
            'Ditylenchus destructor',
            'Ditylenchus myceliophagus',
            'Tylenchulus semipenetrans',
            'Helicotylenchus dihystera',
            'Helicotylenchus multicinctus',
            'Belonolaimus longicaudatus']

Entrez.email= raw_input('Enter email:')
num_seqs = 0
for each_org in get_nematodes():
    h = Entrez.esearch(db="nucest",term=each_org+"[Organism]")
    record = Entrez.read(h)
    num_seqs += int(record['Count'])
    for each_entry in record['IdList']:
        res = Entrez.efetch(db="nucleotide", id=each_entry, rettype="fasta", retmode="text")
        print res.read().strip()
