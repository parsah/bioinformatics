import argparse
from Bio import SeqIO

'''
Trivial script to remove a leading adapter from a FASTQ file.
'''

def remove_adapter(fastq, adapter):
    # Parse FASTQ and remove all adapter instances and respective quality.
    records = SeqIO.parse(fastq, 'fastq')
    for entry in records:
        seq = str(entry.seq)
        qscores = [i+33 for i in entry.letter_annotations["phred_quality"]]
        if str(entry.seq).startswith(adapter):
            seq = str(entry.seq)[len(adapter): ]
            qscores = qscores[len(adapter): ]
        qscores = ''.join([chr(i) for i in qscores])
        print('@'+entry.id+'\n'+seq+'\n+\n'+qscores)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', metavar='FASTQ', required=True,
                        help='Input FASTQ file [na]')
    parser.add_argument('-adapter', metavar='STR', required=True,
                        help='5\' DNA adapter to remove [na]')
    args = vars(parser.parse_args())
    remove_adapter(fastq=args['in'], adapter=args['adapter'])
