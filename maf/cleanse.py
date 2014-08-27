import sys
import argparse

def clean_maf(f):
    for line in open(f):
        line = line.strip()
        if line == '' or line.startswith('a') or line.startswith('s'):
            print(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-maf', metavar='FILE', required=True, help='MAF file [reqd]')
    args = vars(parser.parse_args())
    clean_maf(f=args['maf'])
