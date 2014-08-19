"""
Trivial script that enumerates frequency of properties in a GTF file.
"""

import argparse


def run_summary(gtf):
    """
    Derives a summary of a user-provided GTF file; enumerates only the
    source column of the respective file.
    """
    data = {}
    for i in open(gtf):
        i = i.strip().split('\t')
        source = i[1]
        if source not in data:
            data[source] = 0
        data[source] += 1

    # sort the dictionary by its values.
    data = sorted(data.items(), reverse=True, key=lambda x: x[1])
    for i in data:
        print(i[0], i[1])  # print output

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-gtf', metavar='GTF', required=True,
                        help='GTF annotations file [req]')
    args = vars(parser.parse_args())
    run_summary(args['gtf'])
