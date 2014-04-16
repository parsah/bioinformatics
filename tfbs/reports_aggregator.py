'''
Since there is redundancy is many PWMs, there may be instances whereby one
PWM variant is enriched at one instance but another. This script helps to
identify which PWM variants are the most maximum so that this specific
instance is saved.
'''

import pandas
import argparse

DELIM = '..'  # tfSearch PWM IDs are delimited by 2x periods.


def parse_report(f):
    '''
    Parses a LASSO Reports file generated using the /lasso/sif.R function.

    @param f: LASSO Reports csv file.
    @return: pandas data-frame.
    '''
    df = pandas.read_csv(f)
    df = df.set_index('Unnamed: 0')  # PWMs are under an unnamed column
    return df


def get_ids(df):
    '''
    Gets all the unique PWM IDs since many PWMs exist with multiple variants.

    @param df: Parsed pandas data-frame.
    @return: list of PWM IDs.
    '''
    ids = [a_id for a_id in df.index]
    return sorted(set(ids))  # return sorted list of PWM IDs


def aggregate(df):
    '''
    Identifies the maximum score assigned to each PWM variant and merges their
    scores into one centralized PWM vector. This vector is simply the PWM
    ID that is shared by all its variants. In doing so, the number of
    observations decreases since all variants have been aggregated into one
    entry. Output from each vector is saved locally as a CSV file.

    @param df: Parsed pandas data-frame
    '''
    ids = get_ids(df)
    dfs = []
    for i in ids:
        int_locations = []
        for linenum, idx in enumerate(df.index):
            if str(idx).endswith(i):
                int_locations.append(linenum)
        df_subset = df.loc[int_locations]
        df_max = pandas.DataFrame(df_subset.max(axis=0))
        df_max.columns = [i]
        df_max = df_max.transpose()
        dfs.append(df_max)
    df_concat = pandas.concat(dfs)
    df_concat.to_csv('aggregated.csv', header=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', metavar='FILE', required=True,
                        help='Report file generated using LASSO [req]')
    args = vars(parser.parse_args())
    df = parse_report(f=args['f'])
    aggregate(df)
