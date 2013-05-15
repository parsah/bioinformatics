'''
A helpful script which takes a PWM and adjusts its cell-values so that the
minimum is 1 and the maximum is 100.
'''

import argparse

def parse(fname):
    ''' 
    Parses a user-provided PWM input file.
    @param fname: Input filename referencing the PWM.
    @return: multi-dimensional array containing PWM values.
    '''
    for line in open(fname):
        line = line.strip().split(' ')
        print(line)

if __name__ == '__main__':
    desc = """ 
    Adjusts a PWM so that minimum-cell values are 1 and maximum values are
    set to 100. This ensures a PWM with fractions are scaled for readability.
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-pwm', required=True, metavar='PWM',
                        help='PWM file; 4x columns wide [na]')
    args = vars(parser.parse_args())
    parse(fname=args['pwm'])
    