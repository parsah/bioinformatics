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
    pwm_contents = []
    for line in open(fname):
        line = line.strip().split(' ') # PWMs are space-delimited
        # multiply floats by 100 to represent as percentage.
        line = [round(float(i)*100, 2) for i in line]
        pwm_contents.append(line)
    return pwm_contents

def cleanse(pwm_data):
    ''' 
    Transposes a 2D matrix.
    @param pwm_data: Multi-dimensional matrix of real-numbers
    @return: transposed matrix.
    '''
    tranposed = list(zip(*pwm_data))
    row_names = ['A', 'C', 'G', 'T']
    for num, row in enumerate(tranposed):
        row = [str(i) for i in row] # cast float back to string for writing
        print(row_names[num]+ ' ' + ' '.join(row))

if __name__ == '__main__':
    desc = """ 
    Adjusts a PWM so that minimum-cell values are 1 and maximum values are
    set to 100. This ensures a PWM with fractions are scaled for readability.
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-pwm', required=True, metavar='PWM',
                        help='PWM file; 4x columns wide [na]')
    args = vars(parser.parse_args())
    pwm_data = parse(fname=args['pwm'])
    pwm_data = cleanse(pwm_data) # streams data out to standard-out
    