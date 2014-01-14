'''
Useful script for running tfSearch on a large input FASTA file.
Each FASTA entry in the user-provided file is therefore executed against
a user-provided PWM file using the tfSearch application.
'''

import subprocess
import tempfile
import argparse
from Bio import SeqIO

def _exec_tfsearch(f, p):
    ''' 
    Performs tfSearch on a one-entry FASTA file
    @param f: One-entry FASTA file.
    @param p: TRANSFAC PWMs.
    '''
    proc = subprocess.Popen(["tfSearch", p, f], stdout=subprocess.PIPE)
    stdout, err = proc.communicate()
    if not err:
        lines = stdout.splitlines()
        results = '' # results-string
        record = SeqIO.read(f, 'fasta')
        for line in lines:
            results += record.id + ' ' + str(len(record.seq)) + ' ' + line.decode('utf-8') + '\n'
        return results.strip()
    else:
        raise IOError('Error thrown during running of tfSearch [error]')

def run_tfsearch(f, p):
    ''' 
    Runs tfSearch given a FASTA file and set of PWMs
    @param f: Input FASTA file.
    @param p: Input TRANSFAC PWMs file.
    '''
    records = SeqIO.parse(f, 'fasta')
    for record in records: # iterate through each FASTA entry
        print('*** RUNNING ', record.id, '***')
        temp_file = tempfile.NamedTemporaryFile(delete=True, mode='w') # create temp file
        temp_file.write('>' + record.id + '\n' + str(record.seq))
        temp_file.flush()
        results = _exec_tfsearch(temp_file.name, p) # run tfSearch given FASTA file and PWMs 
        print(results + '\n') # write results to a centralized output file
        temp_file.close()
    print('All operations successfully completed [OK]')

if __name__ == '__main__':
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-fasta', required=True, metavar='FASTA',
                            help='Input FASTA file [req]')
        parser.add_argument('-pwms', required=True, metavar='FILE',
                            help='File containing PWMs [req]')
        args = vars(parser.parse_args())
        run_tfsearch(f = args['fasta'], p = args['pwms']) # run the program
    except OSError as e:
        print(e)
    except KeyboardInterrupt as e:
        print()
