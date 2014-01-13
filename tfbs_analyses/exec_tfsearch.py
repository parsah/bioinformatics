'''
Useful script for running tfSearch on a large input FASTA file.
Each FASTA entry in the user-provided file is therefore executed against
a user-provided PWM file using the tfSearch application.
'''

import sys
import subprocess
import tempfile
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

def run_script(f, p):
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
        if len(sys.argv) < 2:
            raise OSError('Input FASTA and TRANSFAC PWMs are required [error]')
        else:
            fasta = sys.argv[1] # user-provided input FASTA file
            pwms = sys.argv[2] # user-provided input PWMs (TRANSFAC format)
            run_script(fasta, pwms) # run the program
    except OSError as e:
        print(e)
    except KeyboardInterrupt as e:
        print()

