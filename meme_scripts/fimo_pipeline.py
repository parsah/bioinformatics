'''
A pipeline to execute the MEME FIMO application. This script takes
a list of PFM files, converts each to MEME format, maps each onto
user-provided sequences and parses the output txt to yield an
abundance matrix.
'''

import argparse
import os
import subprocess
import shutil

class IterativeFIMOWorker():
    def __init__(self, pfm_dir, fasta, outdir):
        self.pfm_dir = pfm_dir
        self.fasta = fasta
        self.outdir = outdir
        self.pfms = [] # references user-provided PFM filenames
        
    def assign_pfms(self):
        '''
        Get a list of all PFMs within a user-provided directory
        @return: list of sorted fully-qualified PFMs
        '''
        files = os.listdir(path=self.pfm_dir)
        pfms = [self.pfm_dir + '/' + i for i in files if i.endswith('.pfm')]
        print('# PFMs found:', len(pfms))
        self.pfms = sorted(pfms)
        
    def rmdir(self):
        os.rmdir(self.outdir)
    
    def pfm_to_meme(self, tempdir, meme):
        ''' 
        Each PFM is stored in a temporary directory. Thus, running jaspar2meme
        converts this PFM and writes it as a MEME output file
        '''
        out = subprocess.Popen(["jaspar2meme", "-pfm" ,tempdir], 
                                   stdout=subprocess.PIPE).communicate()[0]
        meme_handle = open(tempdir + '/' + meme, 'w')
        meme_handle.write(str(out))
        meme_handle.flush()
        meme_handle.close()
    
    def run(self):
        print('Aligning individual PFMs...')
        for i, pfm in enumerate(self.pfms):
            progress = '[ ' + str(i+1) + ' / ' + str(len(self.pfms)) + ' ]'
            basename = os.path.basename(pfm)
            print(progress, basename)
            tempdir = self.outdir + '/temp/'
            meme_name = 'matrix.meme'
            os.makedirs(tempdir) # temporary directory to run jaspar2meme
            shutil.copy(pfm, tempdir) # move PFM -> temporary directory
            self.pfm_to_meme(tempdir=tempdir, meme=meme_name) # PFM to MEME
            
            shutil.rmtree(tempdir) # at the end, remove the folder

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pfms', required=True, metavar='DIR',
                            help='Folder containing PFMs [na]')
    parser.add_argument('-i', required=True, metavar='FILE',
                            help='Input FASTA file [na]')
    parser.add_argument('-o', required=True, metavar='FILE',
                            help='Output folder [na]')
    args = vars(parser.parse_args())
    ifw = IterativeFIMOWorker(pfm_dir = args['pfms'], fasta = args['i'], 
                              outdir = args['o'])
    ifw.assign_pfms() # set PFMs
    ifw.run()