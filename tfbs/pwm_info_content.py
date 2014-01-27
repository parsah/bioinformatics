''' 
Classes that parse and quantify PWMs used by the MEME Suite.
Such classes first parse a MEME PWM file, with each PWM being modeled as
an object of type MEMEPWMRecord. Each record of this type can have its
information content computed; useful in analyses where you wish to determine
which PWM has the lowest information given another. This very conclusion helps
to filter-out low-information PWMs, ultimately helping to reduce the 
feature-space during TFBS analysis.
'''

import argparse
import math
import numpy
import re

class MEMEPWMRecord():
    ''' 
    Models an individual PWM originally presented in the MEME format.
    Each object references its original weight matrix and a corresponding 
    motif name.
    '''
    def __init__(self, name):
        ''' 
        Build a skeleton PWM record.
        @param name: String referencing the respective object.
        '''
        self.matrix = []
        self.name = name
    
    def set_pwm(self, m):
        ''' 
        Set the actual matrix referencing this record object. This matrix
        must be originally from a MEME PWM file; see file format for details.
        @param m: Multi-dimensional matrix of probabilities.
        '''
        assert m.shape[1] == 4 # assert 4 rows (A, T, G, C) are added
        self.matrix = m
    
    def information(self):
        info = 0.0
        for row in range(len(self.matrix)):
            for col in range(len(self.matrix[row])):
                cell = self.matrix[row][col] # reference matrix cell
                if cell == 0:
                    cell = 1e-4
                info += cell * math.log(cell / 0.25)
        return -info     

class MEMEPWMParser():
    def __init__(self, f):
        self.f = f
        self.records = {}
        self.parse()
        
    def parse(self):
        entries = open(self.f)
        hit_str = '' # save all results to string for easier processing
        for line in entries:
            line = line.strip()
            # match entries beginning with MOTIF and floating points 
            match = re.findall('^MOTIF|^[0-9]', line)
            if match:
                hit_str += line + '\n'
        hit_str = hit_str.split('MOTIF')[1:] # first hit is empty
        
        # iterate over each motif and save into an array
        for i in hit_str:
            i = i.splitlines()
            motif_id = i[0].strip().split()[0]
            probs = i[1: ]
            matrix = []
            for row in probs:
                matrix.append([float(x) for x in row.split()])
            rec = MEMEPWMRecord(motif_id)
            rec.set_pwm(numpy.array(matrix)) # set the respective PWM
            self.records[rec.name] = rec.information()
    
    def get_records(self):
        return self.records        
