#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from astropy.time import Time
import os.path

def _compute_midpoint(head):
    try:
        return Time(head['UT-MID'], format='iso', scale='utc')
    except KeyError:
        #import ipdb;ipdb.set_trace()
        start=Time(head['UT-DATE']+' '+head['UT-TIME'],
                   format='iso', scale='utc')
        end=Time(head['UT-DATE']+' '+head['UT-END'],
                   format='iso', scale='utc')

        return start +.5*(end-start)

class M2FS_Obs_Info(object):
    def __init__(self,file, no_load=False):
    
        self.file=file
        if 'r' in os.path.basename(file):
            self.side='r'
        else:
            self.side='b'
        self.seqno=(os.path.basename(file)[1:5])
    
        if no_load:
            self.header=None
            self.midpoint=None
        else:
            self.load()

    def load(self):
        try:
            self.header=fits.getheader(self.file)
            self.midpoint=_compute_midpoint(self.header)
        except IOError:
            self.header=None
            self.midpoint=None

    def seqno_match(self, seq_id):
        if type(seq_id) in [list, tuple]:
            ret=np.array(map(self.seqno_match, seq_id))
            return ret.any()
        
        if type(seq_id)==str:
            if not seq_id:
                return False
            if seq_id[0] in 'RBrb':
                side=seq_id[0].lower()
                seq_id=seq_id[1:]
            else:
                side=''
        try:
            no=int(seq_id)
        except TypeError:
                raise

        return (not side or side==self.side) and no==int(self.seqno)

    @property
    def is_quadrant(self):
        """Relies on file name, doesn't inspect image"""
        return 'c' in os.path.basename(self.file)


def info(file, no_load=False):
    return M2FS_Obs_Info(file, no_load=no_load)


    