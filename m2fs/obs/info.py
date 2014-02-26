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
    def __init__(self):
        self.side=''
        self.seqno=0
        self.header=None
        self.exp_midpoint=None

def info(file):



    ret=M2FS_Obs_Info()
    
    if 'r' in os.path.basename(file):
        ret.side='r'
    else:
        ret.side='b'
    
    ret.seqno=(os.path.basename(file)[1:5])
    ret.header=fits.getheader(file)
    ret.exp_midpoint=_compute_midpoint(ret.header)

    return ret


    