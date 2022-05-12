#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.coordinates import Latitude, Longitude
from astropy import units as u
import os.path
from jbastro.misc import derangify


def _compute_midpoint(head):
    try:
        time, offset=head['UT-MID'], TimeDelta(0.0,format='sec')
    except KeyError:
        time=head['UT-DATE']+' '+head['UT-TIME']
        offset=TimeDelta(head['EXPTIME']/2.0, format='sec')

    try:
        return Time(time, format='iso', scale='utc', lat=Latitude(-29.01423,unit=u.degree),
                    lon=Longitude(-70.69242,unit=u.degree))+offset
    except:
        return Time(time, format='iso', scale='utc', location=(Longitude(-70.69242,unit=u.degree),
                                                               Latitude(-29.01423,unit=u.degree)))+offset


class M2FS_Obs_Info(object):
    def __init__(self, file, no_load=False):
        self._file=file
        self.side='r' if 'r' in os.path.basename(file) else 'b'
        if no_load:
            self.header=None
        else:
             self.load()

    @property
    def file(self):
        return self._file

    @property
    def seqno(self):
        bf=os.path.basename(self.file)
        try:
            no=int(bf[1:].split('_')[0].split(',')[0].split('-')[0].split('.')[0].split('c')[0])
        except Exception as e:
            print('Fault on ' + self.file)
            raise e
        return '{:04}'.format(no)

    @property
    def seqnos(self):
        bf=os.path.basename(self.file)
        return derangify(bf[1:].split('.')[0].split('c')[0].replace('_',','))

    def load(self):
        try:
            self.header=fits.getheader(self.file)
        except IOError:
            self.header=None

    @property
    def midpoint(self):
        if self.header is None:
            self.load()
        if self.header is None:
            return None
        return _compute_midpoint(self.header)

    def seqno_match(self, seq_id):
        if type(seq_id) in [list, tuple]:
            ret=np.array(list(map(self.seqno_match, seq_id)))
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
        except TypeError as e:
            raise e

        return (not side or side==self.side) and no==int(self.seqno)

    @property
    def is_quadrant(self):
        """Relies on file name, doesn't inspect image"""
        return 'c' in os.path.basename(self.file)


def info(file, no_load=False):
    return M2FS_Obs_Info(file, no_load=no_load)
