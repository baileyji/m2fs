#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import argparse
from collections import namedtuple
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.coordinates import Latitude, Longitude
from astropy import units as u
import os.path
from jbastro.misc import derangify


def _compute_midpoint(head):
    try:
        time, offset = head['UT-MID'], TimeDelta(0.0, format='sec')
    except KeyError:
        time = head['UT-DATE'] + ' ' + head['UT-TIME']
        offset = TimeDelta(head['EXPTIME'] / 2.0, format='sec')

    try:
        return Time(time, format='iso', scale='utc', lat=Latitude(-29.01423, unit=u.degree),
                    lon=Longitude(-70.69242, unit=u.degree)) + offset
    except:
        return Time(time, format='iso', scale='utc', location=(Longitude(-70.69242, unit=u.degree),
                                                               Latitude(-29.01423, unit=u.degree))) + offset


class M2FSObsInfo(object):
    def __init__(self, file, no_load=False):
        self._file = file
        self.side = 'r' if 'r' in os.path.basename(file) else 'b'
        self.instrument = ''
        self.header = None
        if not no_load:
            self.load()

    @property
    def exists(self):
        return os.path.exists(self._file)

    @property
    def utdir(self):
        return os.path.split(self._file)[-1]

    @property
    def temps(self):
        x = [self.header.cards[k] for k in self.header if k.startswith('TEMP')]
        names = [c.comment.lower() for c in x]
        names[names.index('ccd temperature [c]')]='ccd'
        temps = namedtuple('{}Temps'.format(self.instrument.upper()), names)
        return temps(*[c.value for c in x])

    @property
    def file(self):
        return self._file

    @property
    def seqno(self):
        """ returns -1 if unable to parse"""
        try:
            bf = os.path.basename(self.file)
            try:
                no = int(bf[1:].split('_')[0].split(',')[0].split('-')[0].split('.')[0].split('c')[0])
            except Exception as e:
                print('Fault on ' + self.file)
                raise e
            return '{:04}'.format(no)
        except Exception:
            return -1

    @property
    def seqnos(self):
        bf = os.path.basename(self.file)
        return derangify(bf[1:].split('.')[0].split('c')[0].replace('_', ','))

    def load(self):
        try:
            self.header = fits.getheader(self.file)
            self.instrument = self.header['INSTRUME'].lower()
        except IOError:
            self.header = None

    @property
    def midpoint(self):
        if self.header is None:
            self.load()
        if self.header is None:
            return None
        return _compute_midpoint(self.header)

    def seqno_match(self, seq_id):
        if type(seq_id) in [list, tuple]:
            ret = np.array(list(map(self.seqno_match, seq_id)))
            return ret.any()

        if type(seq_id) == str:
            if not seq_id:
                return False
            if seq_id[0] in 'RBrb':
                side = seq_id[0].lower()
                seq_id = seq_id[1:]
            else:
                side = ''
        try:
            no = int(seq_id)
        except TypeError as e:
            raise e

        return (not side or side == self.side) and no == int(self.seqno)

    @property
    def is_quadrant(self):
        """Relies on file name, doesn't inspect image"""
        return 'c' in os.path.basename(self.file)

    @property
    def binning(self):
        b = self.header['BINNING']
        return int(b[0]), int(b[-1])

    @property
    def exposure_id(self):
        """ut dir and file up to first ."""
        return self.utdir, os.bath.basename(self._file).lpartition('.')[0]


M2FS_Obs_Info = M2FSObsInfo


def info(file, no_load=False):
    return M2FSObsInfo(file, no_load=no_load)
