#!/usr/bin/env python
import numpy as np
import sys
from astropy.io import fits
import re
import astropy.stats
import scipy
import scipy.ndimage as ndimage
from jbastro.misc import rangify
from m2fs.seqno import get_seqnos
from jbastro.astroLib import crreject
import argparse
import m2fs.obs
import os
from astropy.time import Time, TimeDelta
from m2fs.obs.stack import *

def parse_cl():
    parser = argparse.ArgumentParser(description='M2FS File Stacker',
                                     add_help=True)
    parser.add_argument('-d', '--dir', dest='dir',
                        action='store', required=False, type=str,
                        help='source dir for files', default='./')
    parser.add_argument('-l', '--listfile', dest='listfile', default='',
                        action='store', required=False, type=str,
                        help='testfile with files to proces')
    parser.add_argument('--crreject', dest='do_cosmic', default=False,
                        action='store_true', required=False,
                        help='Do cosmic ray rejection')
    parser.add_argument('-o', dest='outdir', default='./',
                        action='store', required=False, type=str,
                        help='Out directory')
    parser.add_argument('-z', dest='gzip', default=False,
                        action='store_true', required=False,
                        help='gzip fits files')
    parser.add_argument('--overwrite', dest='overwrite', default=False,
                        action='store_true', required=False,
                        help='Overwite existing files')
    parser.add_argument('-t', dest='dry_run', default=False,
                        action='store_true', required=False,
                        help='Test, dont do anything')

    args = parser.parse_args()
    if args.outdir[-1] != os.path.sep:
        args.outdir += os.path.sep
    if args.dir[-1] != os.path.sep:
        args.dir += os.path.sep
    return args

if __name__ == '__main__':

    args = parse_cl()

    files = [os.path.join(dirpath, f) for dirpath, dirnames, files in os.walk(args.dir)
             for f in files if '.fits' in f and '-' not in f and
             ',' not in f and 'c' not in f and not f.startswith('.')]
    try:
        seqno_stacks = get_seqnos(args.listfile)

        to_stack_lists = [[f for f in files
                           if m2fs.obs.info(f, no_load=True).seqno_match(seqnos)
                           and side == m2fs.obs.info(f, no_load=True).side]
                           for seqnos in seqno_stacks for side in 'rb']

        to_stack_lists = [l for l in to_stack_lists if len(l) > 1]  # only stack 2 or more

    except IOError:
        print('No listfile')
        exit()

    for i, filestack in enumerate(to_stack_lists):
        print("Stacking {} of {} ({} files)".format(i, len(to_stack_lists), len(filestack)))
        if len(filestack)>5:
            print(" First file: "+filestack[0])
        if len(filestack) == 1:
            print("Updating a single frame")
            outf = os.path.join(args.outdir, os.path.basename(filestack[0]).replace('.fits', ''))
            patch_single(filestack[0], outf, gzip=args.gzip, overwrite=args.overwrite)
            continue
        seqnos = list(map(lambda f: int(m2fs.obs.info(f, no_load=True).seqno), filestack))
        filestack = list(zip(*sorted(list(zip(filestack, seqnos)), key=lambda x: x[1])))[0]
        color = m2fs.obs.info(filestack[0], no_load=True).side
        try:
            if (not args.overwrite and (os.path.exists(args.outdir + color + rangify(seqnos, delim='_') + '.fits') or
                                        os.path.exists(args.outdir + color + rangify(seqnos, delim='_') + '.fits.gz'))):
                continue
            print('   Stacking ', color + rangify(seqnos))
            if not args.dry_run:
                stackimage(filestack, args.outdir + color + rangify(seqnos, delim='_'),
                           gzip=args.gzip, do_cosmic=args.do_cosmic, overwrite=args.overwrite)

        except IOError as e:
            print("Couldn't stack {}".format(str(e)))
