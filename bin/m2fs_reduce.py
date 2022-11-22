from datetime import datetime
from multiprocessing.pool import ThreadPool
import multiprocessing as mp
import argparse
import os, sys, shutil
from glob import glob
from m2fs.obs.fetcher import compute_download_paths, fetch_url

from logging import getLogger
import logging

import m2fs.obs
from m2fs.obs.merge import mergequad
from m2fs.seqno import *
from m2fs.obs.stack import update_single, stack_images


# M2FSIFUMListfileDataset('/Users/one/no_backup/ifum/list_sept', rawdir, ddir)
msf, seq_strs = generate_filesets('/Users/one/Desktop/m2fsnov22/list.txt',
                                  '/Volumes/DATA_CLAY/M2FS/',
                                  '/Users/one/Desktop/m2fsnov22/merged/',
                                  '/Users/one/Desktop/m2fsnov22/stacked/')


# files = [mf for a in x for mf in find(a, msf)]

for collection in msf:
    for f in collection:
        if f.listfile_record.bias or f.listfile_record.dark:
            continue
        f.merged()



dry_run=False #args.dry_run
overwrite = False
force_gzip = True
do_cosmic_stack = False
for i, collection in enumerate(msf):
    print(f"Processing {i} of {len(msf)}")
    for side, group in collection.sides:

        files = [x.merged_file() for x in group]
        print(f"  Side {side}, {len(files)} files")
        of = collection.stacked_file(side)

        if len(files) == 1:
            print(f"Updating {files[0]}")
            update_single(files[0], of, overwrite=overwrite)
        else:
            try:
                print(f'   Stacking {of}')
                if dry_run:
                    continue
                stack_images(files, of, do_cosmic=do_cosmic_stack, overwrite=overwrite)

            except IOError as e:
                print("Couldn't stack {}".format(str(e)))




def parse_cl():
    parser = argparse.ArgumentParser(description='Quadrant merger',
                                     add_help=True)
    parser.add_argument('--dl', dest='useweb', type=int, default=0,
                        action='store', required=False,
                        help='Download files from webserver with N threads')
    parser.add_argument('-d', '--dir', dest='dir',
                        action='store', required=False, type=str,
                        help='source dir for files', default='./')
    parser.add_argument('-l', '--listfile', dest='listfile', default='',
                        action='store', required=False, type=str,
                        help='testfile with files to proces')
    parser.add_argument('--crreject', dest='do_cosmic', default=False,
                        action='store_true', required=False,
                        help='Do cosmic ray rejection')
    parser.add_argument('-o', dest='outdir', default='',
                        action='store', required=False, type=str,
                        help='Out directory')
    parser.add_argument('--sigclip', dest='sigclip', default=15.0,
                        action='store', required=False, type=float,
                        help='CR Sigmaclip (sigma limit for flagging as CR)')
    parser.add_argument('--sigfrac', dest='sigfrac', default=0.3,
                        action='store', required=False, type=float,
                        help='CR Sigmafrac (sigclip fract for neighboring pix)')
    parser.add_argument('--objlim', dest='objlim', default=1.4,
                        action='store', required=False, type=float,
                        help='CR Object Limit (raise if normal data clipped)')
    parser.add_argument('--criter', dest='criter', default=5,
                        action='store', required=False, type=int,
                        help='CR Iteration Limit')
    parser.add_argument('--overwrite', dest='overwrite', default=False,
                        action='store_true', required=False,
                        help='overwrite existing output')
    parser.add_argument('--copy', dest='copy', default='',
                        action='store', required=False,
                        help='copy found files and exit (for downloading from sshfs mounts')
    parser.add_argument('-z', dest='gzip', default=False,
                        action='store_true', required=False,
                        help='gzip fits files')
    return parser.parse_args()
