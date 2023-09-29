import os.path
from datetime import datetime
import argparse
from logging import getLogger
import logging
from m2fs.seqno import *
from m2fs.obs.stack import patch_single, stack_images
import numpy as np


def parse_cl():
    parser = argparse.ArgumentParser(description='M2FS Reduction Script', add_help=True)
    parser.add_argument('--dl', dest='useweb', type=int, default=0,
                        action='store', required=False,
                        help='Download files from webserver with N threads')
    parser.add_argument('-r', dest='raw', default='/Volumes/DATA_CLAY/M2FS/',
                        action='store', required=False, type=str,
                        help='source dir for raw quadrant files')
    parser.add_argument('-l', '--listfile', dest='listfile', default='list',
                        action='store', required=False, type=str,
                        help='frame listing with frame numbers to process')
    parser.add_argument('--skipdark', dest='skip_dark', default=False,
                        action='store_true', required=False,
                        help='Skip processing dark frames')
    parser.add_argument('--skipbias', dest='skip_bias', default=False,
                        action='store_true', required=False,
                        help='Skip processing bias frames')
    parser.add_argument('--dry', dest='dry_run', default=False,
                        action='store_true', required=False,
                        help='Do not actually do anything')
    parser.add_argument('-o', dest='out', default='./no_backup/lcoreduction/reduced/',
                        action='store', required=False, type=str,
                        help='Out directory, also used for non-raw inputs')
    parser.add_argument('--overwrite', dest='overwrite', default=False,
                        action='store_true', required=False,
                        help='overwrite existing output')
    # parser.add_argument('--copy', dest='copy', default='',
    #                     action='store', required=False,
    #                     help='copy found files and exit (e.g. for downloading from sshfs mounts')
    parser.add_argument('-z', dest='gzip', default=False,
                        action='store_true', required=False,
                        help='Force gzip of fits file output')
    # parser.add_argument('--sigclip', dest='sigclip', default=15.0,
    #                     action='store', required=False, type=float,
    #                     help='CR Sigmaclip (sigma limit for flagging as CR)')
    # parser.add_argument('--sigfrac', dest='sigfrac', default=0.3,
    #                     action='store', required=False, type=float,
    #                     help='CR Sigmafrac (sigclip fract for neighboring pix)')
    # parser.add_argument('--objlim', dest='objlim', default=1.4,
    #                     action='store', required=False, type=float,
    #                     help='CR Object Limit (raise if normal data clipped)')
    # parser.add_argument('--criter', dest='criter', default=5,
    #                     action='store', required=False, type=int,
    #                     help='CR Iteration Limit')

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_cl()

    logging.basicConfig()
    log = getLogger('m2fs_reduce')
    log.propagate=True
    list_file = args.listfile
    raw_dir = args.raw
    out_dir = args.out

    assert os.path.exists(raw_dir)

    dry_run = False or args.dry_run
    overwrite = False or args.overwrite
    force_gzip = args.gzip
    do_cosmic_stack = False # do CR rejections
    do_flux_weighted_stack = True

    dset = M2FSIFUMListfileDataset(list_file, raw_dir, out_dir)
    msf = dset.msf
    seq_strs = dset.seq_strs

    skip_bias = True or args.skip_bias
    skip_dark = True or args.skip_dark
    # files = [mf for a in x for mf in find(a, msf)]

    # Merge raw quadrants, TODO throw more processes at this we can do CR on as many quadrants as we have cpus
    for collection in msf:
        for f in collection.files:
            if f.listfile_record.bias and skip_bias:
                continue
            if f.listfile_record.dark and skip_dark:
                continue
            f.merged()

    #Stack the merged images
    for i, collection in enumerate(msf):
        print(f"Processing {i} of {len(msf)}")

        for side, group in collection.sides:

            files = [x.merged_file() for x in group]
            print(f"  Side {side}, {len(files)} files")
            of = collection.stacked_file(side)

            if len(files) == 1:
                print(f"Updating {files[0]}")
                patch_single(files[0], of, overwrite=overwrite)
            else:
                try:
                    print(f'Stacking {of}')
                    if dry_run:
                        continue
                    stack_images(files, of, do_cosmic=do_cosmic_stack, overwrite=overwrite,
                                 flux_weighted=do_flux_weighted_stack)

                except IOError as e:
                    print("Couldn't stack {}".format(str(e)))
