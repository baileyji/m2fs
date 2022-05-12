#!/usr/bin/env python
from datetime import datetime
from multiprocessing.pool import ThreadPool
import multiprocessing as mp
import requests
import argparse
import os, sys, shutil
from glob import glob

from logging import getLogger
import logging

import m2fs.obs
from m2fs.obs.merge import mergequad
from m2fs.seqno import *

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



def compute_download_paths(seqnum_info, outbase='./', url='https://m2fs.astro.lsa.umich.edu/data'):
    """ files is list of ('utYYYMMDD', {rb}####) """
    # Dates extend before and after runs by a good bit
    pathmap = {'Aug2013': (datetime(2013, 8, 1), datetime(2013, 8, 31)),
               'Nov2013': (datetime(2013, 11, 1), datetime(2013, 12, 2)),
               'Feb2014': (datetime(2014, 2, 1), datetime(2014, 3, 1)),
               'Jun2014': (datetime(2014, 6, 1), datetime(2014, 6, 30)),
               'Sep2014': (datetime(2014, 9, 1), datetime(2014, 9, 30)),
               'Dec2014': (datetime(2014, 12, 1), datetime(2014, 12, 31)),
               'Feb2015': (datetime(2015, 2, 1), datetime(2015, 3, 10)),
               'Apr2015': (datetime(2015, 4, 1), datetime(2015, 4, 28)),
               'Jul2015': (datetime(2015, 7, 1), datetime(2015, 7, 30)),
               'Sep2015': (datetime(2015, 9, 1), datetime(2015, 9, 30)),
               'Nov2015': (datetime(2015, 11, 1), datetime(2015, 11, 30)),
               'Feb2016': (datetime(2016, 2, 1), datetime(2016, 2, 28)),
               'Jun2016': (datetime(2016, 6, 1), datetime(2016, 7, 2)),
               'AugSep2016': (datetime(2016, 8, 1), datetime(2016, 9, 15)),
               'NovDec2016': (datetime(2016, 11, 1), datetime(2016, 12, 15)),
               'FebMar2017': (datetime(2017, 2, 1), datetime(2017, 3, 10)),
               'MayJun2017': (datetime(2017, 5, 1), datetime(2017, 6, 10)),
               'Jul2017': (datetime(2017, 7, 1), datetime(2017, 7, 31)),
               'Sep2017': (datetime(2017, 9, 1), datetime(2017, 9, 30)),
               'Nov2017': (datetime(2017, 11, 1), datetime(2017, 11, 30)),
               'Feb2018': (datetime(2018, 2, 1), datetime(2018, 2, 28)),
               'May2018': (datetime(2018, 5, 1), datetime(2018, 5, 31)),
               'Aug2018': (datetime(2018, 8, 1), datetime(2018, 8, 31)),
               'NovDec2018': (datetime(2018, 11, 1), datetime(2018, 12, 10)),
               'FebMar2019': (datetime(2019, 2, 1), datetime(2019, 3, 10)),
               'MayJun2019': (datetime(2019, 5, 1), datetime(2019, 6, 10)),
               'AugSep2019': (datetime(2019, 8, 1), datetime(2019, 9, 10))}

    to_fetch = []
    for seq in seqnum_info:
        ut, seqstr, _ = seq
        run = ''
        for k, v in pathmap.items():
            if v[0] < datetime.strptime(ut, 'ut%Y%m%d') < v[1]:
                run = k
                break
        if not run:
            getLogger(__name__).warning("Unable to determine download location for {} with ut '{}'".format(seqstr, ut))
            continue
        to_fetch.extend([('{}/{}/{}/{}c{}.fits'.format(url, run, ut if run != 'AugSep2016' else 'ALL', seqstr, i),
                          os.path.join(outbase, run, ut, '{}c{}.fits'.format(seqstr, i)))
                         for i in range(1, 5)])

    return to_fetch


def fetch_url(entry):
    uri, path = entry

    r = requests.get(uri, stream=True, auth=('m2fs', 'm2fsdata'))
    if r.status_code != 200:
        getLogger(__name__).error('m2fs.astro.lsa.umich.edu returned {} for {}, skipping'.format(r.status_code, uri))
        return

    total_size = int(r.headers.get('content-length', 0))  # Total size in bytes.
    if os.path.exists(path):
        if os.stat(path).st_size == total_size:
            # getLogger(__name__).info('{} exists with same size, skipping'.format(path))
            return
        else:
            getLogger(__name__).info('{} exists with mismatched size, downloading'.format(path))

    getLogger(__name__).debug('Fetching {} to {}'.format(uri, path))

    directory = os.path.dirname(path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(path, 'wb') as f:
        for data in r.iter_content(1024):
            f.write(data)


if __name__ == '__main__':

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(message)s"))
    logger.addHandler(handler)

    args = parse_cl()

    if args.useweb:
        getLogger(__name__).info('Downloading files from m2fs.astro.lsa.umich.edu')
        seqnos = sorted(get_seqnos_advanced(args.listfile))
        paths = compute_download_paths(seqnos, outbase=args.dir)
        if args.useweb>1:
            pool = mp.pool.Pool(args.useweb)
            pool.map(fetch_url, paths)
            pool.close()
            pool.join()
        else:
            for p in paths:
                fetch_url(p)
        getLogger(__name__).info('Downloads complete, run again without --dl to continue')
        exit()

    getLogger(__name__).info('Looking in {} for files'.format(args.dir))
    basedirs = glob(args.dir) if '*' in args.dir else (args.dir,)
    files = [os.path.join(dirpath, f)
             for basedir in basedirs
             for dirpath, dirnames, files in os.walk(basedir)
             for f in files
             if 'c1.fits' in f and f not in ('rc1.fits', 'bc1.fits') and not f.startswith('.')]

    try:
        seqno_data = get_seqnos_advanced(args.listfile)
        # seqno = [s[0][1] for s in seqno_data]
        # seqno=get_seqnos(args.listfile)
        info = [m2fs.obs.info(f, no_load=True) for f in files]

        files, cr_requested = [], []
        for i in info:
            for date, s, docr in seqno_data:
                if i.seqno_match(s):
                    files.append(i.file)
                    cr_requested.append(docr)
                    break
    except IOError:
        print('No listfile, doing all')

    if len(files) == 0:
        print('No files found')

    if args.copy:
        print('Copying:')
        for f in files:
            print(f)
        for f in files:
            d = os.path.join(args.copy, os.path.split(os.path.dirname(f))[1])
            if not os.path.exists(d):
                os.mkdir(d)
            if not os.path.exists(os.path.join(d, os.path.basename(f))):
                print('Copying {} to {}'.format(f, d))
                shutil.copy2(f, d)
        exit()

    if not os.path.exists(os.path.dirname(args.outdir)):
        os.makedirs(os.path.dirname(args.outdir))

    for i, (f, docr) in enumerate(zip(files, cr_requested)):
        print("Merging {} of {}, with{} CR rejection".format(i, len(files), '' if args.do_cosmic or docr else 'out'))
        s = os.statvfs('/')
        if (s.f_bavail * s.f_frsize) / 1024 ** 2 < 5000:
            print("Disk space too low")
            exit()

        cosmic_settings = False
        if args.do_cosmic or docr:
            cosmic_settings = {'sigclip': args.sigclip, 'sigfrac': args.sigfrac,
                               'objlim': args.objlim, 'iter': args.criter}

        mergequad(f, do_cosmic=cosmic_settings, file=True, odir=args.outdir,
                  dogzip=args.gzip, overwrite=args.overwrite)
