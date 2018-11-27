#!/usr/bin/env python
from m2fs.obs.merge import mergequad
import argparse
from jbastro.misc import derangify
import m2fs.obs
import os
from glob import glob

def parse_cl():
    parser = argparse.ArgumentParser(description='Quadrant merger',
                                     add_help=True)
    parser.add_argument('-d','--dir', dest='dir',
                        action='store', required=False, type=str,
                        help='source dir for files',default='./')
    parser.add_argument('-l','--listfile', dest='listfile', default='',
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
    parser.add_argument('--criter', dest='criter', default=10,
                        action='store', required=False, type=int,
                        help='CR Iteration Limit')
    parser.add_argument('--overwrite', dest='overwrite', default=False,
                        action='store_true', required=False,
                        help='clobber existing output')
    parser.add_argument('-z', dest='gzip', default=False,
                 action='store_true', required=False,
                 help='gzip fits files')
    return parser.parse_args()


def get_seqnos(listfile):
    ret=set()
    with open(listfile,'r') as lf:
        for l in lf:
            if l[0] in '1234567890':
                range=l.split()[0]
                ret.update(map(str,derangify(range)))
            elif len(l)>1 and l[0] in 'RBrb' and l[1] in '1234567890':
                range=l[1:].split()[0]
                ret.update(map(lambda x: l[0].lower()+str(x), derangify(range)))
    return list(ret)

if __name__ =='__main__':
    
    args=parse_cl()
    
    print 'Looking in {} for files'.format(args.dir)
    
    if '*' in args.dir:
        basedirs = glob(args.dir)
    else:
        basedirs = (args.dir,)
        
    files = [os.path.join(dirpath, f)
             for basedir in basedirs
             for dirpath, dirnames, files in os.walk(basedir)
             for f in files
             if 'c1.fits' in f and f not in ['rc1.fits','bc1.fits']]
#    files = [os.path.join(dirpath, f)
#             for dirpath, dirnames, files in os.walk(basedir) if base in dirpath
#             for f in files
#             if 'c1.fits' in f and f not in ['rc1.fits','bc1.fits']]

    try:
        seqno=get_seqnos(args.listfile)
        info=[m2fs.obs.info(f,no_load=True) for f in files]

        files=[]
        for i in info:
            for s in seqno:
                if i.seqno_match(s):
                    files.append(i.file)
                    break
    except IOError:
        print 'No listfile, doing all'

    if len(files) == 0:
        print 'No files found'

    for i,f in enumerate(files):
        print "Merging {} of {}".format(i,len(files))
        s=os.statvfs('/')
        if (s.f_bavail * s.f_frsize) / 1024**2 < 5000:
            print "Disk space too low"
            continue
        if args.do_cosmic:
            cosmic_settings={'sigclip': args.sigclip, 'sigfrac': args.sigfrac,
                             'objlim': args.objlim, 'iter':args.criter}
        else:
            cosmic_settings=False
        mergequad(f, do_cosmic=cosmic_settings, file=True, odir=args.outdir,
                  dogzip=args.gzip, clobber=args.overwrite)
