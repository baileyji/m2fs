#!/usr/bin/env python
from m2fs.obs.merge import mergequad
import argparse
from jbastro.misc import derangify
import m2fs.obs
import os

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
    
    files = [os.path.join(dirpath, f)
             for dirpath, dirnames, files in os.walk(args.dir)
             for f in files
             if 'c1.fits' in f and f not in ['rc1.fits','bc1.fits']] #allow fits.gz

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


    for i,f in enumerate(files):
        print "Merging {} of {}".format(i,len(files))
        s=os.statvfs('/')
        if (s.f_bavail * s.f_frsize) / 1024**2 < 5000:
            print "Disk space too low"
            continue
        if args.do_cosmic:
            cosmic_settings={'sigclip': 15.0, 'sigfrac': .3, 'objlim': 1.4,
                'iter':10}
        else:
            cosmic_settings=False
        mergequad(f, do_cosmic=cosmic_settings, file=True, odir=args.outdir)
