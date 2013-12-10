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
            try:
                if l[0] in '1234567890':
                    ret.update(derangify(l.split()[0]))
            except Exception:
                raise
    return list(ret)

if __name__ =='__main__':
    
    args=parse_cl()
    
    files = [os.path.join(dirpath, f)
                 for dirpath, dirnames, files in os.walk(args.dir)
                 for f in files if 'c1.fits' in f] #allow fits.gz
    
    try:
        seqno=get_seqnos(args.listfile)
        files=[f for f in files if int(m2fs.obs.info(f)['seqno']) in seqno]
    except IOError:
        print 'No listfile, doing all'


    for i,f in enumerate(files):
        print "Merging {} of {}".format(i,len(files))
        mergequad(f, do_cosmic=args.do_cosmic, file=True, odir=args.outdir)
