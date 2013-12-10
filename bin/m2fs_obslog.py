#!/usr/bin/env python2.7
from m2fs.obs import summarize
import os
import sys

def filter_files(files):
    """Only show merged frame or first quadrant (if merged not available)"""
    out=[]
    for f in files:
        try:
            if 'c' not in f or ('c1' in f and f[0:-7]+'.fits' not in files):
                out.append(f)
        except IndexError:
                out.append(f)
    return out

def hjobs(header):
    return (header['FILTER']=='HotJupiter' or
            (int(header['FILENAME'][1:-2]) in range(106,116)))

if __name__ =='__main__':
    if len(sys.argv) >1 :
        dir = sys.argv[1]
    else:
        dir ='./'

    files = [os.path.join(dirpath, f)
             for dirpath, dirnames, files in os.walk(dir)
             for f in files if f.endswith('.fits')]
             
    files=filter_files(files)
    summarize.run(files,
                  prepend=['NIGHT', 'RA','DEC'],
                  append=['FF-THNE','FF-THAR','FF-QRTZ'],
                  filterfunc=hjobs)
