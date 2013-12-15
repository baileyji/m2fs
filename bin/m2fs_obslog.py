#!/usr/bin/env python2.7
from m2fs.obs import summarize
import os
import sys
import argparse

def parse_cl():
    parser = argparse.ArgumentParser(description='Observing Log generator',
                                     add_help=True)
    parser.add_argument('-d','--dir', dest='dir',
                        action='store', required=False, type=str,
                        help='source dir for files',default='./')
    return parser.parse_args()

#"<html><body>\n"+'\n'.join(rows)+"</html></body>\n"
#header_row="<b>"+column_fmtstr.format(column_names)+"</b>"
#rows=[]
#for i, r in enumerate(recs):
#    if i % 40 ==0:
#        rows.append(header_row)
#    rows.append(r)


def filter_files(files):
    """Only show merged frame or first quadrant (if merged not available)"""
    out=[]
    for f in files:
        fname=os.path.basename(f)
        try:
            if ('c' not in fname or
                ('c1' in fname and f[0:-7]+'.fits' not in files and
                f[0:-7]+'.fits.gz' not in files)):
                out.append(f)
        except IndexError:
                out.append(f)
    return out

def hjobs(header):
    return True
#    return (header['FILTER']=='HotJupiter' or
#            (int(header['FILENAME'][1:-2]) in range(106,116)))

if __name__ =='__main__':
    args=parse_cl()

    files = [os.path.join(dirpath, f)
             for dirpath, dirnames, files in os.walk(args.dir)
             for f in files if f.endswith('.fits') or f.endswith('.fits.gz')]

    files=filter_files(files)

    summarize.run(files,
              prepend=['NIGHT', 'RA','DEC'],
              append=['FF-THNE','FF-THAR','FF-QRTZ'],
              filterfunc=hjobs)
