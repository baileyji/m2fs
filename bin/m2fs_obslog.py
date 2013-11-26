#!/usr/bin/env python2.7
from m2fs.obs import summarize

if __name__ =='__main__':
    import glob
    files=glob.glob('*.fits')
    summarize.make_table(files)
