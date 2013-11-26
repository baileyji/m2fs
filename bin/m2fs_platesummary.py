#!/usr/bin/python
from m2fs.plate import summarize
import sys, glob

target=summarize.target

extra=[
target('HD223311','23 49 14.1','-06 18 20'),
target('HIP48331','09 51 06.68','-43 30 05.9'),
target('HIP10798','02 18 58.65','-25 56 48.4'),
target('EG21','03 10 30.98','-68 36 02.2')
]


if __name__ == '__main__':
    
    if len(sys.argv) >2 :
        dir=sys.argv[1]
        fname=sys.argv[2]
    else:
        dir='./'
        fname=sys.argv[1]
    sfile=fname+'_summ.txt'
    tfile=fname+'_tlist.txt'
    files = glob.glob(dir+'*.plate')
    trec=summarize.write_summary_file(sfile, files)
    summarize.write_target_list(tfile, trec+extra)