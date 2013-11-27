#!/usr/bin/env python
import m2fs.obs.merge import mergequad

if __name__ =='__main__':
   import glob
   files = glob.glob('*c?.fits')
   seqnos=set([ int(x[1:-7]) for x in files])
   for i in seqnos:
       print "Merging {} of {}".format(i,len(seqnos))
       mergequad(i)
