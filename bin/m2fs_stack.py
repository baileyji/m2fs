#!/usr/bin/env python2.7
import numpy as np
from astropy.io import fits
import re
import pylab as pl
import astropy.stats
import scipy.ndimage as ndimage
from jbastro.misc import rangify, derangify
from jbastro.astroLib import crreject
import argparse
import m2fs.obs
import os
import operator

def parse_cl():
    parser = argparse.ArgumentParser(description='Quadrant merger',
                                     add_help=True)
    parser.add_argument('-d','--dir', dest='dir',
                     action='store', required=False, type=str,
                     help='source dir for files',default='./')
    parser.add_argument('-l','--listfile', dest='listfile', default='',
                     action='store', required=True, type=str,
                     help='testfile with files to proces')
    parser.add_argument('--crreject', dest='do_cosmic', default=False,
                     action='store_true', required=False,
                     help='Do cosmic ray rejection')
    parser.add_argument('-o', dest='outdir', default='./',
                     action='store', required=False, type=str,
                     help='Out directory')
    
    args=parser.parse_args()
    if args.outdir[-1]!=os.path.sep:
        args.outdir+=os.path.sep
    if args.dir[-1]!=os.path.sep:
        args.dir+=os.path.sep
    return args


def get_seqlists(listfile):
    ret=set()
    with open(listfile,'r') as lf:
        for l in lf:
            try:
                if l[0] in '1234567890':
                    ret.add(derangify(l.split()[0]))
            except Exception:
                raise
    return list(ret)




def stackimage(files,  outfile, do_cosmic=False, **crparams):
    """List of files to stack, output file sans extension."""

    
    nfile=len(files)
    #load first bias to get info
    with fits.open(files[0]) as im:
        header=im[0].header
        imcube=np.zeros(im[0].data.shape+(nfile,),dtype=np.float32) #15 frames with float64 would be 3.7GB ram
        if len(im) >2:
            masked=True
            mask=np.zeros_like(imcube, dtype=np.bool)
        else:
            masked=False

    varcube=np.zeros_like(imcube)

    #load data and variance frames
    members=[]
    etimes=[]
    try:
        for i,f in enumerate(files):
            with fits.open(f) as im:
                etimes.append(float(im[0].header['EXPTIME']))
                members.append(im[0].header['FILENAME'])
                if etimes[-1]!=0.0:
                    imcube[:,:,i]=im[0].data/etimes[-1]
                else:
                    imcube[:,:,i]=im[0].data
                varcube[:,:,i]=im[1].data

                if masked:
                    msk=im[2].data
                    if do_cosmic:
                        foo=crreject(im[1].data)
                        msk+=foo
                    msk=ndimage.morphology.binary_dilation(
                        msk.astype(np.bool), structure=np.ones((3,3)),
                        iterations=1, mask=None, output=None,
                        border_value=0, origin=0, brute_force=False)
                
                    mask[:,:,i]=msk
    except ValueError as e:
        print "ValueError while merging:", files,f
        return

    exptime=sum(etimes)

    if not masked:
        #Time to reject cosmic rays
        clipped=crreject(imcube, **crparams)
        im=imcube.mean(axis=2)
        var=np.ma.array(varcube,mask=clipped.mask).sum(axis=2)
    else:
        #Create masked arrays
        #Average the data count rates wieghted by the variance
        imcube_masked=np.ma.array(imcube,mask=mask)
        if exptime != 0.0:
            im=np.ma.average(imcube_masked, axis=2,
                             weights=varcube/np.array(etimes))
        else:
            im=np.ma.average(imcube_masked, axis=2)
        im=im.filled(float(0))
        im*=exptime
        #Sum the variance frames
        var=np.ma.array(varcube,mask=mask).sum(axis=2).filled(float('nan'))

    #Write out the merged file
    header['FILENAME']=os.path.basename(outfile)
    header['EXPTIME']=exptime
    header['COMMENT']=','.join(members)
    hdu = fits.PrimaryHDU(im.astype(np.float32), header=header)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(var.astype(np.float32),name='variance'))
    hdul.append(fits.ImageHDU(np.sum(mask,axis=2).astype(np.uint8),name='crmask'))
    hdul.append(fits.ImageHDU(np.zeros_like(im,dtype=np.uint8), name='bpmask'))
    hdul.writeto(outfile+'.fits')


if __name__ =='__main__':
    
    args=parse_cl()
    
    files=[os.path.join(dirpath, f)
            for dirpath, dirnames, files in os.walk(args.dir)
            for f in files if '.fits' in f and '-' not in f and ',' not in f] #allow fits.gz
    try:
        seqno_stacks=get_seqlists(args.listfile)
        
        to_stack_lists=[[f for f in files
                         if int(m2fs.obs.info(f)['seqno']) in seqnos
                         and side==m2fs.obs.info(f)['side']]
                        for seqnos in seqno_stacks for side in 'rb']
                        
        to_stack_lists=[l for l in to_stack_lists if len(l) >1] #only stack 2 or more
    
    except IOError:
        IOError('No listfile')

    for i,filestack in enumerate(to_stack_lists):
        print "Stacking {} of {} ({} files)".format(i,
                len(to_stack_lists), len(filestack))
        seqnos=map(lambda f: int(m2fs.obs.info(f)['seqno']),filestack)
        filestack=sorted(filestack,
                         key=lambda x: seqnos.__getitem__(filestack.index(x)) )
        color=m2fs.obs.info(filestack[0])['side']
        try:
            if (os.path.exists(args.outdir+color+rangify(seqnos)+'.fits') or
                os.path.exists(args.outdir+color+rangify(seqnos)+'.fits.gz')):
                continue
            
            print '   Stacking ',color+rangify(seqnos)
            stackimage(filestack,args.outdir+color+rangify(seqnos),
                       do_cosmic=args.do_cosmic)
        except IOError as e:
            print "Couldn't stack {}".format(str(e))
