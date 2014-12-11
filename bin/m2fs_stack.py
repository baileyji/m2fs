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
from astropy.time import Time, TimeDelta

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
    parser.add_argument('-z', dest='gzip', default=False,
                 action='store', required=False, type=bool,
                 help='gzip fits files')
    
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
            if l[0] in '1234567890':
                range=l.split()[0]
                ret.add(tuple(map(str,derangify(range))))
            elif len(l)>1 and l[0] in 'RBrb' and l[1] in '1234567890':
                range=l[1:].split()[0]
                ret.add(tuple(map(lambda x: l[0].lower()+str(x),
                                  derangify(range))))
    return list(ret)




def stackimage(files, outfile,  gzip=False, do_cosmic=False, **crparams):
    """
    List of files to stack, output file sans extension.
    
    
    output hdu
    cleaned image
    variance image
    summed cr mask
    bad pixel mask
    
    
    So this is pretty far afield from what it should be
    
    ci= estimate of total frame throughput relative to frame with most
    n = number of frames
    k=number of corrupted frames at given pixel
    
    s = (n+k)/n * sum(pixel_i/ci)
    v = ((n+k)/n)^2 sum((pixel_i+rn^2)/ci^2)
    
    
    """

    
    nfile=len(files)
    #load first bias to get info
    with fits.open(files[0]) as im:
        header=im[0].header
        #15 frames with float64 would be 3.7GB ram
        imcube=np.zeros(im[1].data.shape+(nfile,),dtype=np.float32)
        if len(im) >3:
            masked=True
            mask=np.zeros_like(imcube, dtype=np.bool)
        else:
            masked=False

    varcube=np.zeros_like(imcube)

    #load data and variance frames
    members=[]
    etimes=[]
    midpoint_weights=[]
    midpoints=[]
    airmasses=[]
    first_uttime=None
    last_utend=None
    try:
        for i,f in enumerate(files):
            with fits.open(f) as im:
                airmasses.append(im[0].header['AIRMASS'])
                etimes.append(float(im[0].header['EXPTIME']))
                members.append(im[0].header['FILENAME'])

                imcube[:,:,i]=im[1].data
                varcube[:,:,i]=im[2].data

                if masked:
                    msk=im[3].data
                    if do_cosmic:
                        foo=crreject(im[2].data)
                        msk+=foo
                    msk=ndimage.morphology.binary_dilation(
                        msk.astype(np.bool), structure=np.ones((3,3)),
                        iterations=1, mask=None, output=None,
                        border_value=0, origin=0, brute_force=False)
                
                    mask[:,:,i]=msk
                    
                #compute midpoint
                start=Time(im[0].header['UT-DATE']+' '+im[0].header['UT-TIME'],
                           format='iso', scale='utc')
                end=Time(im[0].header['UT-DATE']+' '+im[0].header['UT-END'],
                         format='iso', scale='utc')
                midpoints.append(start +.5*(end-start))
                
                #keep track of stacked beginning and end
                if type(first_uttime) == type(None) or first_uttime>start:
                    first_uttime=start
#                if first_uttime>start: first_uttime=start
                if type(last_utend) == type(None) or last_utend<end:
                    last_utend=end
#                if last_utend<end: last_utend=end

                #import ipdb;ipdb.set_trace() #midpoints untested


    except ValueError as e:
        print "ValueError while merging:", files,f
        return


    if masked:
        use=~mask.sum(2).astype(bool)
        midpoint_weights=[im.T[use].sum() for im in imcube.T]
    else:
        midpoint_weights=[im.sum() for im in imcube.T]


    mean_airmass=np.mean(airmasses)
    min_midpoint=min(midpoints)
    midpoint=min_midpoint + TimeDelta(np.average([(m - min_midpoint).sec
                                                  for m in midpoints],
                                                 weights=midpoint_weights),
                                      format='sec')

    exptime=sum(etimes)


    if not masked:
        #Time to reject cosmic rays
        clipped=crreject(imcube, **crparams)
        im=imcube.mean(axis=2)
        var=np.ma.array(varcube,mask=clipped.mask).sum(axis=2)
    else:
        
        #Create masked arrays
        imcube_masked=np.ma.array(imcube, mask=mask)
        varcube_masked=np.ma.array(varcube, mask=mask)
        
        duration_corr=max(etimes)/np.array(etimes)
        
        throughput_corr=np.array(midpoint_weights)/np.array(etimes)
        throughput_corr=throughput_corr.max()/throughput_corr
        

        #Have correction factors
        #s = n/(n-k) * sum(pixel_i * duration_corr*throughput_corr)
        
#        mask_corr=nfile/(nfile-mask.sum(axis=2).astype(float))
#
#        im=mask_corr*(imcube_masked*throughput_corr*duration_corr).sum(axis=2)
#        im=im.filled(0.0)
#
#        var=mask_corr**2 * (varcube_masked*
#                            (throughput_corr*duration_corr)**2).sum(axis=2)
#        var=var.filled(1e99)

        ####### Try 3
        
#        the correction factor as defined has a pathological edge case:
#        if high throughput fibers or brighter targets are disproportionalty
#        affected by clouds/field rotation, guiding errors, or the like then it
#        would give unfair weight to that frame.

        corr_fac=throughput_corr*duration_corr
        patch_cube=(imcube_masked*corr_fac).mean(axis=2)
        patch_cube=patch_cube[:,:,np.newaxis]/corr_fac
        imcube[mask]=patch_cube[mask]
        im=imcube.sum(axis=2)
        
        var=varcube_masked.sum(axis=2).filled(1e99)
        var*=(im/imcube_masked.sum(axis=2))**2
        
        bad=((~mask).sum(axis=2)==0) | (~np.isfinite(im)) | (~np.isfinite(var))
        
        im[bad]=0.0
        var[bad]=1e99


        #hack for a couple of bad rows
        if 'r' in header['FILENAME'] and header['BINNING']=='1x1':
            var[:2056,821]=1e99

#        import ipdb;ipdb.set_trace()
        ######




    #Update the header
    header['FILENAME']=os.path.basename(outfile)
    header['EXPTIME']=exptime
    header['COMMENT']=','.join(members)
    header['UT-MID']=str(midpoint)
    header['UT-TIME']=str(first_uttime)
    header['UT-END']=str(last_utend)
    header['AIRMASS']=mean_airmass
    
    print('Setting read noise and to approximate value. '
          'Use variance frame for exact value')
    header['EGAIN']=1.0
    header['ENOISE']=np.sqrt(nfile)*2.5 # average for the 4 amps in slow is 2.5 e
    
    #Create a primary extension with header only
    #   (this seems to be a fits convention)
    hdul = fits.HDUList(fits.PrimaryHDU(header=header))
    
    #Append image extensions
    hdul.append(fits.ImageHDU(im.astype(np.float32),
                              name='science', header=header))
    hdul.append(fits.ImageHDU(var.astype(np.float32),
                              name='variance', header=header))
    hdul.append(fits.ImageHDU(mask.sum(axis=2).astype(np.uint8),
                              name='crmask'))
    hdul.append(fits.ImageHDU(np.zeros_like(im,dtype=np.uint8),
                              name='bpmask'))
    hdul.writeto(outfile+'.fits'+('.gz' if gzip else ''))


if __name__ =='__main__':
    
    args=parse_cl()
    
    files=[os.path.join(dirpath, f)
            for dirpath, dirnames, files in os.walk(args.dir)
            for f in files
            if '.fits' in f and '-' not in f and
                ',' not in f and 'c' not in f]
    try:
        seqno_stacks=get_seqlists(args.listfile)

        to_stack_lists=[[f for f in files
                         if m2fs.obs.info(f, no_load=True).seqno_match(seqnos)
                         and side==m2fs.obs.info(f,no_load=True).side]
                         for seqnos in seqno_stacks for side in 'rb']
                        
        to_stack_lists=[l for l in to_stack_lists if len(l) >1] #only stack 2 or more
    
    except IOError:
        IOError('No listfile')

    for i,filestack in enumerate(to_stack_lists):
        print "Stacking {} of {} ({} files)".format(i,
                len(to_stack_lists), len(filestack))
        seqnos=map(lambda f: int(m2fs.obs.info(f, no_load=True).seqno),
                   filestack)
        filestack=sorted(filestack,
                         key=lambda x: seqnos.__getitem__(filestack.index(x)) )
        color=m2fs.obs.info(filestack[0], no_load=True).side
        try:
            if (os.path.exists(args.outdir+color+rangify(seqnos)+'.fits') or
                os.path.exists(args.outdir+color+rangify(seqnos)+'.fits.gz')):
                continue
            
            print '   Stacking ',color+rangify(seqnos)
            stackimage(filestack,args.outdir+color+rangify(seqnos),
                       gzip=args.gzip, do_cosmic=args.do_cosmic)
        except IOError as e:
            print "Couldn't stack {}".format(str(e))
