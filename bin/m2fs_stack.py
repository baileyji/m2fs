#!/usr/bin/env python
import numpy as np
import sys
from astropy.io import fits
import re
import pylab as pl
import astropy.stats
import scipy
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
                     action='store', required=False, type=str,
                     help='testfile with files to proces')
    parser.add_argument('--crreject', dest='do_cosmic', default=False,
                     action='store_true', required=False,
                     help='Do cosmic ray rejection')
    parser.add_argument('-o', dest='outdir', default='./',
                     action='store', required=False, type=str,
                     help='Out directory')
    parser.add_argument('-z', dest='gzip', default=False,
                 action='store_true', required=False,
                 help='gzip fits files')
    parser.add_argument('-s', dest='single', default='',
                 action='store', required=False, type=str,
                 help='Single file to update')
    parser.add_argument('--overwrite', dest='clobber', default=False,
                 action='store_true', required=False,
                 help='Overwite existing files')
    parser.add_argument('-t', dest='dry_run', default=False,
                 action='store_true', required=False,
                 help='Test, dont do anything')

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

def compute_patch_val(variances):
    patch_val=max(1.05*variances.max(), 1e8)
    if patch_val>=1e38:
        print('WARNING: Real variances {}'.format(patch_val)+
              'exceed maximum bad variance patch value (1e38) for float32!')
    return min(patch_val,1e38)

def mask_problem_spots(var, header, mask_val=1e35):
    if 'r' in header['FILENAME']:
        if header['BINNING']=='1x1':
            var[:2056,821]=mask_val
        if header['BINNING']=='2x2':
            var[:2056,410]=mask_val
    if 'b' in header['FILENAME'] and header['BINNING']=='1x1':
        var[328,:2048]=mask_val
        var[1600:2056,2518]=mask_val

def stackimage(files, outfile,  gzip=False, do_cosmic=False, clobber=False,
               **crparams):
    """
    List of files to stack, output file sans extension.
    
    
    output hdu
    cleaned image
    variance image
    summed cr mask
    bad pixel mask
    
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
                        msk+=ndimage.morphology.binary_dilation(
                            foo.astype(np.bool), structure=np.ones((3,3)),
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


    if not masked and do_cosmic:
        #Time to reject cosmic rays
        clipped=crreject(imcube, **crparams)
        mask=clipped.mask
    elif not masked:
        mask=np.zeros_like(imcube, dtype=bool)

    #this code issues warnings like crazy, but the math is ok because of a
    #final pass collecitng all the garbage

    #Create masked arrays
    imcube_masked=np.ma.array(imcube, mask=mask)
    varcube_masked=np.ma.array(varcube, mask=mask)
    
    duration_corr=max(etimes)/np.array(etimes)
    
    throughput_corr=np.array(midpoint_weights)/np.array(etimes)
    throughput_corr=throughput_corr.max()/throughput_corr

    ####### Try 3
    
    #the correction factor as defined has a pathological edge case:
    #if high throughput fibers or brighter targets are disproportionalty
    #affected by clouds/field rotation, guiding errors, or the like then it
    #would give unfair weight to that frame.

    corr_fac=throughput_corr*duration_corr
    patch_cube=(imcube_masked*corr_fac).mean(axis=2)
    patch_cube=patch_cube[:,:,np.newaxis]/corr_fac
    imcube[mask]=patch_cube[mask]
    im=imcube.sum(axis=2)

    #fill value doesn't matter as mask check below will flag and cause patch
    var=varcube_masked.sum(axis=2).filled(0)
    var*=(im/imcube_masked.sum(axis=2))**2
    
    mask_problem_spots(var, header, mask_val=np.nan)

    im=im.astype(np.float32)
    var=var.astype(np.float32)

    bad=(mask.all(axis=2) | (~np.isfinite(im)) | (~np.isfinite(var)))
    
    patch_val=compute_patch_val(var[~bad])
    im[bad]=0.0
    var[bad]=patch_val

    if (~np.isfinite(var)).any():
        import ipdb;ipdb.set_trace()
    if (~np.isfinite(im)).any():
        import ipdb;ipdb.set_trace()

    #Update the header
    header['FILENAME']=os.path.basename(outfile)
    header['EXPTIME']=exptime
    header['COMMENT']=','.join(members)
    header['UT-MID']=str(midpoint)
    header['UT-TIME']=str(first_uttime)
    header['UT-END']=str(last_utend)
    header['AIRMASS']=mean_airmass
    header['PATCHVAL']=patch_val
    
    print('Setting read noise and to approximate value. '
          'Use variance frame for exact value')
    header['EGAIN']=1.0
    header['ENOISE']=np.sqrt(nfile)*2.5 # average for the 4 amps in slow is 2.5 e
    
    #Create a primary extension with header only
    #   (this seems to be a fits convention)
    hdul = fits.HDUList(fits.PrimaryHDU(header=header))
    
    #Append image extensions
    hdul.append(fits.ImageHDU(im, name='science', header=header))
    hdul.append(fits.ImageHDU(var, name='variance', header=header))
    hdul.append(fits.ImageHDU(mask.sum(axis=2).astype(np.uint8),
                              name='crmask'))
    hdul.append(fits.ImageHDU(np.zeros_like(im,dtype=np.uint8),
                              name='bpmask'))
    hdul.writeto(outfile+'.fits'+('.gz' if gzip else ''), clobber=clobber)



def update_single(file, outfile,  gzip=False, do_cosmic=False, clobber=False,
                  **extra):
    hdul=fits.open(file)
    mv=compute_patch_val(hdul['VARIANCE'].data)
    mask_problem_spots(hdul['VARIANCE'].data,
                       hdul['SCIENCE'].header, mask_val=mv)
    hdul.writeto(outfile+'.fits'+('.gz' if gzip else ''), clobber=clobber)

if __name__ =='__main__':
    
    args=parse_cl()
    
    if args.single:
        outf=os.path.join(args.outdir,
                          os.path.basename(args.single).replace('.fits',''))
        update_single(args.single, outf,
                           gzip=args.gzip, clobber=args.clobber)
        sys.exit()
    
    
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
            if ((os.path.exists(args.outdir+color+
                                rangify(seqnos,delim='_')+'.fits') or
                 os.path.exists(args.outdir+color+
                                rangify(seqnos,delim='_')+'.fits.gz')) and not
                args.clobber):
                continue
            
            print '   Stacking ',color+rangify(seqnos)
            if not args.dry_run:
                stackimage(filestack, args.outdir+color+
                           rangify(seqnos,delim='_'),
                           gzip=args.gzip, do_cosmic=args.do_cosmic,
                           clobber=args.clobber)

        except IOError as e:
            print "Couldn't stack {}".format(str(e))
