#!/usr/bin/env python2.7
import numpy as np
from astropy.io import fits
import astropy.stats
import m2fs.ccd
import re
import os
from jbastro.lacosmics.cosmics import cosmicsimage
import m2fs.obs

from multiprocessing import Pool

def proc_quad(qdata, header, cosmic_settings):
    """cosmic_settings = cosmic_kwargs + 'iter' or just 'iter':0"""
    
    cosmic_iter=cosmic_settings.pop('iter')
    
    
    biassec=header['BIASSEC']
    bias=[int(s) for s in re.findall(r'\d+',biassec)]
    
    trimsec=header['TRIMSEC']
    crop=[int(s) for s in re.findall(r'\d+',trimsec)]
    
    #Compute & subtract mean bias row
    #biasrow=np.mean(qdata[bias[2]:,:],axis=0)
    #qdata-=biasrow
    
    #Compute mean for overscan region rowwise
    biaslevels=np.median(qdata[crop[2]-1:crop[3],bias[0]-1:bias[1]],axis=1)
    
    #Crop image & subtract bias levels row wise
    qdata=(qdata[crop[0]-1:crop[3],crop[0]-1:crop[1]] -
           biaslevels[:,np.newaxis])
           
    #Cosmic ray rejection
    if cosmic_iter>0:
        c=cosmicsimage(qdata, gain=header['EGAIN'],
                      readnoise=header['ENOISE'],
                      satlevel=.95*m2fs.ccd.satlevel,
                      **cosmic_settings)
        c.run(maxiter = cosmic_iter)
        qmask=c.mask
    else:
        qmask=np.zeros_like(qdata, dtype=np.bool)

    #Convert to e-
    qdata*=header['EGAIN']

    #Variance image
    qvari=qdata + header['ENOISE']**2
    
    return (qdata, qvari, qmask)


def mergequad(frameno, side=None, do_cosmic=False, file=False, odir=''):
    """Give a seqno or a path to a quad if file set
    do_cosmic=bool or dict like
    sigma for init clip, fraction for neighbors, how much above background*
    {'sigclip': 7.0, 'sigfrac': 0.428, 'objlim': 1.4, 'iter':10}
    """
    if side == None and not file:
        try:
            mergequad(frameno,side='r', do_cosmic=do_cosmic)
        except IOError:
            print 'Need all quadrants for r{}'.format(frameno)
        try:
            mergequad(frameno,side='b',do_cosmic=do_cosmic)
        except IOError:
            print 'Need all quadrants for b{}'.format(frameno)
        return
    
    if file:
        fname=frameno
        frameno=int(m2fs.obs.info(fname).seqno)
        basename=os.path.basename(fname)
        side=basename[0]
        dir=os.path.dirname(fname)+os.path.sep
        _,_,extension=basename.partition('.')
    else:
        dir=''
        extension='fits'
    
    basename=side+'{frameno:04}'
    f=dir+basename+'c{quad}.'+extension
    
    if (os.path.exists(odir+basename.format(frameno=frameno)+'.fits') or
        os.path.exists(odir+basename.format(frameno=frameno)+'.fits.gz')):
        print ("Skipping "+odir+basename.format(frameno=frameno)+".fits. "
               "Already done.")
        return
    else:
        print "Merging "+odir+basename.format(frameno=frameno)+'.fits'

    #Load the data
    try:
        quadrant1=fits.open(f.format(frameno=frameno, quad=1))
        quadrant2=fits.open(f.format(frameno=frameno, quad=2))
        quadrant3=fits.open(f.format(frameno=frameno, quad=3))
        quadrant4=fits.open(f.format(frameno=frameno, quad=4))
    except IOError, e:
        raise e

    #Cosmic ray removal configuration
    if type(do_cosmic)==bool:
        if do_cosmic:
            cosmic_settings={'sigclip': 7.0, 'sigfrac': 0.428, 'objlim': 1.4,
                    'iter':10}
            do_cosmic=cosmic_settings.copy()
        else:
            cosmic_settings={'iter':0}
    else:
        cosmic_settings=do_cosmic.copy()
    

    #Grab the bias and crop region from the first quadrant
    biassec=quadrant1[0].header['BIASSEC']
    bias=[int(s) for s in re.findall(r'\d+',biassec)]
    
    trimsec=quadrant1[0].header['TRIMSEC']
    crop=[int(s) for s in re.findall(r'\d+',trimsec)]
    
    #Create output arrays
    out=np.zeros((crop[3]*2,crop[1]*2),dtype=np.float32)
    vout=np.zeros((crop[3]*2,crop[1]*2),dtype=np.float32)
    mask=np.zeros((crop[3]*2,crop[1]*2),dtype=np.uint8)
    headerout=quadrant1[0].header.copy()
    
    #Define where the quadrants show go
    quadLoc=[(0, crop[3], 0, crop[1]),
             (0,crop[3], crop[1], 2*crop[1]),
             ( crop[3], 2*crop[3], crop[1], 2*crop[1]),
             (crop[3],2*crop[3],0, crop[1])]
             
    #Create a list of the quadrant's data to loop through
    quadrantData=[quadrant1[0], quadrant2[0], quadrant3[0], quadrant4[0]]


    #Process them in multiple processes
    args=[(q.data.copy(),q.header.copy(), cosmic_settings.copy())
          for q in quadrantData]
    pool = Pool(processes=4)
    res=[pool.apply_async(proc_quad, arg) for arg in args]
    pool.close()
    pool.join()

    for i,r in enumerate(res):
        qdata, qvar, qmask = r.get()
        #Position the quadrants
        if i ==0:
            out[quadLoc[i][0]:quadLoc[i][1],quadLoc[i][2]:quadLoc[i][3]]=qdata
            vout[quadLoc[i][0]:quadLoc[i][1],quadLoc[i][2]:quadLoc[i][3]]=qvar
            mask[quadLoc[i][0]:quadLoc[i][1],quadLoc[i][2]:quadLoc[i][3]]=qmask
        if i==1:
            out[quadLoc[i][0]:quadLoc[i][1],
                quadLoc[i][2]:quadLoc[i][3]]=np.fliplr(qdata)
            vout[quadLoc[i][0]:quadLoc[i][1],
                 quadLoc[i][2]:quadLoc[i][3]]=np.fliplr(qvar)
            mask[quadLoc[i][0]:quadLoc[i][1],
                 quadLoc[i][2]:quadLoc[i][3]]=np.fliplr(qmask)
        if i==2:
            out[quadLoc[i][0]:quadLoc[i][1],
                quadLoc[i][2]:quadLoc[i][3]]=np.rot90(qdata,2)
            vout[quadLoc[i][0]:quadLoc[i][1],
                 quadLoc[i][2]:quadLoc[i][3]]=np.rot90(qvar,2)
            mask[quadLoc[i][0]:quadLoc[i][1],
                 quadLoc[i][2]:quadLoc[i][3]]=np.rot90(qmask,2)
        if i==3:
            out[quadLoc[i][0]:quadLoc[i][1],
                quadLoc[i][2]:quadLoc[i][3]]=np.rot90(np.fliplr(qdata),2)
            vout[quadLoc[i][0]:quadLoc[i][1],
                 quadLoc[i][2]:quadLoc[i][3]]=np.rot90(np.fliplr(qvar),2)
            mask[quadLoc[i][0]:quadLoc[i][1],
                 quadLoc[i][2]:quadLoc[i][3]]=np.rot90(np.fliplr(qmask),2)

    #Flip so it is in agreement with Mario's process
    out=np.flipud(out)
    vout=np.flipud(vout)
    mask=np.flipud(mask)

    #Write out the merged file
    headerout.pop('TRIMSEC')
    headerout.pop('BIASSEC')
    headerout.pop('DATASEC')
    headerout['FILENAME']=basename.format(frameno=frameno)
    headerout['BUNIT']='E-/PIXEL'

    hdu = fits.PrimaryHDU(out, header=headerout)
    hdul = fits.HDUList([hdu])
    hdul.append(fits.ImageHDU(vout,name='variance'))
    if do_cosmic:
        hdul.append(fits.ImageHDU(mask,name='mask'))
    if not os.path.exists(odir+basename.format(frameno=frameno)+'.fits'):
        hdul.writeto(odir+basename.format(frameno=frameno)+'.fits')

def makesuperbias(filenos, side, name):
    
    f=side+'{frameno:04}.fits'
    
    #load first bias to get info
    with fits.open(f.format(frameno=filenos[0])) as bias:
        header=bias[0].header
        superbias_data=np.zeros_like(bias[0].data,dtype=np.float32)
    
    #Number of bias frames
    nbias=len(filenos)
    
    #Sum the bias counts
    for num in filenos:
        with fits.open(f.format(frameno=num)) as bias:
            superbias_data+=bias[0].data
    
    #Merge bias
    superbias_data/=nbias
    
    #Write out the merged file
    hdu = fits.PrimaryHDU(superbias_data)
    hdu.header=header
    hdu.header['FILENAME']=name
    hdu.header['COMMENT']=','.join(map(str,filenos))
    hdu.writeto(name+'.fits')



def makesuperflat(filenos, side, name, superbias=None):
    
    f=side+'{frameno:04}.fits'
    
    #load first flat to get info
    with fits.open(f.format(frameno=filenos[0])) as flat:
        header=flat[0].header
        superflat_data=np.zeros_like(flat[0].data,dtype=np.float32)
        nrow,ncol=superflat_data.shape
    
    #Number of flat frames
    nflat=len(filenos)
    
    
    #Grab the bias
    try:
        biasdata=fits.open(superbias)[0].data
    except IOError:
        biasdata=None
    
    #Make the datacube
    cube=np.zeros((nrow, ncol, nflat), dtype=np.float32)
    
    
    for i,num in enumerate(filenos):
        with fits.open(f.format(frameno=num)) as flat:
            #Remove bias
            if biasdata != None:
                cube[:,:,i]=flat[0].data-biasdata
            else:
                cube[:,:,i]=flat[0].data
    
    #Merge flat
    masked=astropy.stats.sigma_clip(cube,sig=3,axis=2)
    superflat_data=np.ma.MaskedArray.mean(masked,axis=2).data
    
    #Write out the superflat
    hdu = fits.PrimaryHDU(superflat_data)
    hdu.header=header
    hdu.header['FILENAME']=name
    hdu.header['COMMENT']=','.join(map(str,filenos))+',Bias:'+str(superbias)
    hdu.writeto(name+'.fits')

