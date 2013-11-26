#!/usr/bin/env python2.7
import numpy as np
from astropy.io import fits
import re
from pylab import *
import astropy.stats


def mergeimage(frameno,side=None):
    if side == None:
        try:
            mergeimage(frameno,side='r')
        except IOError:
            print 'Need all quadrants for r{}'.format(frameno)
        try:
            mergeimage(frameno,side='b')
        except IOError:
            print 'Need all quadrants for b{}'.format(frameno)
        return

    basename=side+'{frameno:04}'
    f=basename+'c{quad}.fits'

    try:
        quadrant1=fits.open(f.format(frameno=frameno, quad=1))
        quadrant2=fits.open(f.format(frameno=frameno, quad=2))
        quadrant3=fits.open(f.format(frameno=frameno, quad=3))
        quadrant4=fits.open(f.format(frameno=frameno, quad=4))
    except IOError, e:
        raise e

    #Grab the bias and crop region from the first quadrant
    biassec=quadrant1[0].header['BIASSEC']
    bias=[int(s) for s in re.findall(r'\d+',biassec)]

    trimsec=quadrant1[0].header['TRIMSEC']
    crop=[int(s) for s in re.findall(r'\d+',trimsec)]

    #Create an output array
    out=np.zeros((crop[3]*2,crop[1]*2),dtype=np.float)

    #Define where the quadrants show go
    quadLoc=[(0, 2056, 0, 2048),
            (0,2056, 2048, 4096),
             ( 2056, 4112, 2048, 4096),
             (2056, 4112,0, 2048)]

    quadLoc=[(0, crop[3], 0, crop[1]),
             (0,crop[3], crop[1], 2*crop[1]),
             ( crop[3], 2*crop[3], crop[1], 2*crop[1]),
             (crop[3],2*crop[3],0, crop[1])]

    #Create a list of the quadrant's data to loop through
    quadrantData=[quadrant1[0].data, quadrant2[0].data,
                 quadrant3[0].data, quadrant4[0].data]

    for i, qdata in enumerate(quadrantData):
        
        #Compute & subtract mean bias row
        #biasrow=np.mean(qdata[bias[2]:,:],axis=0)
        #qdata-=biasrow
        
        #Compute mean for overscan region rowwise
        biaslevels=np.median(qdata[crop[2]-1:crop[3],bias[0]-1:bias[1]],axis=1)

        #Crop image & subtract bias levels row wise
        qdata=(qdata[crop[0]-1:crop[3],crop[0]-1:crop[1]] -
                          biaslevels[:,np.newaxis])
                          
        #Convert to e-
        qdata*=header['GAIN']
        vimage=qdata**2+header['ENOISE']**2

        if i ==0:
            out[quadLoc[i][0]:quadLoc[i][1],quadLoc[i][2]:quadLoc[i][3]]=qdata
        if i==1:
            out[quadLoc[i][0]:quadLoc[i][1],
                quadLoc[i][2]:quadLoc[i][3]]=np.fliplr(qdata)
        if i==2:
            out[quadLoc[i][0]:quadLoc[i][1],
                quadLoc[i][2]:quadLoc[i][3]]=np.rot90(qdata,2)
        if i==3:
            out[quadLoc[i][0]:quadLoc[i][1],
                quadLoc[i][2]:quadLoc[i][3]]=np.rot90(np.fliplr(qdata),2)



    #Flip so it is in agreement with Mario's process
    out=np.flipud(out)
    
    #Write out the merged file
    hdu = fits.PrimaryHDU(out) TODO
    new_hdul = pyfits.HDUList()
    new_hdul.append(pyfits.ImageHDU())
    new_hdul.append(pyfits.ImageHDU())

    hdu.header=quadrant1[0].header
    hdu.header.pop('TRIMSEC')
    hdu.header.pop('BIASSEC')
    hdu.header.pop('DATASEC')
    hdu.header['FILENAME']=basename.format(frameno=frameno)
    hdu.writeto(basename.format(frameno=frameno)+'.fits')


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

if __name__ =='__main__':
    import glob
    files = glob.glob('*c?.fits')
    seqnos=set([ int(x[1:-7]) for x in files])
    for i in seqnos:
        mergeimage(i)
