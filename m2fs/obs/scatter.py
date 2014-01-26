#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from scipy.ndimage import median_filter
from matplotlib.pyplot import *


if __name__ =='__main__':
    args=parse_cl()
    trfile='stacked/r856-860.fits.gz'
    hdul=fits.open(trfile)

    im=hdul[0].data
    scatter_regions=[(2,75), (480,575), (1010,1100), (1525,1600), (2035,2120),
                     (2540,2620), (3030,3110), (3525, 3600), (3995 ,4109)]

    scat_disp=m2fs.obs.scatterxfirst(im,scatter_regions,prof_order=8)

    #scat_prof=m2fs.obs.scatter(im,scatter_regions,prof_order=8)
    #
    plt.figure(1)
    plt.imshow(scat_disp,origin='lower')
    plt.title('disp: max={} min={} mean={}'.format(
                scat_disp.max(),
                scat_disp.min(),
                scat_disp.mean()))

    plt.figure(2)
    plt.imshow(scat_disp-scat_prof,origin='lower')
    plt.title('disp-prof: max={} min={} mean={}'.format(
              (scat_disp-scat_prof).max(),
              (scat_disp-scat_prof).min(),
              (scat_disp-scat_prof).mean()))

    plt.figure(3)
    plt.imshow(im-scat_disp,origin='lower')
    plt.title('im-disp: max={} min={} mean={}'.format(
               (im-scat_disp).max(),
               (im-scat_disp).min(),
               (im-scat_disp).mean()))

    plt.figure(4)
    x=2000
    plt.plot((im-scat_disp)[:,x:x+100].sum(axis=1)/100)
    plt.title('im-disp')
    plt.figure(5)
    plt.plot((im-scat_prof)[:,x:x+100].sum(axis=1)/100)
    plt.title('im-prof')

    plt.show()


def scatter(im_in, scatter_regions, plotalot=False,
            prof_order=6, disp_order=11, disp_bin=64, dispFitst=False):
    """Create a polynomial surface of the scattered light in the frame"""
    msk=np.ones_like(im_in,dtype=np.bool)
    for r0,r1 in scatter_regions:
        msk[r0:r1+1,:]=False
    #Useing a masked array doesn't seem to work with median_filter
    #im=np.ma.array(im_in,mask=msk,fill_value=np.nan)
    im=im_in.copy()
    im[msk]=np.nan

    if plotalot:
        figure(1)
        foo=im_in.copy()
        foo[~msk]=np.nan
        imshow(foo, origin='lower')
        show()
    nrow,ncol=im.shape

    bin_col=np.arange(disp_bin/2, ncol, disp_bin)
    rows=np.arange(nrow)
    cols=np.arange(ncol)

    binned=median_filter(im, size=(1,disp_bin) )
    binned[msk]=np.nan
    binned=binned[:,bin_col]

    bin_poly_coeff=np.zeros((ncol/disp_bin, prof_order+1))
    for i,column in enumerate(binned.T):
        bin_poly_coeff[i,:]=np.polyfit(rows[~np.isnan(column)],
                                       column[~np.isnan(column)],
                                       prof_order)
        if plotalot:
            figure(2)
            line,=plt.plot(rows, column)
            plt.plot(rows, np.poly1d(bin_poly_coeff[i,:])(rows),
                     color=line.get_color())

    if plotalot:
        show()


    #Fit polynomials to each term
    coeff_poly=[]
    for i, coeff in enumerate(bin_poly_coeff.T):
        coeff_poly.append(np.poly1d(np.polyfit(bin_col, coeff, disp_order)))

    if plotalot:
        figure(3)
        fig, subplots = plt.subplots(3,3)
        subplots=subplots.flatten()
        for i,coeff in enumerate(bin_poly_coeff.T):
            subplots[i].plot(bin_col,coeff)
            #subplots[i].title('Coeff {}'.format(i))
            subplots[i].plot(cols,coeff_poly[i](cols))
        show()

    #Use polynomials to generate the image
    scatterim=np.zeros_like(im)
    for j in range(ncol):
        scatterim[:,j]=np.poly1d([poly(j) for poly in coeff_poly])(rows)

    return scatterim

def scatterxfirst(im_in, scatter_regions, plotalot=False, disp_bin=32,
            prof_order=8, disp_order=11, prof_bin=1):
    """Create a polynomial surface of the scattered light in the frame"""
    
    im=im_in.copy()

    nrow,ncol=im.shape
    rows=np.arange(nrow)
    cols=np.arange(ncol)

    #Create mask for scatterd light regions
    #True=bad
    msk=np.ones_like(im_in,dtype=np.bool)
    row_msk=np.ones(nrow,dtype=np.bool)
    for r0,r1 in scatter_regions:
        row_msk[r0:r1+1]=False
        msk[r0:r1+1,:]=False
    
    #Using a masked array doesn't seem to work with median_filter
    #im=np.ma.array(im_in,mask=msk,fill_value=np.nan)
    im[msk]=np.nan

    #Show the regions not being used
    if plotalot:
        figure(1)
        foo=im_in.copy()
        foo[~msk]=np.nan
        imshow(foo, origin='lower')
        show()

    #Configure for dispersion or profile fit first
    xx=np.arange(disp_bin/2, ncol, disp_bin)
    poly_xx=rows

    poly_ndxs=rows[~row_msk]

    ord=disp_order
    poly_ord=prof_order

    #Bin the image
    binned=median_filter(im, size=(prof_bin, disp_bin) )
    binned[msk]=np.nan

    binned=binned[:,xx]


    #Fit polynomials to the the scattered light along xx

    poly_coeff=np.zeros( (len(poly_ndxs), ord+1) )

    for i,ndx in enumerate(poly_ndxs):
        yy=binned[ndx]
        poly_coeff[i,:]=np.polyfit(xx[~np.isnan(yy)], yy[~np.isnan(yy)], ord)
        if plotalot:
           figure(2)
           line,=plt.plot(xx, yy)
           plt.plot(xx, np.poly1d(poly_coeff[i,:])(xx), color=line.get_color())

    if plotalot:
        show()

    #Fit polynomials to each term of the polynomials found
    polys=[]
    for coeff in poly_coeff.T:
        polys.append(np.poly1d(np.polyfit(poly_ndxs, coeff, poly_ord)))


    if plotalot:
        fig, subplots = plt.subplots(3,4)
        subplots=subplots.flatten()
        for i,coeff in enumerate(poly_coeff.T):
            subplots[i].plot(poly_ndxs, coeff,'.')
            subplots[i].plot(poly_xx, polys[i](poly_xx))
        show()
    
    #Use polynomials to generate the image
    scatterim=np.zeros_like(im)
    for x in poly_xx:
        scatterim[x,:]=np.poly1d([poly(x) for poly in polys])(cols)
    
    return scatterim




