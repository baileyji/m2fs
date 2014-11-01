#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from scipy.ndimage import median_filter
from matplotlib.pyplot import *


from scipy.ndimage.filters import gaussian_filter as gsmooth2d
import numpy as np
import matplotlib.pyplot as plt
from jbastro.astroLib import gauss2D, gaussfit2D


#if __name__ =='__main__':
#    args=parse_cl()
#    trfile='stacked/r856-860.fits.gz'
#    hdul=fits.open(trfile)
#
#    im=hdul[0].data
#    scatter_regions=[(2,75), (480,575), (1010,1100), (1525,1600), (2035,2120),
#                     (2540,2620), (3030,3110), (3525, 3600), (3995 ,4109)]
#
#    scat_disp=m2fs.obs.scatterxfirst(im,scatter_regions,prof_order=8)
#
#    #scat_prof=m2fs.obs.scatter(im,scatter_regions,prof_order=8)
#    #
#    plt.figure(1)
#    plt.imshow(scat_disp,origin='lower')
#    plt.title('disp: max={} min={} mean={}'.format(
#                scat_disp.max(),
#                scat_disp.min(),
#                scat_disp.mean()))
#
#    plt.figure(2)
#    plt.imshow(scat_disp-scat_prof,origin='lower')
#    plt.title('disp-prof: max={} min={} mean={}'.format(
#              (scat_disp-scat_prof).max(),
#              (scat_disp-scat_prof).min(),
#              (scat_disp-scat_prof).mean()))
#
#    plt.figure(3)
#    plt.imshow(im-scat_disp,origin='lower')
#    plt.title('im-disp: max={} min={} mean={}'.format(
#               (im-scat_disp).max(),
#               (im-scat_disp).min(),
#               (im-scat_disp).mean()))
#
#    plt.figure(4)
#    x=2000
#    plt.plot((im-scat_disp)[:,x:x+100].sum(axis=1)/100)
#    plt.title('im-disp')
#    plt.figure(5)
#    plt.plot((im-scat_prof)[:,x:x+100].sum(axis=1)/100)
#    plt.title('im-prof')
#
#    plt.show()
#
#
#def scatter(im_in, scatter_regions, plotalot=False,
#            prof_order=6, disp_order=11, disp_bin=64, dispFitst=False):
#    """Create a polynomial surface of the scattered light in the frame"""
#    msk=np.ones_like(im_in,dtype=np.bool)
#    for r0,r1 in scatter_regions:
#        msk[r0:r1+1,:]=False
#    #Useing a masked array doesn't seem to work with median_filter
#    #im=np.ma.array(im_in,mask=msk,fill_value=np.nan)
#    im=im_in.copy()
#    im[msk]=np.nan
#
#    if plotalot:
#        figure(1)
#        foo=im_in.copy()
#        foo[~msk]=np.nan
#        imshow(foo, origin='lower')
#        show()
#    nrow,ncol=im.shape
#
#    bin_col=np.arange(disp_bin/2, ncol, disp_bin)
#    rows=np.arange(nrow)
#    cols=np.arange(ncol)
#
#    binned=median_filter(im, size=(1,disp_bin) )
#    binned[msk]=np.nan
#    binned=binned[:,bin_col]
#
#    bin_poly_coeff=np.zeros((ncol/disp_bin, prof_order+1))
#    for i,column in enumerate(binned.T):
#        bin_poly_coeff[i,:]=np.polyfit(rows[~np.isnan(column)],
#                                       column[~np.isnan(column)],
#                                       prof_order)
#        if plotalot:
#            figure(2)
#            line,=plt.plot(rows, column)
#            plt.plot(rows, np.poly1d(bin_poly_coeff[i,:])(rows),
#                     color=line.get_color())
#
#    if plotalot:
#        show()
#
#
#    #Fit polynomials to each term
#    coeff_poly=[]
#    for i, coeff in enumerate(bin_poly_coeff.T):
#        coeff_poly.append(np.poly1d(np.polyfit(bin_col, coeff, disp_order)))
#
#    if plotalot:
#        figure(3)
#        fig, subplots = plt.subplots(3,3)
#        subplots=subplots.flatten()
#        for i,coeff in enumerate(bin_poly_coeff.T):
#            subplots[i].plot(bin_col,coeff)
#            #subplots[i].title('Coeff {}'.format(i))
#            subplots[i].plot(cols,coeff_poly[i](cols))
#        show()
#
#    #Use polynomials to generate the image
#    scatterim=np.zeros_like(im)
#    for j in range(ncol):
#        scatterim[:,j]=np.poly1d([poly(j) for poly in coeff_poly])(rows)
#
#    return scatterim
#
#def scatterxfirst(im_in, scatter_regions, plotalot=False, disp_bin=32,
#            prof_order=8, disp_order=11, prof_bin=1):
#    """Create a polynomial surface of the scattered light in the frame"""
#    
#    im=im_in.copy()
#
#    nrow,ncol=im.shape
#    rows=np.arange(nrow)
#    cols=np.arange(ncol)
#
#    #Create mask for scatterd light regions
#    #True=bad
#    msk=np.ones_like(im_in,dtype=np.bool)
#    row_msk=np.ones(nrow,dtype=np.bool)
#    for r0,r1 in scatter_regions:
#        row_msk[r0:r1+1]=False
#        msk[r0:r1+1,:]=False
#    
#    #Using a masked array doesn't seem to work with median_filter
#    #im=np.ma.array(im_in,mask=msk,fill_value=np.nan)
#    im[msk]=np.nan
#
#    #Show the regions not being used
#    if plotalot:
#        figure(1)
#        foo=im_in.copy()
#        foo[~msk]=np.nan
#        imshow(foo, origin='lower')
#        show()
#
#    #Configure for dispersion or profile fit first
#    xx=np.arange(disp_bin/2, ncol, disp_bin)
#    poly_xx=rows
#
#    poly_ndxs=rows[~row_msk]
#
#    ord=disp_order
#    poly_ord=prof_order
#
#    #Bin the image
#    binned=median_filter(im, size=(prof_bin, disp_bin) )
#    binned[msk]=np.nan
#
#    #Pixk out every disp_bin column starting at disp_bin/2
#    binned=binned[:,xx]
#
#
#    #Fit polynomials to the the scattered light along xx
#    #fitting every row in scattered region at binned columns
#    poly_coeff=np.zeros( (len(poly_ndxs), ord+1) )
#
#    for i,ndx in enumerate(poly_ndxs):
#        yy=binned[ndx]
#        poly_coeff[i,:]=np.polyfit(xx[np.isfinite(yy)],
#                                   yy[np.isfinite(yy)], ord)
#        if plotalot:
#           figure(2)
#           line,=plt.plot(xx, yy)
#           plt.plot(xx, np.poly1d(poly_coeff[i,:])(xx), color=line.get_color())
#
#    if plotalot:
#        show()
#
#    #Fit polynomials to each term of the polynomials found
#    polys=[]
#    for coeff in poly_coeff.T:
#        polys.append(np.poly1d(np.polyfit(poly_ndxs, coeff, poly_ord)))
#
#
#    if plotalot:
#        fig, subplots = plt.subplots(3,4)
#        subplots=subplots.flatten()
#        for i,coeff in enumerate(poly_coeff.T):
#            subplots[i].plot(poly_ndxs, coeff,'.')
#            subplots[i].plot(poly_xx, polys[i](poly_xx))
#        show()
#    
#    #Use polynomials to generate the image
#    scatterim=np.zeros_like(im)
#    for x in poly_xx:
#        scatterim[x,:]=np.poly1d([poly(x) for poly in polys])(cols)
#    
#    return scatterim



def make_glow_model(im_in, bin_x=1, bin_y=1):
    """ Model the amplifier glow """
    im=im_in.copy()
    im[0]=im[2]
    im[1]=im[2]
    im[-1]=im[-2]
    
    #glow image
    glow=np.zeros_like(im)
    
    #meshgrid
    x, y = np.meshgrid(np.arange(im.shape[1]), np.arange(im.shape[0]))
    
    
    def model_corner(im, x0, y0, xw, yw, iparams, std_clip=0):
        """ std_clip is the y height of the small corner to use to exclude
        spectra in the large corner,
        
        (iparams=(glow amp, x center, y center, xwid, ywid, xy amount)
        
        positions and initial params adjusted automatically for binning
        pass coordinates in 4k positions
        """
        x0/=bin_x
        y0/=bin_y
        xw/=bin_x
        yw/=bin_y
        iparams=list(iparams)
        iparams[1]/=bin_x
        iparams[2]/=bin_y
        iparams[3]/=bin_x
        iparams[4]/=bin_y
        
        corner=im[y0:y0+yw,x0:x0+xw].copy()
        if std_clip:
            small_corner=im[y0:y0+std_clip,x0:x0+xw].copy()
            patch_locs=corner>2*small_corner.std()
            patch_locs[:y0+std_clip,:]=False
            corner[patch_locs]=np.median(small_corner)
        cim, param= gaussfit2D(corner, iparams)
        param=list(param)
        param[-1]=0
        param[1]+=x0
        param[2]+=y0
        return gauss2D(( x,y), *param)
    
    #Lower R
    try:
        tmp=model_corner(im, 3996, 2, 100, 100,
                           (150, 58, -7, 30.0, 20.0, 0, 0))
        if tmp.min() < 0:
            raise RuntimeError('Glow model has negative values')
        else:
            glow+=tmp

    except RuntimeError, e:
        print 'Lower R glow model failed: {}'.format(str(e))

    #Lower L
    try:
        tmp=model_corner(im, 0, 2, 100, 100,
                           (150, 40, 0, 30.0, 20.0, 0, 0),
                           std_clip=50)
        if tmp.min() < 0:
            raise RuntimeError('Glow model has negative values')
        else:
            glow+=tmp

    except RuntimeError, e:
        print 'Lower L glow model failed: {}'.format(str(e))
    

    #Upper L
    try:
        tmp=model_corner(im, 0, 4012, 100, 100,
                           (150, 40, 100, 30.0, 20.0, 0, 0))
        if tmp.min() < 0:
            raise RuntimeError('Glow model has negative values')
        else:
            glow+=tmp

    except RuntimeError, e:
        print 'Upper L glow model failed: {}'.format(str(e))

    #Upper R
    try:
        tmp=model_corner(im, 3996, 4000, 100, 100,
                           (150, 58, 100, 30.0, 20.0, 0, 0))
        if tmp.min() < 0:
            raise RuntimeError('Glow model has negative values')
        else:
            glow+=tmp

    except RuntimeError, e:
        print 'Upper R glow model failed: {}'.format(str(e))
            
    return glow



def mkscatter(im_in, plot=False, scatter_thresh=.7, debug=False, header=None,
              do_glow=True, offset=0, glowOnly=False):
    """Make a scattered light map"""

    do_glow|=glowOnly

    im_in=im_in.copy()

    if debug:
        import ipdb;ipdb.set_trace()

    if header!=None:
        if header['BINNING'][0]!=header['BINNING'][-1]:
            raise ValueError('Rectangular binning, enable debug'
                             ' and see which of the lines is correct.')
        bin_x, bin_y=int(header['BINNING'][0]), int(header['BINNING'][-1])
    else:
        bin_x, bin_y=1, 1


    #model amp glow
    if do_glow:
        glow=make_glow_model(im_in, bin_x=bin_x, bin_y=bin_y)
    else:
        glow=np.zeros_like(im_in)


    #first two rows and last row are abberations
    im_in[0]=im_in[2]
    im_in[1]=im_in[2]
    im_in[-1]=im_in[-2]

    #regions to estimate the scattered light
    scatter_regions=[(0,50), (480,575), (1010,1100), (1525,1600), (2035,2120),
                     (2540,2620), (3030,3110), (3525, 3600), (3995 ,4109)]

    scatter_regions=np.array(scatter_regions)/bin_y

    scatter_regions[0,1]+=offset
    scatter_regions[-1,0]+=offset
    scatter_regions[1:-1,:]+=offset

    #Select all pixels in the scatter regions
    msk=np.ones_like(im_in, dtype=np.bool)
    for r0,r1 in scatter_regions: msk[r0:r1+1,:]=False
    s_pix=(im_in-glow)[~msk]


    #Compute a ~robust stddev for the scattered light
    s_std=s_pix[abs(s_pix-s_pix.mean())< 2*s_pix.std()].std()


    #from astropy.modeling import models, fitting
    #p_init = models.Polynomial2D(degree=5)
    #f = fitting.LinearLSQFitter()
    #y, x = np.mgrid[:im_masked.shape[0], :im_masked.shape[1]]
    #p = f(p_init, x, y, im_masked)

    #Create a masked image that includes all pixels within 1 scattered
    #light stddev of the mean scattered light level. Subtracting off the glow
    im_masked=np.ma.array(im_in-glow, copy=True,
                          mask=abs(im_in-glow-s_pix.mean())>
                                scatter_thresh*s_std)

    s_model=im_in.copy()

    nrow,ncol=im_in.shape
    rows=np.arange(nrow)
    cols=np.arange(ncol)

    #Model each column with a polynomial fit to the masked image
    for c in cols:
        s_model[:,c]=np.poly1d(np.ma.polyfit(rows, im_masked[:,c],8))(rows)

    if not glowOnly:
        #Smooth the model & add the glow to make the scattered light image
        s_im=gsmooth2d(s_model,(8,32)) + glow

        err=(im_masked-s_im).std()
    else:
        s_im=glow
        err=0.0

    #Debugging plots
    if plot or debug:
        plt.figure(1)
        cmap=plt.cm.jet
        cmap.set_bad('black')
        plt.imshow(im_masked, vmin=im_masked.min(),
                   vmax=im_masked.mean()+2*im_masked.std(),
                   cmap=cmap,origin='lower')
        plt.title('scatter data')
        plt.colorbar()
        
#        plt.figure(2)
#        plt.imshow(s_model, vmin=s_model.min(),
#                   vmax=s_model.mean()+2*s_model.std())
#        plt.title('c')
#        plt.colorbar()

        plt.figure(3)
        vmax_scatter=s_im.mean()+2*s_im.std()
        plt.imshow(s_im, vmin=s_im.min(), vmax=vmax_scatter, origin='lower')
        plt.title('model')
        plt.colorbar()
        
        plt.figure(4)
        clean=im_in-s_im
        plt.imshow(clean, vmin=0, vmax=vmax_scatter, origin='lower')
        plt.title('cleaned')
        plt.colorbar()
        
        for sr in scatter_regions:
            for x in sr:
                axhline(x, color='black')

        plt.show(block=False)

        if debug:
            import ipdb;ipdb.set_trace()

        return s_im, err,  (s_model, im_masked, glow)

    return s_im, err


