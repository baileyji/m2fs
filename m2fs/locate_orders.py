import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter1d

def locate_orders(filename, peakcol=2048, thresh=None,
                  bin=51, sigkern=1.5, nthreshsig=1.0):
    if isinstance(filename, (fits.hdu.hdulist.HDUList, str)):
        if type(filename)==str:
            try:
                hdul = fits.open(filename)
            except:
                print "Could not open filename %s" % filename
                return None
        else:
            hdul = filename
        im=hdul[0].data
        var=hdul[1].data
        mask=hdul[3].data
    else:
        #asssume ndarray
        im=filename.copy()
        var=im
        mask=np.zeros_like(im,dtype=np.uint8)


    nx,ny=im.shape


    #Take a swath at dispersion midpoint

    xlo = max(peakcol-bin/2, 0)
    xhi = min(peakcol+bin/2, nx-1)
    mid = np.median(im[:,xlo:xhi+1], axis=1)

    #Get approx Sigma of MID
    sig = np.median(im[~mask.astype(np.bool)])/np.sqrt(bin)  # sigma for MID

    #Convolve with a symmetric kernel
    mid2=gaussian_filter1d(mid, sigkern, mode='constant')

    # Find peaks, must be greater than neighbors and above threshold
    if not thresh:
        thresh=nthreshsig*sig

    peak , = np.where(np.r_[False, mid2[1:] > mid2[:-1]] &
                      np.r_[mid2[:-1] > mid2[1:], False] &
                      (mid2 > thresh))
    npeak=len(peak)

    print '{} peaks found'.format(npeak)
#    if npeak > 128:
#        import ipdb
#        ipdb.set_trace()

    return peak



