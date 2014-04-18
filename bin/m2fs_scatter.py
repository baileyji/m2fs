#!/usr/bin/python
import scipy.interpolate
import numpy as np
from astropy.stats import sigma_clip
from astropy.io import fits
from m2fs.obs import scatterxfirst
import matplotlib.pyplot as plt

if __name__ == '__main__':
    import sys
    file=sys.argv[1]
    if len(sys.argv) >2:
        ofile=sys.argv[2]
    else:
        print ' Output file not specified'
        raise ValueError

    hdul=fits.open(file)

    im=hdul[0].data

    scatter_regions=[(2,75), (480,575), (1010,1100), (1525,1600), (2035,2120),
                     (2540,2620), (3030,3110), (3525, 3600), (3995 ,4109)]

    scat_disp=scatterxfirst(im,scatter_regions,prof_order=8)

    plt.figure(1)
    plt.imshow(scat_disp,origin='lower')
    plt.title('disp: max={} min={} mean={}'.format(
                scat_disp.max(),
                scat_disp.min(),
                scat_disp.mean()))


    plt.figure(2)
    plt.imshow(im-scat_disp,origin='lower')
    plt.title('im-disp: max={} min={} mean={}'.format(
               (im-scat_disp).max(),
               (im-scat_disp).min(),
               (im-scat_disp).mean()))

    plt.figure(3)
    x=2000
    plt.plot((im-scat_disp)[:,x:x+100].sum(axis=1)/100)
    plt.title('im-disp')

    plt.show()

    im-=scat_disp
    if ofile==file:
        hdul.append(fits.ImageHDU(scat_disp.astype(np.float64), name='scatter'))
        hdul.writeto(ofile)
    else:
        hdu=fits.PrimaryHDU(data=scat_disp)
        hdu.writeto(ofile)
