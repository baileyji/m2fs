#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
from scipy.optimize import curve_fit
from skimage.feature import peak_local_max
import argparse
import ipdb
import os.path

FITTING_TOL=1e-2

def parse_cl():
    parser = argparse.ArgumentParser(description='S/N Check',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
                                 
    parser.add_argument('file', metavar='FILE', type=str,
                     help='File to process')
    
    #Finding settings
    parser.add_argument('-d', dest='sep', default=10,
                        action='store', required=False, type=int,
                        help='Minimum delta between order peaks')
    parser.add_argument('-m', dest='min', default=50,
                        action='store', required=False, type=int,
                        help='Minimum peak intensity')
    parser.add_argument('-w', dest='xw', default=25,
                        action='store', required=False, type=int,
                        help='Width around peak to use for fitting gaussian')

    #Measurement selection settings
    parser.add_argument('-s', dest='hwid', default=8.0,
                        action='store', required=False, type=int,
                        help='Width to either side of smash')

    return parser.parse_args()


def gaussfit(xdata, ydata, initialp, ftol=1e-5, maxfev=5000):
    def gauss_quad(x, a0, a1, a2, a3, a4, a5):
        z = (x - a1)
        y = a0 * np.exp(-z**2 / 2*a2**2) + a3 + a4 * x + a5 * x**2
        return y

    from scipy.optimize import curve_fit
    parameters, covariance = curve_fit(gauss_quad, xdata, ydata, p0=initialp,
                                       ftol=ftol, maxfev=maxfev)

    return (parameters, gauss_quad(xdata, *parameters))

def find_peaks(im, min_sep=10, minv=25, maxv=5000):
    """find points to use for psf measurements"""
    points=peak_local_max(im, min_distance=min_sep,
                          threshold_abs=minv,
                          threshold_rel=0, indices=False)
    points[((im>maxv) & points)]=False
    
    return np.where(points)[0]

def fit_peaks(image, points, wx, xw=25, plot=False, maxfev=5000):
    """fit the each peak at each point"""
    ret={}
    for i,x in enumerate(points):

        try:
            if (x-xw/2 <0 or x+xw/2 >= image.size): continue

            xx=np.arange(x-xw/2,x+xw/2+1)
            sim=image[x-xw/2:x+xw/2+1]
            offset=sim.min()
            fit=gaussfit(xx, sim,
                         (image[x]-offset, xw/2, wx/2.35482, offset, 0, 0),
                         maxfev=maxfev, ftol=FITTING_TOL)

            ret[x]={'param':fit[0],'model':fit[1]}
            plt.plot(xx, fit[1],'b')

        except RuntimeError as e:
            print 'Point {} failed. {}'.format(i, e)
        except (ValueError, IndexError) as e:
            import ipdb;ipdb.set_trace()

    return ret

def compute_sn(file, args):
    
    #Settings
    min_sep=args.sep
    minv=args.min
    xw=args.xw #Width of window around each peak to model
    smashw=args.hwid
    maxv=1e99
    
    #Fitting settings
    fwhm_guess=xw/2

    #fetch image
    hdu=fits.open(file)
    header=hdu[0].header
    side= os.path.basename(file)[0].lower()
    if side not in 'rb': side='?'

    binning=int(header['BINNING'][2])
    min_sep/=binning
    xw/=binning
    smashw/=binning

    fim=hdu['SCIENCE'].data
    vfim=hdu['VARIANCE'].data
    
    ncol=fim.shape[1]
    c0=ncol/2 - 2*.02*ncol
    c1=ncol/2 - .02*ncol
    im=fim[:,c0:c1].mean(1)
    vim=vfim[:,c0:c1].mean(1)

    #Find the peaks
    peaks=find_peaks(im, min_sep=min_sep, minv=minv, maxv=maxv)

    print('Found {} peaks between {} and {} separated by at least {}.'.format(
          len(peaks), minv, maxv, min_sep))


    #Plot image and detected peaks
    plt.clf()
    plt.plot(im)
    for x in peaks: plt.plot(x,im[x],'w*')

    #Fit the peaks
#    ret=fit_peaks(im, peaks, fwhm_guess, xw, plot=args.plot or args.debug)


    #Smash peaks
    sns=[]
    for i,x in enumerate(peaks):
#        cent=ret[x]['param'][1]
        cent=x
        sl=slice(cent-smashw,(cent+1)+smashw)
        sn=im[sl].sum()/np.sqrt(vim[sl].sum())
        sns.append(sn)
        
        plt.text(x ,1.01*im[x], '{:.0f}'.format(sn))

        plt.gca().add_patch(matplotlib.patches.Rectangle((x-smashw,0),
                                                   2*smashw+1, im[x],
                                                   alpha=.2))


    plt.title('Average of {:.0f} pixels at at {:.0f}'.format(c1-c0,.5*(c0+c1))+
              '  Median S/N {:.0f}'.format(np.median(sns)))
    plt.xlabel('Row (R1-01/B8-01 to left)')
    plt.ylabel('electrons')
    plt.axhline(minv)
    plt.xlim(0,fim.shape[0])
    plt.show(0)
    print 'Median S/N: {:.0f}'.format(np.median(sns))

#    if args.debug:
#        input=raw_input('Continue(any), Abort(a), Debug(db)>')
#        if input=='db':
#           ipdb.set_trace()
#        elif input=='a':
#           exit(1)



if __name__ =='__main__':
    
    import matplotlib as mpl
    mpl.rcParams['font.size']=16
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['ytick.major.size']=6
    mpl.rcParams['ytick.major.width']=2
    mpl.rcParams['xtick.major.size']=6
    mpl.rcParams['xtick.major.width']=2

    args=parse_cl()
    
    compute_sn(args.file,args)

    raw_input('Any key to exit')
