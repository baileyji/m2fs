#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from jbastro.astroLib import gauss2D, gaussfit2D, gauss2dmodel
from skimage.feature import peak_local_max
import argparse
from glob import glob
import ipdb
import os.path

FITTING_TOL=1e-2

def parse_cl():
    parser = argparse.ArgumentParser(description='Focus Series Analyzer',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #Finding settings
    parser.add_argument('-sep', dest='min_sep', default=10,
                        action='store', required=False, type=int,
                        help='Spot separation threshold')
    parser.add_argument('-min', dest='minv', default=50,
                        action='store', required=False, type=int,
                        help='Minimum spot intensity')
    parser.add_argument('-max', dest='maxv', default=2000,
                        action='store', required=False, type=int,
                        help='Maximum spot intensity')
    parser.add_argument('-npeak', dest='max_peaks', default=4000,
                        action='store', required=False, type=int,
                        help='Maximum number of points to measure')
    #Fitting Settings
    parser.add_argument('-xwin', dest='xw', default=25,
                        action='store', required=False, type=int,
                        help='width of x fitting window')
    parser.add_argument('-ywin', dest='yw', default=25,
                        action='store', required=False, type=int,
                        help='Width of y fitting window')
    parser.add_argument('-x_wid', dest='x_wid_guess', default=4,
                        action='store', required=False, type=float,
                        help='PSF FWHM x guess')
    parser.add_argument('-y_wid', dest='y_wid_guess', default=4,
                        action='store', required=False, type=float,
                        help='PSF FWHM y guess')

    #Measurement selection settings
    parser.add_argument('-x_min', dest='x_min', default=1.0,
                        action='store', required=False, type=float,
                        help='Minimum PSFx')
    parser.add_argument('-x_max', dest='x_max', default=10.0,
                        action='store', required=False, type=float,
                        help='Maximum PSFx')
    parser.add_argument('-y_min', dest='y_min', default=1.0,
                        action='store', required=False, type=float,
                        help='Minimum PSFx')
    parser.add_argument('-y_max', dest='y_max', default=10.0,
                        action='store', required=False, type=float,
                        help='Maximum PSFy')

    #Program Settings
    parser.add_argument('-f', dest='files', default='',
                        action='store', required=False, type=str,
                        help='Input files')
    parser.add_argument('-l', dest='listfile', default='',
                        action='store', required=False, type=str,
                        help='Input listfile')
    parser.add_argument('-d', dest='dir', default='./',
                        action='store', required=False, type=str,
                        help='Input file directory')
    parser.add_argument('--plot', dest='plot', default=False,
                        action='store_true', required=False,
                        help='Show many plots')
    parser.add_argument('--debug', dest='debug', default=False,
                        action='store_true', required=False,
                        help='Debugging prompts')
                        
    return parser.parse_args()

def derangify(s):
    """
    Takes a range in form of "a-b" and generate a list of numbers between 
    a and b inclusive.
    Also accepts comma separated ranges like "a-b,c-d,f" will build a 
    list which will include
    Numbers from a to b, a to d and f
    http://code.activestate.com/recipes/577279-generate-list-of-
    numbers-from-hyphenated-and-comma/
    """
    s="".join(s.split())#removes white space
    r=set()
    for x in s.split(','):
        t=x.split('-')
        if len(t) not in [1,2]:
            raise SyntaxError("hash_range is given its "
                              "arguement as "+s+" which seems not "
                              "correctly formated.")
        if len(t)==1:
            r.add(int(t[0]))
        else:
            r.update(set(range(int(t[0]),int(t[1])+1)))
    l=list(r)
    l.sort()
    return tuple(l)

def gauss2dmodel(xw, yw, amp, xo, yo, sigma_x, sigma_y, covar, offset):
    xo = float(xo)
    yo = float(yo)
    
    x = np.arange(xw+1,dtype=np.float)-xw/2 + xo
    y = np.arange(yw+1,dtype=np.float)-yw/2 + yo
    x, y = np.meshgrid(x, y)
    
    rho=covar/(sigma_x*sigma_y)
    
    z=(((x-xo)/sigma_x)**2 - 2*rho*(x-xo)*(y-yo)/(sigma_x/sigma_y) +
       ((y-yo)/sigma_y)**2)
    
    g=amp*np.exp(-z/(2*(1-rho**2)))+ offset
    
    return g

def gauss2D((x, y), amp, xo, yo, sigma_x, sigma_y, covar, offset):
    xo = float(xo)
    yo = float(yo)
    
    rho=covar/(sigma_x*sigma_y)
    
    z=(((x-xo)/sigma_x)**2 - 2*rho*(x-xo)*(y-yo)/(sigma_x/sigma_y) +
        ((y-yo)/sigma_y)**2)
    
    g=amp*np.exp(-z/(2*(1-rho**2)))+ offset

    return g

def gaussfit2D(im, initialp, ftol=1e-5, maxfev=5000):
    def g2d((x, y), amp, xo, yo, sigma_x, sigma_y, covar, offset):
        return gauss2D((x, y), amp, xo, yo, sigma_x,
                       sigma_y, covar, offset).ravel()
    
    x = np.arange(im.shape[0],dtype=np.float)
    y = np.arange(im.shape[1],dtype=np.float)
    x, y = np.meshgrid(x, y)
    

    popt, pcov = curve_fit(g2d, (x, y), im.ravel(), p0=initialp,
                           ftol=ftol, maxfev=maxfev)
    model=gauss2D((x, y), *popt)
                           
    return model, popt

def find_peaks(im, min_sep=10, minv=25, maxv=5000):
    """find points to use for psf measurements"""
    points=peak_local_max(im, min_distance=min_sep,
                          threshold_abs=minv,
                          threshold_rel=0, indices=False)
    points[((im>maxv) & points)]=False
    
    return zip(*np.where(points))

def filter_peaks(im, ret, x_min=2.0, x_max=15.0, y_min=2.0, y_max=15.0,
                 xw=25, yw=25):
    """
    recover data after fitting peaks and make cuts on it
        
    returns x,y,xw,yw,covar, modelim, modeldat
    """
    
    modeldat=np.zeros(im.shape)
    modelim=np.zeros(im.shape)
    x,y=[],[]
    xwid,ywid,covar=[],[],[]
    for peak,(psf,mod) in ret.items():
        yy,xx=peak
        xsli=slice(xx-xw/2, xx+xw/2+1)
        ysli=slice(yy-yw/2, yy+yw/2+1)
        modeldat[ysli, xsli]=im[ysli, xsli]
        modelim[ysli, xsli]=mod
        x.append(xx)
        y.append(yy)
        xwid.append(psf[0])
        ywid.append(psf[1])
        covar.append(psf[2])
    
    x=np.array(x)
    y=np.array(y)
    xwid=np.array(xwid)
    ywid=np.array(ywid)
    covar=np.array(covar)
    
    #Filter the data
    keep=(xwid<=x_max) & (xwid>=x_min) & (ywid>=y_min) & (ywid<=y_max)
    x=x[keep]
    y=y[keep]
    xwid=xwid[keep]
    ywid=ywid[keep]
    covar=covar[keep]
    
    return x,y,xwid,ywid,covar, modelim, modeldat

def fit_peaks(image, points, wx, wy, xw=25, yw=25, covar=.05, debug=False,
              maxfev=5000):
    """fit the psf at each point"""
    ret={}
    for i,(y,x) in enumerate(points):

        try:
            if (y-yw/2 <0 or y+yw/2 >= image.shape[0] or
                x-xw/2 <0 or x+xw/2 >= image.shape[1]):
                continue

            sim=image[y-yw/2:y+yw/2+1,x-xw/2:x+xw/2+1]
            offset=sim.min()
            fit=gaussfit2D(sim,
                           (image[y,x]-offset,
                            yw/2, xw/2,
                            wy/2.35482, wx/2.35482,
                            covar,
                            offset),maxfev=maxfev, ftol=FITTING_TOL)
            psf=np.array(fit[1][3:6])
            psf[0:2]*=2.35482
            model=fit[0]
            ret[(y,x)]=psf,model

            if debug:
                plt.figure(10)
                plt.subplot(1,2,1)
                plt.imshow(sim, vmin=sim.min(), vmax=image[y,x],
                           interpolation='nearest')
                plt.title('{},{}'.format(x,y))
                plt.subplot(1,2,2)
                plt.imshow(model, vmin=sim.min(), vmax=image[y,x],
                           interpolation='nearest')
                plt.title('FWHMx: {:.1f} FWHMy:{:.1f} Covar:{:.1f}'.format(*psf))
                plt.show(0)
                input=raw_input('Debug peak fitter? (n/a/db)>')
                if input=='n':
                    debug=False
                elif input=='a':
                    return ret
                elif input=='db':
                    import ipdb;ipdb.set_trace()

        except RuntimeError as e:
            print 'Point {} failed. {}'.format(i, e)
        except (ValueError, IndexError) as e:
            import ipdb;ipdb.set_trace()
    
    return ret

def find_focus(files, args):
    
    #Finding settings
    min_sep=args.min_sep
    minv=args.minv
    maxv=args.maxv
    max_peaks=args.max_peaks

    #Fitting settings
    xw, yw=args.xw, args.yw #Width of window around each line to model
    x_wid_guess=args.x_wid_guess #FWHM
    y_wid_guess=args.y_wid_guess

    #Fit cut settings
    x_min,x_max=args.x_min,args.x_max
    y_min,y_max=args.y_min,args.y_max

    focus_data={}
    for f in files:

        print('Processing {}'.format(f))

        #fetch image
        header=fits.getheader(f)
        focus=header['FOCUS']
        filt=header['FILTER']
        if os.path.basename(f)[0].lower()=='b':
            side='b'
            tempi=[1,3,5,7,9]
        else:
            side='r'
            tempi=[1,2,4,6,8]
        temps=tuple(header['TEMP{:02}'.format(i)] for i in tempi)
        im=fits.getdata(f)

        #Find the peaks
        peaks=find_peaks(im, min_sep=min_sep, minv=minv, maxv=maxv)

        print('Found {} peaks between {} and {} separated by at least {}.'.format(
              len(peaks), minv, maxv, min_sep))

        if len(peaks) > max_peaks:
            np.random.seed(1)
            peaks=np.array(peaks)[np.random.randint(0, len(peaks), max_peaks)].tolist()

        #Plot image and detected peaks
        if args.plot or args.debug:
            plt.figure(1)
            plt.clf()
            plt.imshow(im, vmin=10, vmax=100, origin='lower')
            for y,x in peaks: plt.plot(x,y,'w*')
            plt.title('Image & Peaks')
            plt.xlim(0,im.shape[1])
            plt.ylim(0,im.shape[0])
            plt.show(0)


        #Fit the peaks
        ret=fit_peaks(im, peaks, x_wid_guess, y_wid_guess, xw=xw, yw=yw,
                      debug=False)#args.debug)

        #Filter fits and construct results
        x, y, xwid, ywid, covar, modelim, modeldat=filter_peaks(im, ret,
                      x_min=x_min, x_max=x_max,
                      y_min=y_min, y_max=y_max,
                      xw=xw, yw=yw)

        #Plot fit data & models
        if args.plot or args.debug:
            plt.figure(2)
            plt.clf()
            plt.subplot(1,2,1)
            plt.imshow(modelim, vmin=0, vmax=100, origin='lower')
            plt.title('Line fits')
            plt.xlim(800,1200)
            plt.ylim(0,375)
            plt.subplot(1,2,2)
            plt.imshow(modeldat, vmin=0, vmax=100, origin='lower')
            plt.title('Line fit data')
            plt.xlim(800,1200)
            plt.ylim(0,375)
            plt.show(0)


            #Plot width data
            plt.figure(3)
            plt.subplot(1,3,1)
            plt.hist(xwid,100,histtype='stepfilled')
            plt.title('X FWHM: Avg {:.2f} Focus:{}'.format(xwid.mean(), focus))
            plt.xlim(0,15)
            plt.subplot(1,3,2)
            plt.hist(ywid,100,histtype='stepfilled')
            plt.title('Y FWHM: Avg {:.2f}'.format(ywid.mean()))
            plt.xlim(0,15)
            plt.subplot(1,3,3)
            plt.hist(covar,100,histtype='stepfilled')
            plt.title('Covar Avg {:.2f}'.format(covar.mean()))
            plt.show(0)
            plt.draw()

        if args.debug:
            input=raw_input('Continue(any), Abort(a), Debug(db)>')
            if input=='db':
               ipdb.set_trace()
            elif input=='a':
               exit(1)

#        print('Focus: {}. FWHM {:.1f}, {:.1f}. Covariance {:.1f}'.format(
#               focus, xwid.mean(), ywid.mean(),covar.mean()))
        focus_data[focus]=(xwid.mean(), ywid.mean(),covar.mean(),temps)


    foc=focus_data.keys()
    xvals=[focus_data[f][0] for f in foc]
    yvals=[focus_data[f][1] for f in foc]
    cvals=[focus_data[f][2] for f in foc]
    
    if len(foc) >2:
        cx=np.polyfit(foc, xvals, 2)
        cy=np.polyfit(foc, yvals, 2)
        min_x=-cx[1]/cx[0]/2
        min_y=-cy[1]/cy[0]/2
        print 'Processed {}'.format(files)
        print('Filter: {}'.format(filt))
        print('Best x focus @ {:.1f} with value of {:.2f}'.format(min_x,
              np.poly1d(cx)(min_x)))
        print('Best y focus @ {:.1f} with value of {:.2f}'.format(min_y,
              np.poly1d(cy)(min_y)))
        print '{:.1f} {:.1f} {} {} {} {} {}'.format(min_x, min_y,
                                                    *focus_data[f][3])
    for i,f in enumerate(foc):
        print('Focus: {} FWHM: {:.1f}, {:.1f} Covar: {:.1f} Temps:{}'.format(
              foc[i],xvals[i],yvals[i],cvals[i], focus_data[f][3]))

    temps=np.array([k[3][1] for k in focus_data.values()])

    if len(foc) >2:
        plt.figure(6)
        plt.plot(foc, xvals, 'bo', label='PSF x')
        plt.plot(foc, yvals, 'ro', label='PSF y')
        plt.plot(foc, cvals, 'go', label='Covar')

        extent=np.abs(np.array(foc)-min_x).max()
        xx=np.linspace(min(min(foc),min_x-extent), max(max(foc),min_x+extent),
                       100)
        plt.plot(xx, np.poly1d(cx)(xx),'b')
        plt.plot(xx, np.poly1d(cy)(xx),'r')
        plt.text(min_x,1.5,
                 'Temp: {:.3f}\nFilter: {}'.format(temps.mean(),filt),
                 color=side)
        plt.legend()
        plt.show(0)


    if args.debug:
        ipdb.set_trace()

def get_seqnos(listfile):
    ret=set()
    with open(listfile,'r') as lf:
        for l in lf:
            if l[0] in '1234567890':
                range=l.split()[0]
                ret.update(map(str,derangify(range)))
            elif len(l)>1 and l[0] in 'RBrb' and l[1] in '1234567890':
                range=l[1:].split()[0]
                ret.update(map(lambda x: l[0].lower()+str(x), derangify(range)))
    return list(ret)

if __name__ =='__main__':
    
    args=parse_cl()
    
    if not args.files:
        files=(glob(os.path.join(args.dir,'*.fits'))+
               glob(os.path.join(args.dir,'*.fits.gz')))
    
        try:
            seqno=get_seqnos(args.listfile)
            info=[m2fs.obs.info(f,no_load=True) for f in files]

            files=[]
            for i in info:
                for s in seqno:
                    if i.seqno_match(s):
                        files.append(i.file)
                        break
        except IOError:
            print 'No listfile, doing all'
    else:
        files=[os.path.join(args.dir, x) for x in args.files.split(' ')]

    print('Running on {}'.format(files))

    rfiles=[f for f in files if os.path.basename(f)[0] == 'r']
    bfiles=[f for f in files if os.path.basename(f)[0] == 'b']
    if rfiles: find_focus(rfiles,args)
    if bfiles: find_focus(bfiles,args)

    raw_input('Any key to exit')

    #from scipy.interpolate import griddata
    #yi,xi=np.mgrid[0:im.shape[0], 0:im.shape[1]]
    #xwidi=griddata((x,y), xwid, (xi, yi), method='nearest')
    #ywidi=griddata((x,y), ywid, (xi, yi), method='nearest')
    #covari=griddata((x,y), covar, (xi, yi), method='nearest')
    #plt.imshow(xwidi)
    #
    #
    #for xx,yy, w in zip(x,y,xwid):
    #    gcf().gca().add_artist(plt.Circle((xx,yy), w*10,fill=False, color='w'))
    #
    #
    #
    #
    #plt.figure(3)
    #plt.clf()
    #plt.subplot(1,3,1)
    #plt.imshow(xwidi, origin='lower')
    #plt.title('X Width')
    #plt.colorbar()
    #plt.subplot(1,3,2)
    #plt.imshow(ywidi, origin='lower')
    #plt.title('Y Width')
    #plt.colorbar()
    #plt.subplot(1,3,3)
    #plt.imshow(covari, origin='lower')
    #plt.title('Width Covariance')
    #plt.colorbar()

