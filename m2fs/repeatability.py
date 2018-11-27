import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from jbastro.astroLib import gauss2D, gaussfit2D, gauss2dmodel
import argparse
from glob import glob
import ipdb
import os.path

FITTING_TOL=1e-2

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
    for peak,param,mod in ret:
        psf=param[2:]
        yy,xx=param[:2]
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
    ret=[]
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

            param=np.array(fit[1][1:6])
            param[2:4]*=2.35482
            param[:2]+=(y-yw/2,x-xw/2)
            model=fit[0]
            ret.append(((y,x),param,model))

            if debug:
                plt.figure(10)
                plt.subplot(1,2,1)
                plt.imshow(sim, vmin=sim.min(), vmax=image[y,x],
                           interpolation='nearest')
                plt.title('{},{}'.format(x,y))
                plt.subplot(1,2,2)
                plt.imshow(model, vmin=sim.min(), vmax=image[y,x],
                           interpolation='nearest')
                plt.title('FWHMx: {:.1f} FWHMy:{:.1f} Covar:{:.1f}'.format(*param[2:]))
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



def refine_peaks(im, peaks0, recenter=False,wid=9):
    """
    im, list of (x,y), wether to refine the xy to the nearest local maximum
    """
    min_sep=10 #minimum meak seperation
    num_peaks=1 #num to find at each input point
    minv=40 #minimum intensity
    maxv=55000 #maximum intensity
    
     #width around the input point to search

    xw, yw=25,25 #Width of window around each line to model

    x_wid_guess,y_wid_guess=4,4 #FWHM guess
    
    #min/max gaussin widths
    x_min,y_min=1.75,3
    x_max,y_max=10,10
    
    if recenter:
        peaks=[]
        plt.imshow(im,vmax=np.median(im)*3)
        ax=plt.gca()
        for pt in peaks0:
            x0=max(pt[0]-wid/2,0)
            y0=max(pt[1]-wid/2,0)
            subim=im[y0:pt[1]+wid/2+1,x0:pt[0]+wid/2+1]
            amax=np.unravel_index(subim.argmax(),subim.shape)
            if subim[amax]<maxv and subim[amax]>minv:
                npt=tuple((np.array(amax)+(y0,x0)).ravel().tolist())
                peaks.append(npt)
                ax.add_patch(mpl.patches.Circle(pt,10,facecolor='None',
                                    edgecolor='k',linewidth=4))
                ax.add_patch(mpl.patches.Arrow(pt[0],pt[1],npt[0]-pt[1],npt[1]-pt[0],
                                               color='r',linewidth=1))

    else:
        peaks=[x[::-1] for x in peaks0]

    #Fit the peaks
    ret=fit_peaks(im, peaks, x_wid_guess, y_wid_guess, xw=xw, yw=yw,
                  debug=True)#args.debug)

    #Filter fits and construct results
    x, y, xwid, ywid, covar, modelim, modeldat=filter_peaks(im, ret,
                  x_min=x_min, x_max=x_max,
                  y_min=y_min, y_max=y_max,
                  xw=xw, yw=yw)

    return x, y, xwid, ywid, covar, modelim, modeldat



dat=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0380.fits')
dat2=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0383.fits')
radius=10
markerlinewid=4

#step one plot and click a bunch

fig=plt.figure(figsize=(20.5,10))
plt.subplots_adjust(.01,.01,.99,.99)
refax=plt.subplot(1,2,1)

plt.imshow(dat,vmin=0,vmax=np.median(dat)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

ax2=plt.subplot(1,2,2,sharex=refax,sharey=refax)
plt.imshow(dat2,vmin=0,vmax=np.median(dat2)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

#now attach a picker
clicks=[]
def _click(event):
    tb = plt.get_current_fig_manager().toolbar
    global clicks
    if event.inaxes and tb.mode == '':
        clicks.append((event.xdata,event.ydata))
fig.canvas.mpl_connect('button_press_event', _click)


#step two plot circles where the clicks are in both plots
circles=[refax.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
circles=[ax2.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
plt.draw()


#Step three run the core of the focus program at those points
peaks0=list(clicks)
im=dat.copy()
wid=7
recenter=True
peaks_dat1=refine_peaks(dat,peaks0,recenter=True)
peaks1=np.array([x for x in zip(peaks_dat1[0],peaks_dat1[1])])

peaks_dat2=refine_peaks(dat2, peaks1,recenter=True)
peaks2=np.array([x for x in zip(peaks_dat2[0],peaks_dat2[1])])

#Now look at home much things shifted
#If the shifts were large then the second round of fits might not have worked
print 'Shift from clicks to 1 {:.2f} pix, clicks to 2 {:.2f} pix, 1 to 2 {:.2f} pix'.format(np.sqrt(((peaks1-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean())


print 'Y shift {:.2f}pm{:.2f} pix, X {:.2f}pm{:.2f} pix, Net shift {:.2f}pm{:.2f} pix'.format((peaks2-peaks1)[:,1].mean(),(peaks2-peaks1)[:,1].std(),(peaks2-peaks1)[:,0].mean(),(peaks2-peaks1)[:,0].std(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean(),np.sqrt(((peaks2-peaks1)**2).sum(1)).std())

"""
45-180-45
Y shift 0.37pm1.48 pix, X -0.44pm1.42 pix, Net shift 1.71pm1.27 pix
"""



dat=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0383.fits')
dat2=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0385.fits')
radius=10
markerlinewid=4

#step one plot and click a bunch

fig=plt.figure(figsize=(20.5,10))
plt.subplots_adjust(.01,.01,.99,.99)
refax=plt.subplot(1,2,1)

plt.imshow(dat,vmin=0,vmax=np.median(dat)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

ax2=plt.subplot(1,2,2,sharex=refax,sharey=refax)
plt.imshow(dat2,vmin=0,vmax=np.median(dat2)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

#step two plot circles where the clicks are in both plots
circles=[refax.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
circles=[ax2.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
plt.draw()


#Step three run the core of the focus program at those points
peaks0=list(clicks)
im=dat.copy()
wid=7
recenter=True
peaks_dat1=refine_peaks(dat,peaks0,recenter=True)
peaks1=np.array([x for x in zip(peaks_dat1[0],peaks_dat1[1])])

peaks_dat2=refine_peaks(dat2, peaks1,recenter=True)
peaks2=np.array([x for x in zip(peaks_dat2[0],peaks_dat2[1])])

#Now look at home much things shifted
#If the shifts were large then the second round of fits might not have worked
print 'Shift from clicks to 1 {:.2f} pix, clicks to 2 {:.2f} pix, 1 to 2 {:.2f} pix'.format(np.sqrt(((peaks1-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean())


print 'Y shift {:.2f}pm{:.2f} pix, X {:.2f}pm{:.2f} pix, Net shift {:.2f}pm{:.2f} pix'.format((peaks2-peaks1)[:,1].mean(),(peaks2-peaks1)[:,1].std(),(peaks2-peaks1)[:,0].mean(),(peaks2-peaks1)[:,0].std(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean(),np.sqrt(((peaks2-peaks1)**2).sum(1)).std())
"""
hires-lores-hires
Y shift -0.37pm1.26 pix, X 0.03pm1.22 pix, Net shift 1.38pm1.15 pix
"""



dat=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0385.fits')
dat2=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0388.fits')
radius=10
markerlinewid=4

#step one plot and click a bunch

fig=plt.figure(figsize=(20.5,10))
plt.subplots_adjust(.01,.01,.99,.99)
refax=plt.subplot(1,2,1)

plt.imshow(dat,vmin=0,vmax=np.median(dat)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

ax2=plt.subplot(1,2,2,sharex=refax,sharey=refax)
plt.imshow(dat2,vmin=0,vmax=np.median(dat2)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

#step two plot circles where the clicks are in both plots
circles=[refax.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
circles=[ax2.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
plt.draw()


#Step three run the core of the focus program at those points
peaks0=list(clicks)
im=dat.copy()
wid=7
recenter=True
peaks_dat1=refine_peaks(dat,peaks0,recenter=True)
peaks1=np.array([x for x in zip(peaks_dat1[0],peaks_dat1[1])])

peaks_dat2=refine_peaks(dat2, peaks1,recenter=True)
peaks2=np.array([x for x in zip(peaks_dat2[0],peaks_dat2[1])])

#Now look at home much things shifted
#If the shifts were large then the second round of fits might not have worked
print 'Shift from clicks to 1 {:.2f} pix, clicks to 2 {:.2f} pix, 1 to 2 {:.2f} pix'.format(np.sqrt(((peaks1-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean())


print 'Y shift {:.2f}pm{:.2f} pix, X {:.2f}pm{:.2f} pix, Net shift {:.2f}pm{:.2f} pix'.format((peaks2-peaks1)[:,1].mean(),(peaks2-peaks1)[:,1].std(),(peaks2-peaks1)[:,0].mean(),(peaks2-peaks1)[:,0].std(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean(),np.sqrt(((peaks2-peaks1)**2).sum(1)).std())
"""
hiaz by -.7 and back
Y shift -0.46pm1.59 pix, X -0.45pm1.72 pix, Net shift 2.13pm1.17 pix
"""




dat=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0388.fits')
dat2=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0393.fits')
radius=10
markerlinewid=4

#step one plot and click a bunch

fig=plt.figure(figsize=(20.5,10))
plt.subplots_adjust(.01,.01,.99,.99)
refax=plt.subplot(1,2,1)

plt.imshow(dat,vmin=0,vmax=np.median(dat)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

ax2=plt.subplot(1,2,2,sharex=refax,sharey=refax)
plt.imshow(dat2,vmin=0,vmax=np.median(dat2)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

#step two plot circles where the clicks are in both plots
circles=[refax.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
circles=[ax2.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
plt.draw()


#Step three run the core of the focus program at those points
peaks0=list(clicks)
im=dat.copy()
wid=7
recenter=True
peaks_dat1=refine_peaks(dat,peaks0,recenter=True)
peaks1=np.array([x for x in zip(peaks_dat1[0],peaks_dat1[1])])

peaks_dat2=refine_peaks(dat2, peaks1,recenter=True)
peaks2=np.array([x for x in zip(peaks_dat2[0],peaks_dat2[1])])

#Now look at home much things shifted
#If the shifts were large then the second round of fits might not have worked
print 'Shift from clicks to 1 {:.2f} pix, clicks to 2 {:.2f} pix, 1 to 2 {:.2f} pix'.format(np.sqrt(((peaks1-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean())


print 'Y shift {:.2f}pm{:.2f} pix, X {:.2f}pm{:.2f} pix, Net shift {:.2f}pm{:.2f} pix'.format((peaks2-peaks1)[:,1].mean(),(peaks2-peaks1)[:,1].std(),(peaks2-peaks1)[:,0].mean(),(peaks2-peaks1)[:,0].std(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean(),np.sqrt(((peaks2-peaks1)**2).sum(1)).std())
"""
hiel by -.2 and back
Y shift 0.16pm1.57 pix, X 0.16pm1.48 pix, Net shift 1.77pm1.26 pix
"""




dat=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0393.fits')
dat2=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0394.fits')
radius=10
markerlinewid=4

#step one plot and click a bunch

fig=plt.figure(figsize=(20.5,10))
plt.subplots_adjust(.01,.01,.99,.99)
refax=plt.subplot(1,2,1)

plt.imshow(dat,vmin=0,vmax=np.median(dat)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

ax2=plt.subplot(1,2,2,sharex=refax,sharey=refax)
plt.imshow(dat2,vmin=0,vmax=np.median(dat2)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

#step two plot circles where the clicks are in both plots
circles=[refax.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
circles=[ax2.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
plt.draw()


#Step three run the core of the focus program at those points
peaks0=list(clicks)
im=dat.copy()
wid=7
recenter=True
peaks_dat1=refine_peaks(dat,peaks0,recenter=True)
peaks1=np.array([x for x in zip(peaks_dat1[0],peaks_dat1[1])])

peaks_dat2=refine_peaks(dat2, peaks1,recenter=True)
peaks2=np.array([x for x in zip(peaks_dat2[0],peaks_dat2[1])])

#Now look at home much things shifted
#If the shifts were large then the second round of fits might not have worked
print 'Shift from clicks to 1 {:.2f} pix, clicks to 2 {:.2f} pix, 1 to 2 {:.2f} pix'.format(np.sqrt(((peaks1-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean())


print 'Y shift {:.2f}pm{:.2f} pix, X {:.2f}pm{:.2f} pix, Net shift {:.2f}pm{:.2f} pix'.format((peaks2-peaks1)[:,1].mean(),(peaks2-peaks1)[:,1].std(),(peaks2-peaks1)[:,0].mean(),(peaks2-peaks1)[:,0].std(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean(),np.sqrt(((peaks2-peaks1)**2).sum(1)).std())
"""
swap filter and back
Y shift -0.24pm1.54 pix, X 0.19pm1.39 pix, Net shift 1.86pm0.97 pix
"""



dat=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0379.fits')
dat2=fits.getdata('/Users/one/obsdata/sep2015/repeatability/b0380.fits')
radius=10
markerlinewid=4

#step one plot and click a bunch

fig=plt.figure(figsize=(20.5,10))
plt.subplots_adjust(.01,.01,.99,.99)
refax=plt.subplot(1,2,1)

plt.imshow(dat,vmin=0,vmax=np.median(dat)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

ax2=plt.subplot(1,2,2,sharex=refax,sharey=refax)
plt.imshow(dat2,vmin=0,vmax=np.median(dat2)*3)
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.setp(plt.gca().get_yticklabels(), visible=False)

#step two plot circles where the clicks are in both plots
circles=[refax.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
circles=[ax2.add_patch(mpl.patches.Circle(pt,radius,facecolor='None',
                                         edgecolor='k',linewidth=markerlinewid))
                                         for pt in clicks]
plt.draw()


#Step three run the core of the focus program at those points
peaks0=list(clicks)
im=dat.copy()
wid=7
recenter=True
peaks_dat1=refine_peaks(dat,peaks0,recenter=True)
peaks1=np.array([x for x in zip(peaks_dat1[0],peaks_dat1[1])])

peaks_dat2=refine_peaks(dat2, peaks1,recenter=True)
peaks2=np.array([x for x in zip(peaks_dat2[0],peaks_dat2[1])])

#Now look at home much things shifted
#If the shifts were large then the second round of fits might not have worked
print 'Shift from clicks to 1 {:.2f} pix, clicks to 2 {:.2f} pix, 1 to 2 {:.2f} pix'.format(np.sqrt(((peaks1-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks0)**2).sum(1)).mean(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean())


print 'Y shift {:.2f}pm{:.2f} pix, X {:.2f}pm{:.2f} pix, Net shift {:.2f}pm{:.2f} pix'.format((peaks2-peaks1)[:,1].mean(),(peaks2-peaks1)[:,1].std(),(peaks2-peaks1)[:,0].mean(),(peaks2-peaks1)[:,0].std(),
           np.sqrt(((peaks2-peaks1)**2).sum(1)).mean(),np.sqrt(((peaks2-peaks1)**2).sum(1)).std())
"""
swap filter and back
Y shift -0.24pm1.54 pix, X 0.19pm1.39 pix, Net shift 1.86pm0.97 pix
"""