from __future__ import print_function

try:
    input = raw_input
except NameError:
    pass

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from jbastro.astroLib import gauss2D, gaussfit2D, gauss2dmodel
from skimage.feature import peak_local_max

FITTING_TOL = 1e-2

import warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.utils.exceptions import AstropyUserWarning

# Fit the data using astropy.modeling
fit_p = fitting.LevMarLSQFitter()

import lmfit


def gaussdit2d_lmfit(im, initialp, xwid_bound=(None, None), ywid_bound=(None, None),
                     bound_pos_by=None, maxfev=100, ftol=FITTING_TOL, bounds=True):
    model = lmfit.models.Gaussian2dModel()
    z = im.ravel()
    x, y = np.mgrid[:im.shape[0], :im.shape[1]]
    params = model.guess(z, x.ravel(), y.ravel())
    result = model.fit(z, x=x.ravel(), y=y.ravel(), params=params)
    mim = model.func(x, y, **result.best_values)
    amp, xo, yo, sigma_x, sigma_y = result.best_values
    offset, covar = 0
    return mim, amp, xo, yo, sigma_x, sigma_y, covar, offset


def gaussfit2D_apy(im, initialp, xwid_bound=(None, None), ywid_bound=(None, None),
                   bound_pos_by=None):
    # Fit the data using astropy.modeling

    amp, xo, yo, sigma_x, sigma_y, offset = initialp
    p_init = models.Gaussian2D(x_mean=xo, y_mean=yo, amplitude=amp, x_stddev=sigma_x, y_stddev=sigma_y, theta=-.05)
    p_init.amplitude.bounds = (0, None)
    p_init.x_stddev.bounds = xwid_bound
    p_init.y_stddev.bounds = ywid_bound
    if bound_pos_by is not None:
        p_init.x_mean.bounds = (xo - bound_pos_by[0], xo + bound_pos_by[0])
        p_init.y_mean.bounds = (yo - bound_pos_by[1], yo + bound_pos_by[1])
    # p_init += models.Polynomial2D(0)
    # p_init.c0_0_1 = offset
    x, y = np.mgrid[:im.shape[0], :im.shape[1]]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Model is linear in parameters', category=AstropyUserWarning)
        model = fit_p(p_init, x, y, im, acc=FITTING_TOL)

    return model(x, y), model.parameters


def gaussfit2D(im, initialp, ftol=1e-5, maxfev=5000, maxv=np.inf,
               xwid_bound=(-np.inf, np.inf), ywid_bound=(-np.inf, np.inf),
               bound_pos_by=None, bounds=None):
    x = np.arange(im.shape[1], dtype=np.float)
    y = np.arange(im.shape[0], dtype=np.float)
    x, y = np.meshgrid(x, y)

    amp, xo, yo, sigma_x, sigma_y, covar, offset = initialp
    initialp = list(initialp)

    if bounds is not None:
        if isinstance(bounds, bool) and bounds:
            bounds = [(0, maxv),
                      (-np.inf, np.inf) if bound_pos_by is None else (xo - bound_pos_by[0], xo + bound_pos_by[0]),
                      (-np.inf, np.inf) if bound_pos_by is None else (yo - bound_pos_by[1], yo + bound_pos_by[1]),
                      tuple(np.asarray(xwid_bound)),
                      tuple(np.asarray(ywid_bound)),
                      (-np.inf, np.inf),
                      (im.min(), im.max())]
            bounds = list(zip(*bounds))

    popt, pcov = curve_fit(gauss2D, (x.ravel(), y.ravel()), im.ravel(), p0=initialp,
                           ftol=ftol, maxfev=maxfev, bounds=bounds)
    model = gauss2D((x, y), *popt)

    return model, np.asarray(popt)


def find_peaks(im, min_sep=10, minv=25, maxv=5000):
    """find points to use for psf measurements"""
    points = peak_local_max(im, min_distance=min_sep, threshold_abs=minv, threshold_rel=0)
    return points[im[points[:, 0], points[:, 1]] < maxv]


def plot_stamp(im, xy, wid, **imargs):
    try:
        xw,yw=wid
    except TypeError:
        yw, xw = [wid]*2
    if xw % 2 == 0:
        xw += 1
    if yw % 2 == 0:
        yw += 1
    plt.imshow(im[xy[0] - yw // 2: xy[0] + yw // 2 + 1, xy[1] - xw // 2: xy[1] + xw // 2 + 1], origin='lower',
               interpolation='nearest', **imargs)


def fit_stamp(im, xy, wid=25, sigma=None, fwhm=None, bounds=True, maxv=np.inf, fitter='scipy',
              xwid_bound=(-np.inf, np.inf), ywid_bound=(-np.inf, np.inf), plot=False, **fitargs):
    """fwhm takes precedence
    bounds pos by sigma/3

    uses std along axes for sig est if not provided.
    """
    if fwhm:
        sigma = np.asarray(fwhm) / 2.35482
    try:
        xw, yw = wid
    except TypeError:
        xw, yw = [wid] * 2

    y, x = xy

    sim = im[y - yw // 2:y + yw // 2 + 1, x - xw // 2:x + xw // 2 + 1]
    offset = np.median(sim)
    sig = np.array([sim.mean(1).std(), sim.mean(0).std()])  # for mspec lamp images y/xdispersion, x/spectral
    if sigma is None:
        sigma = sig
    else:
        sigma = np.asarray(sigma)
        # print("User sig {}. Sig Est: {}".format(sigma, sig))

    if sigma[0] < xwid_bound[0] or sigma[0] > xwid_bound[1] or sigma[1] < ywid_bound[0] or sigma[1] > ywid_bound[1]:
        raise RuntimeError('Bad sigma {} / bounds: {},{}'.format(sigma, xwid_bound, ywid_bound))

    p0 = (sim.max() - offset, xw / 2, yw / 2, sigma[0], sigma[1], .05, offset)
    if fitter == 'astropy':
        model, params = gaussfit2D_apy(sim - offset, p0, xwid_bound=xwid_bound, ywid_bound=ywid_bound,
                                       bound_pos_by=sigma / 3)
    elif fitter == 'lmfit':
        model, params = gaussdit2d_lmfit(sim - offset)
    else:
        model, params = gaussfit2D(sim, p0, xwid_bound=xwid_bound, ywid_bound=ywid_bound, maxv=maxv,
                                   bound_pos_by=sigma / 3, bounds=bounds, **fitargs)
    psf = params[3:6]
    psf[0:2] *= 2.35482

    if plot:
        fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
        plt.sca(axes.flat[0])
        plot_stamp(im, xy, wid, vmin=sim.min(), vmax=sim.max())
        plt.title('{},{}'.format(x, y))
        plt.sca(axes.flat[1])
        plt.imshow(model, vmin=sim.min(), vmax=sim.max(), interpolation='nearest', origin='lower')
        plt.title('FWHMx: {:.1f} FWHMy:{:.1f} Covar:{:.1f}'.format(*psf))
        plt.sca(axes.flat[2])
        plt.imshow(sim - model, vmin=sim.min(), vmax=sim.max(), interpolation='nearest', origin='lower')
        plt.title('Difference')
        plt.colorbar()
        plt.show()

    return psf, model


def filter_peaks(im, ret, x_min=2.0, x_max=15.0, y_min=2.0, y_max=15.0,
                 xw=25, yw=25):
    """
    recover data after fitting peaks and make cuts on it

    returns x,y,xw,yw,covar, modelim, modeldat
    """

    modeldat = np.zeros(im.shape)
    modelim = np.zeros(im.shape)
    psfs, modelims = map(np.array, zip(*ret.values()))
    good = (psfs[:, 0] < x_max) & (psfs[:, 0] > x_min) & (psfs[:, 1] < y_max) & (psfs[:, 1] > y_min)
    y, x = np.array(list(ret.keys())).T
    xwid, ywid, covar = psfs[good].T

    for yy, xx, mod in zip(y, x, modelims):
        xsli = slice(xx - xw // 2, xx + xw // 2 + 1)
        ysli = slice(yy - yw // 2, yy + yw // 2 + 1)
        modeldat[ysli, xsli] = im[ysli, xsli]
        modelim[ysli, xsli] = mod

    return x, y, xwid, ywid, covar, modelim, modeldat


def fit_peaks_v2(image, points, window, sigma_guess=None, xwid_bound=(None, None), ywid_bound=(None, None), debug=False):
    ret = {}
    for i, point in enumerate(points):
        if len(points) % (i + 1) == len(points) // 10:
            print('Fitting {} of {}'.format(i, len(points)))
        try:
            if (point-window//2 <0).any() or (point+window//2 > np.asarray(image.shape)).any():
                continue

            if debug and i>0:
                plt.clf()
            psf, model = fit_stamp(image, point, window, sigma=sigma_guess, xwid_bound=xwid_bound,ywid_bound=ywid_bound, plot=debug)
            ret[tuple(point)] = psf, model

            if debug:
                val = input('Debug peak fitter? (n/a/db)>')
                if val == 'n':
                    debug = False
                elif val == 'a':
                    return ret
        except RuntimeError as e:
            print('Point {} failed. {}'.format(i, e))

    return ret

def fit_peaks(image, points, wx, wy, xw=25, yw=25, covar=.05, debug=False,
              maxfev=5000, xwid_bound=(None, None), ywid_bound=(None, None)):
    """fit the psf at each point"""
    ret = {}
    xw = int(round(xw))
    yw = int(round(yw))
    sigu = wy / 2.35482, wx / 2.35482
    for i, (y, x) in enumerate(points):

        try:
            if (y - yw / 2 < 0 or y + yw / 2 >= image.shape[0] or
                    x - xw / 2 < 0 or x + xw / 2 >= image.shape[1]):
                continue

            sim = image[y - yw // 2:y + yw // 2 + 1, x - xw // 2:x + xw // 2 + 1]
            offset = sim.min()
            sig = sim.mean(0).std(), sim.mean(1).std()
            if sig[0] > sigu[0] or sig[1] > sigu[1]:
                continue

            covar = sig[1] / sig[0]
            if (sim.mean(0) * np.arange(sim.shape[1])).sum() / sim.mean(0).sum() > sim.shape[1] / 2:
                covar *= -1

            # print("User sig. {}. Sig. est.: {}".format(sig, sige))
            fit = gaussdit2d_lmfit(sim, (
                image[y, x] - offset, sim.shape[0] / 2, sim.shape[1] / 2, sig[0], sig[1], covar, offset),
                                   maxfev=maxfev, ftol=FITTING_TOL,
                                   bounds=True, xwid_bound=xwid_bound, ywid_bound=xwid_bound, bound_pos_by=(2, 2))
            psf = np.array(fit[1][3:6])
            psf[0:2] *= 2.35482
            model = fit[0]
            ret[(y, x)] = psf, model

            if debug:
                plt.figure(10)
                plt.subplot(1, 2, 1)
                plt.imshow(sim, vmin=sim.min(), vmax=image[y, x], interpolation='nearest')
                plt.title('{},{}'.format(x, y))
                plt.subplot(1, 2, 2)
                plt.imshow(model, vmin=sim.min(), vmax=image[y, x], interpolation='nearest')
                plt.title('FWHMx: {:.1f} FWHMy:{:.1f} Covar:{:.1f}'.format(*psf))
                plt.show()
                val = input('Debug peak fitter? (n/a/db)>')
                if val == 'n':
                    debug = False
                elif val == 'a':
                    return ret
                elif val == 'db':
                    pass

        except RuntimeError as e:
            print('Point {} failed. {}'.format(i, e))

    return ret


def fit_peaks_apy(image, points, wx, wy, xw=25, yw=25, covar=.05, debug=False, maxfev=5000,
                  xwid_bound=(None, None), ywid_bound=(None, None)):
    """fit the psf at each point"""
    ret = {}
    xw = int(round(xw))
    yw = int(round(yw))
    for i, (y, x) in enumerate(points):
        if len(points) % (i + 1) == len(points) // 10:
            print('Fitting {} of {}'.format(i, len(points)))

        try:
            if y < yw // 2 or y + yw // 2 >= image.shape[0] or x < xw // 2 or x + xw // 2 >= image.shape[1]:
                continue

            sim = image[y - yw // 2:y + yw // 2 + 1, x - xw // 2:x + xw // 2 + 1]
            offset = sim.min()
            sig = sim.mean(0).std(), sim.mean(1).std()
            # sig = wy / 2.35482, wx / 2.35482
            # print("User sig {}. Sig Est: {}".format(sig, sige))
            amp = image[y, x] - offset
            model, params = gaussfit2D_apy(sim - offset, (amp, yw / 2, xw / 2, sig[0], sig[1], offset),
                                           xwid_bound=xwid_bound, ywid_bound=ywid_bound,
                                           bound_pos_by=(2, 2))
            psf = params[3:6]
            psf[0:2] *= 2.35482
            ret[(y, x)] = psf, model

            if debug:
                plt.figure(10)
                plt.subplot(1, 2, 1)
                plt.imshow(sim, vmin=sim.min(), vmax=image[y, x], interpolation='nearest')
                plt.title('{},{}'.format(x, y))
                plt.subplot(1, 2, 2)
                plt.imshow(model, vmin=sim.min(), vmax=image[y, x], interpolation='nearest')
                plt.title('FWHMx: {:.1f} FWHMy:{:.1f} Covar:{:.1f}'.format(*psf))
                plt.show(block=True)
                val = input('Continue? (any)/s(top prompt)/(a)bort/d(ebug))>')
                if val == 's':
                    debug = False
                elif val == 'a':
                    return ret
                elif val == 'd':
                    pass

        except RuntimeError as e:
            print('Point {} failed. {}'.format(i, e))

    return ret
