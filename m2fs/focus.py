from __future__ import print_function

try:
    input = raw_input
except NameError:
    pass

import warnings
import numpy as np
import matplotlib.pyplot as plt
from m2fs.gaussfit import gaussfit2D, gaussfit2D_apy, gaussdit2d_lmfit

FITTING_TOL = 1e-2
R_FOCUS_PARAM = [-2.1543, 261.5]  # ifum run 1 zp
B_FOCUS_PARAM = [1.7494, 178.0]

RFOCUS = np.poly1d(R_FOCUS_PARAM)
BFOCUS = np.poly1d(B_FOCUS_PARAM)


def plot_stamp(im, xy, wid, **imargs):
    try:
        xw, yw = wid
    except TypeError:
        yw, xw = [wid] * 2
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
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='overflow encountered in exp')
            warnings.filterwarnings('ignore', message='overflow encountered in multiply')
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


def filter_peaks(im, ret, xwid_bound, ywid_bound, xw=25, yw=25):
    """
    recover data after fitting peaks and make cuts on it

    returns x,y,xw,yw,covar, modelim, modeldat
    """
    x_min, x_max = xwid_bound
    y_min, y_max = ywid_bound

    modeldat = np.zeros(im.shape)
    modelim = np.zeros(im.shape)
    psfs, modelims = map(np.array, zip(*ret.values()))
    good = (psfs[:, 0] < x_max) & (psfs[:, 0] > x_min) & (psfs[:, 1] < y_max) & (psfs[:, 1] > y_min)
    print('Excluding {} of {} peaks due to FWHM bounds'.format(good.size - good.sum(), good.size))
    y, x = np.array(list(ret.keys())).T
    xwid, ywid, covar = psfs[good].T

    for yy, xx, mod in zip(y, x, modelims):
        xsli = slice(xx - xw // 2, xx + xw // 2 + 1)
        ysli = slice(yy - yw // 2, yy + yw // 2 + 1)
        modeldat[ysli, xsli] = im[ysli, xsli]
        modelim[ysli, xsli] = mod

    return x, y, xwid, ywid, covar, modelim, modeldat


def fit_peaks(image, points, window, sigma_guess=None, xwid_bound=(None, None), ywid_bound=(None, None), debug=False):
    ret = {}
    for i, point in enumerate(points):
        if len(points) % (i + 1) == len(points) // 10:
            print('Fitting {} of {}'.format(i, len(points)))
        try:
            if (point - window // 2 < 0).any() or (point + window // 2 > np.asarray(image.shape)).any():
                continue

            if debug and i > 0:
                plt.clf()
            psf, model = fit_stamp(image, point, window, sigma=sigma_guess, xwid_bound=xwid_bound,
                                   ywid_bound=ywid_bound, plot=debug)
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
