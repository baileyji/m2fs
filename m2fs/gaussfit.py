from astropy.modeling import models, fitting
from jbastro.astroLib import gauss2D, gaussfit2D
from scipy.optimize import curve_fit
import lmfit
import warnings
import numpy as np
from astropy.utils.exceptions import AstropyUserWarning


FITTING_TOL = 1e-2


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
    fit_p = fitting.LevMarLSQFitter()
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
