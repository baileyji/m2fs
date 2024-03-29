#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from scipy.ndimage import median_filter
from matplotlib.pyplot import *

import scipy.signal as sp_sig
from scipy.interpolate import UnivariateSpline
from scipy.ndimage.filters import gaussian_filter as gsmooth2d
import numpy as np
import matplotlib.pyplot as plt
from jbastro.astroLib import gauss2D, gaussfit2D


def make_glow_model(im_in, bin_x=1, bin_y=1):
    """ Model the amplifier glow """
    im = im_in.copy()
    im[0] = im[2]
    im[1] = im[2]
    im[-1] = im[-2]

    # glow image
    glow = np.zeros_like(im)

    # meshgrid
    x, y = np.meshgrid(np.arange(im.shape[1]), np.arange(im.shape[0]))

    def model_corner(im, x0, y0, xw, yw, iparams, std_clip=0):
        """ std_clip is the y height of the small corner to use to exclude
        spectra in the large corner,

        (iparams=(glow amp, x center, y center, xwid, ywid, xy amount)

        positions and initial params adjusted automatically for binning
        pass coordinates in 4k positions
        """
        x0 //= bin_x
        y0 //= bin_y
        xw //= bin_x
        yw //= bin_y
        iparams = list(iparams)
        iparams[1] //= bin_x
        iparams[2] //= bin_y
        iparams[3] //= bin_x
        iparams[4] //= bin_y

        corner = im[y0:y0 + yw, x0:x0 + xw].copy()
        if std_clip:
            small_corner = im[y0:y0 + std_clip, x0:x0 + xw].copy()
            patch_locs = corner > 2 * small_corner.std()
            patch_locs[:y0 + std_clip, :] = False
            corner[patch_locs] = np.median(small_corner)
        cim, param = gaussfit2D(corner, iparams)
        param = list(param)
        param[-1] = 0
        param[1] += x0
        param[2] += y0
        return gauss2D((x, y), *param)

    # Lower R
    try:
        tmp = model_corner(im, 3996, 2, 100, 100, (150, 58, -7, 30.0, 20.0, 0, 0))
        if tmp.min() < 0:
            raise RuntimeError('Glow model has negative values')
        else:
            glow += tmp
    except RuntimeError as e:
        print('Lower R glow model failed: {}'.format(str(e)))

    # Lower L
    try:
        tmp = model_corner(im, 0, 2, 100, 100, (150, 40, 0, 30.0, 20.0, 0, 0), std_clip=50)
        if tmp.min() < 0:
            raise RuntimeError('Glow model has negative values')
        else:
            glow += tmp
    except RuntimeError as e:
        print('Lower L glow model failed: {}'.format(str(e)))

    # Upper L
    try:
        tmp = model_corner(im, 0, 4012, 100, 100, (150, 40, 100, 30.0, 20.0, 0, 0))
        if tmp.min() < 0:
            raise RuntimeError('Glow model has negative values')
        else:
            glow += tmp
    except RuntimeError as e:
        print('Upper L glow model failed: {}'.format(str(e)))

    # Upper R
    try:
        tmp = model_corner(im, 3996, 4000, 100, 100, (150, 58, 100, 30.0, 20.0, 0, 0))
        if tmp.min() < 0:
            raise RuntimeError('Glow model has negative values')
        else:
            glow += tmp
    except RuntimeError as e:
        print('Upper R glow model failed: {}'.format(str(e)))

    return glow


def mkscatter(im_in, var_im, plot=False, scatter_thresh=.7, debug=False, header=None,
              do_glow=True, offset=0, glowOnly=False):
    """Make a scattered light map"""

    do_glow |= glowOnly
    im_in = im_in.copy()

    if debug:
        import ipdb
        ipdb.set_trace()

    if header is not None:
        bin_x, bin_y = int(header['BINNING'][0]), int(header['BINNING'][-1])
        if bin_x != 1 and bin_y != 1:
            print('WARNING: Non 1x1 image, check results carefully.')
    else:
        print('WARNING: No header assuming 1x1 image, check results carefully.')
        bin_x, bin_y = 1, 1

    # model amp glow
    glow = make_glow_model(im_in, bin_x=bin_x, bin_y=bin_y) if do_glow else np.zeros_like(im_in)

    # first two rows and last row are abberations
    im_in[0] = im_in[2]
    im_in[1] = im_in[2]
    im_in[-1] = im_in[-2]

    # regions to estimate the scattered light
    scatter_regions = [(0, 50), (480, 575), (1010, 1100), (1525, 1600), (2035, 2120),
                       (2540, 2620), (3030, 3110), (3525, 3600), (3995, 4109)]

    scatter_regions = np.array(scatter_regions) // bin_y

    scatter_regions[0, 1] += offset
    scatter_regions[-1, 0] += offset
    scatter_regions[1:-1, :] += offset

    # Select all pixels in the scatter regions
    msk = np.ones_like(im_in, dtype=np.bool)
    for r0, r1 in scatter_regions:
        msk[r0:r1 + 1, :] = False
    s_pix = (im_in - glow)[~msk]

    # Compute a ~robust stddev for the scattered light
    s_std = s_pix[abs(s_pix - s_pix.mean()) < 2 * s_pix.std()].std()

    # from astropy.modeling import models, fitting
    # p_init = models.Polynomial2D(degree=5)
    # f = fitting.LinearLSQFitter()
    # y, x = np.mgrid[:im_masked.shape[0], :im_masked.shape[1]]
    # p = f(p_init, x, y, im_masked)

    # median filter to figure out patch values for bad pixels
    # (e.g. excessive variance)
    patch_im = sp_sig.medfilt2d(im_in, 5)
    patch_spot = (var_im > 1e6)
    im_in[patch_spot] = patch_im[patch_spot]

    # Create a masked image that includes all pixels within X stddevs
    # of the mean scattered light level, excluding the glow
    im_mask = abs(im_in - glow - s_pix.mean()) > scatter_thresh * s_std
    im_masked = np.ma.array(im_in - glow, copy=True, mask=im_mask)

    s_model = im_in.copy()

    nrow, ncol = im_in.shape
    rows = np.arange(nrow)
    cols = np.arange(ncol)

    # Model each column with a polynomial fit to the masked image
    for c in cols:
        s_model[:, c] = np.poly1d(np.ma.polyfit(rows, im_masked[:, c], 8))(rows)

    s_modelraw = s_model.copy()
    s_im = glow
    err = 0.0
    if not glowOnly:
        # Smooth the model & add the glow to make the scattered light image
        s_model = gsmooth2d(s_model, (32, 64))  # 32 ~1.5 fiber separations

        # Model each row with a polynomial fit to the modeled columns
        #col_weights=1/var_im.sum(0)
        #for r in rows:
        #   s_model[r,:]=np.poly1d(np.ma.polyfit(cols, s_model[r,:],8, w=col_weights))(cols)

        s_im = s_model + glow
        err = (im_masked - s_im).std()

    # Debugging plots
    if plot or debug:
        plt.figure(1)
        cmap = plt.cm.jet
        cmap.set_bad('black')
        im_vmin = im_masked.min()
        im_vmax = im_masked.mean() + 2 * im_masked.std()
        plt.imshow(im_masked, vmin=im_vmin, vmax=im_vmax, cmap=cmap, origin='lower')
        plt.title('scatter data')
        plt.colorbar()

        plt.figure(2)
        plt.imshow(s_modelraw, vmin=0, vmax=s_modelraw.max(), origin='lower')
        plt.title('raw polys')
        plt.colorbar()

        plt.figure(3)
        vmax_scatter = s_im.mean() + 2 * s_im.std()
        plt.imshow(s_im, vmin=0, vmax=vmax_scatter, origin='lower')
        plt.title('model')
        plt.colorbar()

        plt.figure(4)
        clean = im_in - s_im
        # im_vmin=0
        # im_vmax=vmax_scatter
        plt.imshow(clean, vmin=im_vmin, vmax=im_vmax, origin='lower')
        plt.title('cleaned')
        plt.colorbar()

        for sr in scatter_regions:
            for x in sr:
                axhline(x, color='black')

        plt.show(block=False)

        if debug:
            import ipdb
            ipdb.set_trace()

        return s_im, err, (s_model, im_masked, glow)

    return s_im, err
