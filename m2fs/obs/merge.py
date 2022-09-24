#!/usr/bin/env python
import numpy as np
from astropy.io import fits
import m2fs.ccd
import re
import os
import jbastro
import jbastro.astroLib as astroLib
import m2fs.obs
from multiprocessing import Pool

from typing import List, Tuple, Dict


def _proc_quad(qdata, header, cosmic_settings):
    """cosmic_settings = cosmic_kwargs + 'iter' or just 'iter':0"""

    cosmic_iter = cosmic_settings.pop('iter')

    biassec = header['BIASSEC']
    bias = [int(s) for s in re.findall(r'\d+', biassec)]

    trimsec = header['TRIMSEC']
    crop = [int(s) for s in re.findall(r'\d+', trimsec)]

    qdata = qdata.astype(float)
    # Compute & subtract median for overscan region rowwise
    biaslevels = np.median(qdata[crop[2] - 1:, bias[0] - 1:bias[1]].astype(float), axis=1)
    qdata -= biaslevels[:, np.newaxis]

    # Compute & subtract median bias row
    qdata -= np.median(qdata[bias[2]:, :], axis=0)

    # Crop the image
    qdata = qdata[crop[2] - 1:crop[3], crop[0] - 1:crop[1]]

    # Cosmic ray rejection
    if cosmic_iter > 0:
        qmask = astroLib.crreject(qdata, dialate=True, gain=header['EGAIN'], readnoise=header['ENOISE'],
                                  satlevel=.95 * m2fs.ccd.satlevel, **cosmic_settings)

        # c = cosmicsimage(qdata, gain=header['EGAIN'],
        #                  readnoise=header['ENOISE'],
        #                  satlevel=.95 * m2fs.ccd.satlevel,
        #                  **cosmic_settings)
        # c.run(maxiter=cosmic_iter)
        # qmask = ndimage.morphology.binary_dilation(
        #     c.mask.astype(np.bool), structure=np.ones((3, 3)),
        #     iterations=1, mask=None, output=None,
        #     border_value=0, origin=0, brute_force=False)

    else:
        qmask = np.zeros_like(qdata, dtype=np.bool)

    # Convert to e-
    qdata *= header['EGAIN']

    # Variance image
    qvari = qdata + header['ENOISE'] ** 2

    return qdata, qvari, qmask


def mergequad(frameno, side=None, do_cosmic=False, file=False, odir='', idir='',
              repair=False, dogzip=False, overwrite=False):
    """Give a seqno or a path to a quad if file set
    do_cosmic=bool or dict like
    sigma for init clip, fraction for neighbors, how much above background*
    {'sigclip': 7.0, 'sigfrac': 0.428, 'objlim': 1.4, 'iter':10}
    """
    if side is None and not file:
        try:
            mergequad(frameno, side='r', do_cosmic=do_cosmic, idir=idir, odir=odir,
                      repair=repair, dogzip=dogzip, overwrite=overwrite)
        except IOError:
            print('Need all quadrants for r{}'.format(frameno))
        try:
            mergequad(frameno, side='b', do_cosmic=do_cosmic, idir=idir, odir=odir,
                      repair=repair, dogzip=dogzip, overwrite=overwrite)
        except IOError:
            print('Need all quadrants for b{}'.format(frameno))
        return

    if file:
        fname = frameno
        frameno = int(m2fs.obs.info(fname, no_load=True).seqno)
        basename = os.path.basename(fname)
        side = basename[0]
        idir = os.path.dirname(fname) + os.path.sep if not idir else idir
        _, _, extension = basename.partition('.')
    else:
        idir = idir
        extension = 'fits'

    basename = side + '{frameno:04}'
    f = idir + basename + 'c{quad}.' + extension

    if ((os.path.exists(odir + basename.format(frameno=frameno) + '.fits') or
         os.path.exists(odir + basename.format(frameno=frameno) + '.fits.gz'))
            and not overwrite):
        print("Skipping " + odir + basename.format(frameno=frameno) + ".fits. Already done.")
        return
    else:
        print("Merging " + odir + basename.format(frameno=frameno) + '.fits')

    gzip = '.gz' if dogzip else ''
    ofile = odir + basename.format(frameno=frameno) + '.fits' + gzip
    if not os.path.exists(os.path.dirname(odir)):
        os.makedirs(os.path.dirname(odir))

    _mergequad((f.format(frameno=frameno, quad=1),
                f.format(frameno=frameno, quad=2),
                f.format(frameno=frameno, quad=3),
                f.format(frameno=frameno, quad=4)),
               ofile, do_cosmic=do_cosmic,
               pool=4, repair=repair, overwrite=overwrite)


def _mergequad(quadrants: (Tuple, List), ofile, do_cosmic: (bool, Dict)=False, pool: (Pool, int)=4,
               repair=False, overwrite=False):
    """quadrants: (Tuple, List), ofile, do_cosmic: bool, dict}=False, pool: {Pool, int}=4,
               repair=False, overwrite=False"""
    if not overwrite and os.path.exists(ofile):
        return
    # Load the data
    try:
        quadrant1 = fits.open(quadrants[0])
        quadrant2 = fits.open(quadrants[1])
        quadrant3 = fits.open(quadrants[2])
        quadrant4 = fits.open(quadrants[3])
    except IOError as e:
        raise e

    # Cosmic ray removal configuration
    if isinstance(do_cosmic, bool):
        cosmic_settings = {'sigclip': 7.0, 'sigfrac': 0.428, 'objlim': 1.4, 'iter': 10} if do_cosmic else {'iter': 0}
    else:
        assert isinstance(do_cosmic, dict)
        cosmic_settings = do_cosmic.copy()

    # Grab the bias and crop region from the first quadrant
    trimsec = quadrant1[0].header['TRIMSEC']
    crop = [int(s) for s in re.findall(r'\d+', trimsec)]

    # Create output arrays
    out = np.zeros((crop[3] * 2, crop[1] * 2), dtype=np.float32)
    vout = np.zeros((crop[3] * 2, crop[1] * 2), dtype=np.float32)
    mask = np.zeros((crop[3] * 2, crop[1] * 2), dtype=np.uint8)
    headerout = quadrant1[0].header.copy()

    # Define where the quadrants show go
    quadLoc = [(0, crop[3], 0, crop[1]),
               (0, crop[3], crop[1], 2 * crop[1]),
               (crop[3], 2 * crop[3], crop[1], 2 * crop[1]),
               (crop[3], 2 * crop[3], 0, crop[1])]

    # Create a list of the quadrant's data to loop through
    quadrantData = [quadrant1[0], quadrant2[0], quadrant3[0], quadrant4[0]]

    # Process them in multiple processes
    args = [(q.data.copy(), q.header.copy(), cosmic_settings.copy())
            for q in quadrantData]

    _pool = Pool(processes=pool) if isinstance(pool, int) else pool

    res = [_pool.apply_async(_proc_quad, arg) for arg in args]
    if isinstance(pool, int):
        _pool.close()
        _pool.join()

    for i, r in enumerate(res):
        qdata, qvar, qmask = r.get()
        # Position the quadrants
        if i == 0:
            out[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = qdata
            vout[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = qvar
            mask[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = qmask
        if i == 1:
            out[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.fliplr(qdata)
            vout[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.fliplr(qvar)
            mask[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.fliplr(qmask)
        if i == 2:
            out[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.rot90(qdata, 2)
            vout[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.rot90(qvar, 2)
            mask[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.rot90(qmask, 2)
        if i == 3:
            out[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.rot90(np.fliplr(qdata), 2)
            vout[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.rot90(np.fliplr(qvar), 2)
            mask[quadLoc[i][0]:quadLoc[i][1], quadLoc[i][2]:quadLoc[i][3]] = np.rot90(np.fliplr(qmask), 2)

    if repair:
        jbastro.median_patch(out, mask)

    # Flip so it is in agreement with Mario's process and convert to 32bit float
    out = np.flipud(out).astype(np.float32)
    vout = np.flipud(vout).astype(np.float32)
    mask = np.flipud(mask)

    # make CRs, NaNs, & Infs count as nothing
    bad = ((mask.astype(bool)) | (~np.isfinite(out)) | (~np.isfinite(vout)))
    out[bad] = 0
    patch_val = min(max(vout[~bad].max() * 1.05, 1e9), 1e35)
    vout[bad] = patch_val

    # Write out the merged file
    headerout.pop('TRIMSEC')
    headerout.pop('BIASSEC')
    headerout.pop('DATASEC')
    fname = os.path.splitext(os.path.splitext(os.path.basename(ofile))[0])[0]
    headerout['FILENAME'] = fname
    headerout['BUNIT'] = 'E-/PIXEL'
    headerout['EGAIN'] = 1.0
    headerout['ENOISE'] = 2.5  # average for the 4 amps in slow is 2.5 e
    headerout['PATCHVAL'] = patch_val
    if cosmic_settings['iter'] == 0:
        headerout['LACR'] = 'False'
    else:
        headerout['LACR'] = 'True'
        headerout['LACRIT'] = cosmic_settings['iter']
        headerout['LACRSC'] = cosmic_settings['sigclip']
        headerout['LACRSF'] = cosmic_settings['sigfrac']
        headerout['LACROL'] = cosmic_settings['objlim']

    hdul = fits.HDUList(fits.PrimaryHDU(header=headerout))
    hdul.append(fits.ImageHDU(out, name='science', header=headerout))
    hdul.append(fits.ImageHDU(vout, name='variance', header=headerout))

    if do_cosmic:
        hdul.append(fits.ImageHDU(mask, name='crmask'))

    hdul.writeto(ofile, overwrite=True)

    return hdul
