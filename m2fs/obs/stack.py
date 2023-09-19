import numpy as np
from astropy.io import fits
import scipy.ndimage as ndimage
from jbastro.astroLib import crreject
import os
from astropy.time import Time, TimeDelta
from logging import getLogger
from m2fs.seqno import MSpecFits

def compute_patch_val(variances):
    patch_val = max(1.05 * variances.max(), 1e8)
    if patch_val >= 1e38:
        print('WARNING: Real variances {}'.format(patch_val) +
              'exceed maximum bad variance patch value (1e38) for float32!')
    return min(patch_val, 1e38)


def mask_problem_spots(var, header, mask_val=1e35):
    if 'r' in header['FILENAME']:
        if header['BINNING'] == '1x1':
            var[:2056, 821] = mask_val
        if header['BINNING'] == '2x2':
            var[:2056, 410] = mask_val
    if 'b' in header['FILENAME'] and header['BINNING'] == '1x1':
        var[328, :2048] = mask_val
        var[1600:2056, 2518] = mask_val


def stack_images(files, outfile, do_cosmic=False, overwrite=False, flux_weighted=True, **crparams):
    """
    List of files to stack, output file sans extension.


    output hdu
    cleaned image
    variance image
    summed cr mask
    bad pixel mask

    """
    if os.path.exists(outfile) and not overwrite:
        return

    nfile = len(files)
    # load first to get info
    msf = MSpecFits(files[0])
    if msf.bias or msf.dark:
        flux_weighted=False
        
    with fits.open(files[0]) as im:
        header = im[0].header
        # 15 frames with float64 would be 3.7GB ram
        imcube = np.zeros(im[1].data.shape + (nfile,), dtype=np.float32)
        if len(im) > 3:
            masked = True
            mask = np.zeros_like(imcube, dtype=bool)
        else:
            masked = False


    varcube = np.zeros_like(imcube)

    # load data and variance frames
    members = []
    etimes = []
    midpoint_weights = []
    midpoints = []
    airmasses = []
    first_uttime = None
    last_utend = None
    try:
        for i, f in enumerate(files):
            with fits.open(f) as im:
                airmasses.append(im[0].header['AIRMASS'])
                etimes.append(float(im[0].header['EXPTIME']))
                members.append(im[0].header['FILENAME'])

                imcube[:, :, i] = im[1].data
                varcube[:, :, i] = im[2].data

                if masked:
                    msk = im[3].data
                    if do_cosmic:
                        foo = crreject(im[2].data)
                        msk += ndimage.binary_dilation(
                            foo.astype(np.bool), structure=np.ones((3, 3)),
                            iterations=1, mask=None, output=None,
                            border_value=0, origin=0, brute_force=False)

                    mask[:, :, i] = msk

                # compute midpoint
                start = Time(im[0].header['UT-DATE'] + ' ' + im[0].header['UT-TIME'],
                             format='iso', scale='utc')
                end = Time(im[0].header['UT-DATE'] + ' ' + im[0].header['UT-END'],
                           format='iso', scale='utc')
                midpoints.append(start + .5 * (end - start))

                # keep track of stacked beginning and end
                if first_uttime is None or first_uttime > start:
                    first_uttime = start
                if last_utend is None or last_utend < end:
                    last_utend = end

    # import ipdb;ipdb.set_trace() #midpoints untested

    except ValueError as e:
        print("ValueError while merging:", files, f)
        return

    if masked:
        use = ~mask.sum(2).astype(bool)
        midpoint_weights = [im.T[use].sum() for im in imcube.T]
    else:
        midpoint_weights = [im.sum() for im in imcube.T]

    mean_airmass = np.mean(airmasses)
    min_midpoint = min(midpoints)
    midpoint = min_midpoint + TimeDelta(np.average([(m - min_midpoint).sec
                                                    for m in midpoints],
                                                   weights=midpoint_weights),
                                        format='sec')

    exptime = sum(etimes)

    if not masked and do_cosmic:
        # Time to reject cosmic rays
        clipped = crreject(imcube, **crparams)
        mask = clipped.mask
    elif not masked:
        mask = np.zeros_like(imcube, dtype=bool)

    # this code issues warnings like crazy, but the math is ok because of a
    # final pass collecitng all the garbage
    import warnings
    with warnings.catch_warnings():
        # Create masked arrays
        imcube_masked = np.ma.array(imcube, mask=mask)
        varcube_masked = np.ma.array(varcube, mask=mask)

        duration_corr = max(etimes) / np.array(etimes)

        throughput_corr = np.array(midpoint_weights) / np.array(etimes)
        throughput_corr[:] = throughput_corr.max() / throughput_corr
        if not flux_weighted:
            getLogger(__name__).info(f'Straigt stacking')
            throughput_corr[:]=1

        ####### Try 3
        # the correction factor as defined has a pathological edge case:
        # if high throughput fibers or brighter targets are disproportionatly
        # affected by clouds/field rotation, guiding errors, or the like then it
        # would give unfair weight to that frame.

        corr_fac = throughput_corr * duration_corr
        patch_cube = (imcube_masked * corr_fac).mean(axis=2)
        patch_cube = patch_cube[:, :, np.newaxis] / corr_fac
        imcube[mask] = patch_cube[mask]
        im = imcube.sum(axis=2)

        # fill value doesn't matter as mask check below will flag and cause patch
        var = varcube_masked.sum(axis=2).filled(0)
        var *= (im / imcube_masked.sum(axis=2)) ** 2

        mask_problem_spots(var, header, mask_val=np.nan)

        im = im.astype(np.float32)
        var = var.astype(np.float32)

        bad = (mask.all(axis=2) | (~np.isfinite(im)) | (~np.isfinite(var)))

        patch_val = compute_patch_val(var[~bad])
        im[bad] = 0.0
        var[bad] = patch_val

    if (~np.isfinite(var)).any():
        import ipdb
        ipdb.set_trace()
    if (~np.isfinite(im)).any():
        import ipdb
        ipdb.set_trace()

    # Update the header
    header['FILENAME'] = os.path.basename(outfile)
    header['EXPTIME'] = exptime
    header['COMMENT'] = ','.join(members)
    header['UT-MID'] = str(midpoint)
    header['UT-TIME'] = str(first_uttime)
    header['UT-END'] = str(last_utend)
    header['AIRMASS'] = mean_airmass
    header['PATCHVAL'] = patch_val

    print('Setting read noise and to approximate value. Use variance frame for exact value')
    header['EGAIN'] = 1.0
    header['ENOISE'] = np.sqrt(nfile) * 2.5  # average for the 4 amps in slow is 2.5 e

    # Create a primary extension with header only
    #   (this seems to be a fits convention)
    hdul = fits.HDUList(fits.PrimaryHDU(header=header))

    # Append image extensions
    hdul.append(fits.ImageHDU(im, name='science', header=header))
    hdul.append(fits.ImageHDU(var, name='variance', header=header))
    hdul.append(fits.ImageHDU(mask.sum(axis=2).astype(np.uint8), name='crmask'))
    hdul.append(fits.ImageHDU(np.zeros_like(im, dtype=np.uint8), name='bpmask'))
    hdul.writeto(outfile, overwrite=overwrite)


def patch_single(file, outfile, do_cosmic=False, overwrite=False, **extra):
    if os.path.exists(outfile) and not overwrite:
        return
    if do_cosmic:
        raise NotImplementedError('No support for single frame CR rejection at present.')
    hdul = fits.open(file)
    mv = compute_patch_val(hdul['VARIANCE'].data)
    mask_problem_spots(hdul['VARIANCE'].data, hdul['SCIENCE'].header, mask_val=mv)
    hdul.writeto(outfile, overwrite=overwrite)
