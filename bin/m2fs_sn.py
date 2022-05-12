#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import argparse
import ipdb
import os.path
from m2fs.focus import find_peaks


def parse_cl():
    parser = argparse.ArgumentParser(description='S/N Check',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', metavar='FILE', type=str,
                        help='File to process')

    # Finding settings
    parser.add_argument('-d', dest='sep', default=10,
                        action='store', required=False, type=int,
                        help='Minimum delta between order peaks')
    parser.add_argument('-m', dest='min', default=12,
                        action='store', required=False, type=int,
                        help='Minimum peak intensity')

    # Measurement selection settings
    parser.add_argument('-s', dest='hwid', default=8.0,
                        action='store', required=False, type=int,
                        help='Width to either side of smash')

    parser.add_argument('-c', dest='min_sn', default=0,
                        action='store', required=False, type=int,
                        help='Minimum S/N cut')

    parser.add_argument('--tpdat', dest='tpdat', default=False,
                        action='store_true', required=False,
                        help='give numbers needed for throughput checking')

    return parser.parse_args()


# def find_peaks(im, min_sep=10, minv=25, maxv=5000):
#     """find points to use for psf measurements"""
#     points=peak_local_max(im, min_distance=min_sep,
#                           threshold_abs=minv,
#                           threshold_rel=0, indices=False)
#     points[((im>maxv) & points)]=False
#
#     return np.where(points)[0]


def compute_sn(file, args):
    """ a merged frame
    args is an object with attrs
        sep - a minimum allowable seperation
        min - a minimum peak value to consider
        hwid - half window width to use when smashing the gaussian or data
        min_sn - a minimun sn value to consider
        tpdat - wether to print rotator info
        plot - do some plotting
        debug - use debugging mode

    """

    # Settings
    min_sep = args.sep
    minv = args.min
    smashw = args.hwid
    maxv = 1e99

    # fetch image
    hdu = fits.open(file)
    header = hdu[0].header
    side = os.path.basename(file)[0].lower()
    if side not in 'rb':
        side = '?'

    binning = int(header['BINNING'][2])
    min_sep //= binning
    smashw //= binning

    fim = hdu['SCIENCE'].data
    vfim = hdu['VARIANCE'].data
    etime = hdu[0].header['EXPTIME']
    rotv = hdu[0].header['ROTANGLE']

    ncol = fim.shape[1]
    c0 = int(round(ncol / 2 - 2 * .02 * ncol))
    c1 = int(round(ncol / 2 - .02 * ncol))

    im = np.median(fim[:, c0:c1], 1)
    vim = np.median(vfim[:, c0:c1], 1)

    # Find the peaks
    peaks = find_peaks(im, min_sep=min_sep, minv=minv, maxv=maxv)

    print('Found {} peaks between {} and {} separated by at least {}.'.format(
        len(peaks), minv, maxv, min_sep))

    # Plot image and detected peaks
    plt.figure(figsize=(18, 9))
    plt.plot(im)
    for x in peaks: plt.plot(x, im[x], 'w*')

    # Smash peaks
    sns = []
    peakdat = []
    for i, x in enumerate(peaks):

        cent = x
        sl = slice(int(round(cent - smashw)), int(round(cent + 1 + smashw)))
        sn = im[sl].sum() / np.sqrt(vim[sl].sum())

        # Need to add some clipping here for major outliers in the variance

        if sn > args.min_sn:
            peakdat.append(x)
            sns.append(sn)
            plt.text(x, 1.01 * im[x], '{:.0f}'.format(sn))

            plt.gca().add_patch(matplotlib.patches.Rectangle((x - smashw, 0),
                                                             2 * smashw + 1, im[x],
                                                             alpha=.2))

    plt.title('{}\nMedian of {:.0f} pixels at col {:.0f}'.format(os.path.basename(file),
                                                                 c1 - c0, .5 * (c0 + c1)) +
              '  Median S/N {:.0f}'.format(np.median(sns)))
    plt.xlabel('Row (R1-01/B8-01 to left)')
    plt.ylabel('electrons')
    plt.axhline(minv)
    plt.xlim(0, fim.shape[0])
    plt.show()

    for i, sn in enumerate(sns):
        print(i, sn)

    print('Mean/Median/Std S/N: {:.0f}/{:.0f}/{:.0f}'.format(np.mean(sns),
                                                             np.median(sns),
                                                             np.std(sns)))

    if args.tpdat:
        print('Rotator value: {:.2f}'.format(rotv))
        print('{:4}   {}'.format('Row', 'e-/s'))
        for x, sn in zip(peakdat, sns):
            print('{:<4} {:.3f}'.format(x, sn ** 2 / etime))

    return zip(peakdat, sns)


if __name__ == '__main__':
    import matplotlib as mpl

    mpl.rcParams['font.size'] = 16
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['ytick.major.size'] = 6
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['xtick.major.width'] = 2

    args = parse_cl()

    compute_sn(args.file, args)

    print('Close plot to exit')

    plt.tight_layout()
    plt.show(False)
