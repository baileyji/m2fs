#!/usr/bin/python
from __future__ import print_function
import argparse
import pickle
from glob import glob
import ipdb
import os.path

import matplotlib.pyplot as plt
import scipy.ndimage

try:
    input = raw_input
except NameError:
    pass
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import m2fs
import m2fs.obs
from m2fs.util import find_peaks
from m2fs.focus import *
import m2fs.focus
from m2fs.seqno import *


def parse_cl():
    parser = argparse.ArgumentParser(description='Focus Series Analyzer',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-sbin', dest='skip_binning', default=False,
                        action='store', required=False, type=bool,
                        help='Skip factoring binning into args')
    # Finding settings
    parser.add_argument('-sep', dest='min_sep', default=10,
                        action='store', required=False, type=int,
                        help='Spot separation threshold')
    parser.add_argument('-min', dest='minv', default=None,
                        action='store', required=False, type=int,
                        help='Minimum spot intensity (defaults to 3x std the robust median)')
    parser.add_argument('-max', dest='maxv', default=None,
                        action='store', required=False, type=int,
                        help='Maximum spot intensity (defaults to 10x std of robust mean peaks)')
    parser.add_argument('-npeak', dest='max_peaks', default=4000,
                        action='store', required=False, type=int,
                        help='Maximum number of points to measure')
    # Fitting Settings
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

    # Measurement selection settings
    parser.add_argument('-x_min', dest='x_min', default=2.2,
                        action='store', required=False, type=float,
                        help='Minimum PSFx')
    parser.add_argument('-x_max', dest='x_max', default=22.0,
                        action='store', required=False, type=float,
                        help='Maximum PSFx')
    parser.add_argument('-y_min', dest='y_min', default=2.2,
                        action='store', required=False, type=float,
                        help='Minimum PSFx')
    parser.add_argument('-y_max', dest='y_max', default=22.0,
                        action='store', required=False, type=float,
                        help='Maximum PSFy')

    # Program Settings
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


def measure_focus(files, args):
    # Finding settings
    min_sep = args.min_sep
    minv = args.minv
    maxv = args.maxv
    max_peaks = args.max_peaks

    focus_data = {}
    for f in files:

        try:
            with open(f + '.focusdata.pickle','rb') as pf:
                d = pickle.load(pf)
            if d['args'] == vars(args):
                print('Using cache for {}'.format(f))
                focus_data[f] = d
                continue
        except IOError:
            pass
        except UnicodeDecodeError:
            pass

        print('Processing {}'.format(f))

        # fetch image
        ob = m2fs.obs.info(f)
        im = fits.getdata(ob.file)
        if args.minv is None:
            stats=sigma_clipped_stats(im)
            minv = max(stats[0] + 2*stats[2], ob.header['ENOISE']*2)
            print('Using a threshold value of {:.1f}'.format(minv))
        if args.maxv is None:
            p = find_peaks(im, min_sep=args.min_sep, minv=minv, maxv=None)
            allpeaks = im[p[:,0],p[:,1]]
            aps = sigma_clipped_stats(allpeaks)
            maxv = aps[0]+10*aps[2]
            print('Using a maximum value of {:.1f}'.format(maxv))
        xbin, ybin = ob.binning

        if args.skip_binning:
            xbin = ybin = 1

        # Fitting settings
        xw, yw = args.xw // xbin, args.yw // ybin  # Width of window around each line to model
        x_wid_guess = args.x_wid_guess / xbin / 2.35482  # FWHM
        y_wid_guess = args.y_wid_guess / ybin / 2.35482

        # Fit cut settings
        x_min, x_max = args.x_min / xbin / 2.35482, args.x_max / xbin / 2.35482
        y_min, y_max = args.y_min / ybin / 2.35482, args.y_max / ybin / 2.35482

        # Find the peaks
        immed = scipy.ndimage.median_filter(im, 5)
        peaks = find_peaks(immed, min_sep=min_sep, minv=minv, maxv=maxv)

        # peaks are y coord x coord
        print('Found {} peaks between {:.1f} and {:.1f} separated by at least {}.'.format(
            peaks.shape[0], minv, maxv, min_sep))

        if peaks.shape[0] > max_peaks:
            print('Using {} predictably random peaks of {}'.format(max_peaks, peaks.shape[0]))
            np.random.seed(0)
            peaks = peaks[np.random.choice(peaks.shape[0], max_peaks, replace=False)]

        vmin = max(10, minv * .9)
        vmax = min(100, 1.1 * maxv)

        if args.plot or args.debug:
            # Plot image and detected peaks
            plt.figure(1)
            plt.clf()
            plt.imshow(im, vmin=vmin, vmax=vmax, origin='lower')
            plt.plot(peaks[:, 1], peaks[:, 0], 'w*', alpha=.5)
            plt.title('Image & Peaks')

        # Fit the peaks
        ret = fit_peaks(im, peaks, np.array((xw, yw), dtype=int), sigma_guess=(x_wid_guess, y_wid_guess),
                        debug=False, xwid_bound=(x_min, x_max), ywid_bound=(y_min, y_max))

        # Filter fits and construct results
        x, y, xwid, ywid, covar, modelim, modeldat = filter_peaks(im, ret, (x_min*2.35482, x_max*2.35482),
                                                                  (y_min * 2.35482, y_max * 2.35482),
                                                                  xw=xw, yw=yw)

        # Plot fit data & models
        if args.plot or args.debug:
            plt.close(2)
            fig, ax = plt.subplots(1, 3, num=2, sharex=True, sharey=True, figsize=(14, 8))
            plt.sca(ax.flat[0])
            plt.imshow(modelim, vmin=vmin, vmax=vmax, origin='lower')
            plt.title('Line fits')
            # plt.colorbar()
            plt.sca(ax.flat[1])
            plt.imshow(modeldat, vmin=vmin, vmax=vmax, origin='lower')
            plt.title('Line fit data')
            # plt.colorbar()
            plt.sca(ax.flat[2])
            pltim=plt.imshow(modelim - modeldat, vmin=vmin, vmax=vmax, origin='lower')
            plt.title('Difference')
            plt.xlim(800, 1200)
            plt.ylim(0, 375)
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax.flat[2])
            plt.colorbar(pltim, cax=divider.append_axes("right", size="5%", pad=0.05))
            # plt.colorbar()
            plt.tight_layout()
            plt.suptitle(f)

            # Plot width data
            plt.close(3)
            fig, ax = plt.subplots(1, 3, num=3, figsize=(14, 8))
            plt.sca(ax.flat[0])
            plt.hist(xwid, 100, histtype='stepfilled')
            plt.title('FWHMx: {:.2f}'.format(xwid.mean()))
            # plt.xlim(0, 15)
            plt.sca(ax.flat[1])
            plt.hist(ywid, 100, histtype='stepfilled')
            plt.title('FWHMy: {:.2f}'.format(ywid.mean()))
            # plt.xlim(0, 15)
            plt.sca(ax.flat[2])
            plt.hist(covar, 100, histtype='stepfilled')
            plt.title('Covar {:.2f}'.format(covar.mean()))
            plt.suptitle(f+' {}x{} bin'.format(xbin,ybin))
            plt.tight_layout()
            plt.savefig(f+'fwhm_dist.pdf')

        if args.debug:
            val = input('Continue(any), Abort(a)>')
            if val == 'a':
                exit(1)

        print('FWHM {:.1f}, {:.1f}. Covariance {:.1f}'.format(xwid.mean(), ywid.mean(), covar.mean()))
        d = dict(file=f, xwid=xwid, ywid=ywid, covar=covar, args=vars(args))
        with open(f + '.focusdata.pickle', 'wb') as pf:
            pickle.dump(d, pf, protocol=pickle.HIGHEST_PROTOCOL)
        focus_data[f] = d
        plt.show()
    return focus_data


def find_focus(focus_data):
    foc, xvals, yvals, cvals, temps = [],[],[],[],[]

    for file, foc_data in focus_data.items():
        ob = m2fs.obs.info(file)
        foc.append(ob.header['FOCUS'])
        filt = ob.header['FILTER']
        side = ob.side
        xvals.append(foc_data['xwid'].mean())
        yvals.append(foc_data['ywid'].mean())
        cvals.append(foc_data['covar'].mean())
        temps.append(getattr(ob.temps, 'cradle_'+ob.side))

    foc, xvals, yvals, cvals, temps = map(np.array, (foc, xvals, yvals, cvals, temps))
    if len(foc) < 2:
        return

    cx = np.polyfit(foc, xvals, 2)
    cy = np.polyfit(foc, yvals, 2)
    min_x = -cx[1] / cx[0] / 2
    min_y = -cy[1] / cy[0] / 2

    print('Processed {}'.format(files))
    print('Filter: {}'.format(filt))
    print('Best x focus @ {:.1f} with value of {:.2f} at temp {:3f}'.format(
        min_x, np.poly1d(cx)(min_x), temps.mean()))
    print('Best y focus @ {:.1f} with value of {:.2f}'.format(min_y, np.poly1d(cy)(min_y)))

    plt.figure(6+(side=='r'))
    plt.plot(foc, xvals, 'o', color='C0', label='PSF x')
    plt.plot(foc, yvals, 'o', color='C1', label='PSF y')
    plt.plot(foc, cvals, 'o', color='C2', label='Covar')

    extent = np.abs(foc - min_x).max()
    xx = np.linspace(min(foc.min(), min_x - extent), max(foc.max(), min_x + extent), 100)
    plt.plot(xx, np.poly1d(cx)(xx), 'C0')
    plt.plot(xx, np.poly1d(cy)(xx), 'C1')
    plt.text(min_x, 1.5, 'Temp: {:.3f}\nFilter: {}'.format(temps.mean(), filt), color=side)
    if side=='r':
        plt.axvline(m2fs.focus.RFOCUS(temps.mean()))
    else:
        plt.axvline(m2fs.focus.BFOCUS(temps.mean()))
    plt.legend()


    print('Focus, FWHMx, FWHMy, Covar, Temp')
    for f, x, y, c, t in zip(foc, xvals, yvals, cvals, temps):
        print('{}, {:.1f}, {:.1f}, {:.1f}, {}'.format(f, x, y, c, t))


if __name__ == '__main__':

    args = parse_cl()

    if not args.files:
        files = (glob(os.path.join(args.dir, '*.fits')) +
                 glob(os.path.join(args.dir, '*.fits.gz')))

        try:
            seqno = get_seqnos(args.listfile)
            info = [m2fs.obs.info(f, no_load=True) for f in files]

            files = []
            for i in info:
                for s in seqno:
                    if i.seqno_match(s):
                        files.append(i.file)
                        break
        except IOError:
            print('No listfile, doing all')
    else:
        files = [os.path.join(args.dir, x) for x in args.files.split(' ')]

    print('Running on {}'.format(files))

    rfiles = [f for f in files if os.path.basename(f)[0] == 'r']
    bfiles = [f for f in files if os.path.basename(f)[0] == 'b']
    if rfiles:
        data = measure_focus(rfiles, args)
        find_focus(data)
        plt.savefig('r_focus.pdf')
        plt.show(block=False)
    if bfiles:
        data = measure_focus(bfiles, args)
        find_focus(data)
        plt.savefig('b_focus.pdf')

    plt.show(block=True)
    input('Any key to exit')
