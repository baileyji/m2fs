#!/usr/bin/python
from __future__ import print_function
import argparse
from glob import glob
import ipdb
import os.path
import scipy.ndimage

try:
    input = raw_input
except NameError:
    pass

import m2fs
from m2fs.focus import *
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
    parser.add_argument('-min', dest='minv', default=50,
                        action='store', required=False, type=int,
                        help='Minimum spot intensity')
    parser.add_argument('-max', dest='maxv', default=2000,
                        action='store', required=False, type=int,
                        help='Maximum spot intensity')
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
    parser.add_argument('-x_min', dest='x_min', default=2,
                        action='store', required=False, type=float,
                        help='Minimum PSFx')
    parser.add_argument('-x_max', dest='x_max', default=25.0,
                        action='store', required=False, type=float,
                        help='Maximum PSFx')
    parser.add_argument('-y_min', dest='y_min', default=5.5,
                        action='store', required=False, type=float,
                        help='Minimum PSFx')
    parser.add_argument('-y_max', dest='y_max', default=25.0,
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




def find_focus(files, args):
    # Finding settings
    min_sep = args.min_sep
    minv = args.minv
    maxv = args.maxv
    max_peaks = args.max_peaks

    focus_data = {}
    for f in files:

        print('Processing {}'.format(f))

        # fetch image
        header = fits.getheader(f)
        focus = header['FOCUS']
        filt = header['FILTER']
        xbin=int(header['BINNING'][0])
        ybin=int(header['BINNING'][-1])

        if args.skip_binning:
            xbin=ybin=1

        # Fitting settings
        xw, yw = args.xw//xbin, args.yw//ybin  # Width of window around each line to model
        x_wid_guess = args.x_wid_guess/xbin/2.35482  # FWHM
        y_wid_guess = args.y_wid_guess/ybin/2.35482

        # Fit cut settings
        x_min, x_max = args.x_min/xbin/2.35482, args.x_max/xbin/2.35482
        y_min, y_max = args.y_min/ybin/2.35482, args.y_max/ybin/2.35482

        if os.path.basename(f)[0].lower() == 'b':
            side = 'b'
            tempi = [1, 3, 5, 7, 9]
        else:
            side = 'r'
            tempi = [1, 2, 4, 6, 8]
        temps = tuple(header['TEMP{:02}'.format(i)] for i in tempi)
        im = fits.getdata(f)

        # Find the peaks
        immed = scipy.ndimage.median_filter(im, 5)
        peaks = find_peaks(immed, min_sep=min_sep, minv=minv, maxv=maxv)

        # peaks are y coord x coord
        print('Found {} peaks between {} and {} separated by at least {}.'.format(
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
            plt.plot(peaks[:,1], peaks[:,0], 'w*', alpha=.5)
            plt.title('Image & Peaks')
            plt.show()

        # Fit the peaks
        ret = fit_peaks_v2(im, peaks, np.array((xw, yw),dtype=int), sigma_guess=(x_wid_guess, y_wid_guess),
                           debug=False, xwid_bound=(x_min, x_max), ywid_bound=(y_min, y_max))

        # Filter fits and construct results
        x, y, xwid, ywid, covar, modelim, modeldat = filter_peaks(im, ret,
                                                                  x_min=x_min, x_max=x_max,
                                                                  y_min=y_min, y_max=y_max,
                                                                  xw=xw, yw=yw)

        # Plot fit data & models
        if args.plot or args.debug:
            plt.figure(2)
            plt.clf()
            fig, ax = plt.subplots(1, 3, num=2, sharex=True, sharey=True)
            plt.sca(ax.flat[0])
            plt.imshow(modelim, vmin=vmin, vmax=vmax, origin='lower')
            plt.title('Line fits')
            plt.sca(ax.flat[1])
            plt.imshow(modeldat, vmin=vmin, vmax=vmax, origin='lower')
            plt.title('Line fit data')
            plt.sca(ax.flat[2])
            plt.imshow(modelim - modeldat, vmin=vmin, vmax=vmax, origin='lower')
            plt.title('Difference')
            plt.xlim(800, 1200)
            plt.ylim(0, 375)
            plt.colorbar()
            plt.show(block=False)

            # Plot width data
            plt.figure(3)
            plt.subplot(1, 3, 1)
            plt.hist(xwid, 100, histtype='stepfilled')
            plt.title('X FWHM: Avg {:.2f} Focus:{}'.format(xwid.mean(), focus))
            plt.xlim(0, 15)
            plt.subplot(1, 3, 2)
            plt.hist(ywid, 100, histtype='stepfilled')
            plt.title('Y FWHM: Avg {:.2f}'.format(ywid.mean()))
            plt.xlim(0, 15)
            plt.subplot(1, 3, 3)
            plt.hist(covar, 100, histtype='stepfilled')
            plt.title('Covar Avg {:.2f}'.format(covar.mean()))
            plt.show(block=True)

        if args.debug:
            val = input('Continue(any), Abort(a)>')
            if val == 'a':
                exit(1)

        print('Focus: {}. FWHM {:.1f}, {:.1f}. Covariance {:.1f}'.format(
              focus, xwid.mean(), ywid.mean(),covar.mean()))
        focus_data[focus] = (xwid.mean(), ywid.mean(), covar.mean(), temps)

    foc = list(focus_data.keys())
    xvals = [focus_data[f][0] for f in foc]
    yvals = [focus_data[f][1] for f in foc]
    cvals = [focus_data[f][2] for f in foc]

    temps = np.array([k[3][1] for k in focus_data.values()])

    if len(foc) > 2:
        cx = np.polyfit(foc, xvals, 2)
        cy = np.polyfit(foc, yvals, 2)
        min_x = -cx[1] / cx[0] / 2
        min_y = -cy[1] / cy[0] / 2
        print('Processed {}'.format(files))
        print('Filter: {}'.format(filt))
        print('Best x focus @ {:.1f} with value of {:.2f} at temp {:3f}'.format(
            min_x, np.poly1d(cx)(min_x), temps.mean()))
        print('Best y focus @ {:.1f} with value of {:.2f}'.format(min_y,
                                                                  np.poly1d(cy)(min_y)))
        print('{:.1f} {:.1f} {} {} {} {} {}'.format(min_x, min_y,
                                                    *focus_data[f][3]))
    for i, f in enumerate(foc):
        print('Focus: {} FWHM: {:.1f}, {:.1f} Covar: {:.1f} Temps:{}'.format(
            foc[i], xvals[i], yvals[i], cvals[i], focus_data[f][3]))

    if len(foc) > 2:
        plt.figure(6)
        plt.plot(foc, xvals, 'bo', label='PSF x')
        plt.plot(foc, yvals, 'ro', label='PSF y')
        plt.plot(foc, cvals, 'go', label='Covar')

        extent = np.abs(np.array(foc) - min_x).max()
        xx = np.linspace(min(min(foc), min_x - extent), max(max(foc), min_x + extent),
                         100)
        plt.plot(xx, np.poly1d(cx)(xx), 'b')
        plt.plot(xx, np.poly1d(cy)(xx), 'r')
        plt.text(min_x, 1.5, 'Temp: {:.3f}\nFilter: {}'.format(temps.mean(), filt),
                 color=side)
        plt.legend()

    if args.debug:
        ipdb.set_trace()


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
        find_focus(rfiles, args)
        plt.show(block=True)
    if bfiles:
        find_focus(bfiles, args)
        plt.show(block=True)

    input('Any key to exit')
