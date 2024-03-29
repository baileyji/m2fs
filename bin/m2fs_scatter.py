#!/usr/bin/env python
from __future__ import print_function
from astropy.io import fits
from m2fs.obs.scatter import mkscatter
import numpy as np
import argparse


def parse_cl():
    parser = argparse.ArgumentParser(description='Scattered light removal', add_help=True)
    parser.add_argument('-d', '--debug', dest='debug',
                        action='store_true', required=False,
                        help='Debug mode', default=False)
    #    parser.add_argument('-l','--listfile', dest='listfile', default='',
    #                     action='store', required=True, type=str,
    #                     help='testfile with files to proces')
    parser.add_argument('-f', dest='file', default='',
                        action='store', required=True, type=str,
                        help='File to proces')
    parser.add_argument('-a', dest='add', default=False,
                        action='store_true', required=False,
                        help='Add scatter')
    parser.add_argument('-s', dest='silent', default=False,
                        action='store_true', required=False,
                        help='Work quietly')
    parser.add_argument('-t', dest='threshold', default=1.5,
                        action='store', required=False, type=float,
                        help='Scattering threshold')
    parser.add_argument('-g', '--noglow', dest='glow', default=True,
                        action='store_false', required=False,
                        help='Do not model amplifier glow')
    parser.add_argument('-o', '--offset', dest='offset', default=0,
                        action='store', required=False, type=int,
                        help='Shift scatter regions up/down')
    parser.add_argument('--glowonly', dest='glowonly', default=False,
                        action='store_true', required=False,
                        help='Only model the glow')
    #    parser.add_argument('-o', dest='outdir', default='./',
    #                     action='store', required=False, type=str,
    #                     help='Out directory')

    args = parser.parse_args()
    return args


def add_scatter(f, debug=False, plot=True, thresh=.7, glow=True, offset=0,
                glowOnly=False, noprompt=False):
    hdul = fits.open(f, mode='update')
    s_im, err, (s_model, im_masked, glowout) = mkscatter(hdul[1].data, hdul[2].data,
                                                         plot=plot,
                                                         debug=debug,
                                                         scatter_thresh=thresh,
                                                         header=hdul[1].header,
                                                         do_glow=glow,
                                                         offset=offset,
                                                         glowOnly=glowOnly)

    if not noprompt:
        try:
            user_resp = raw_input('Ok? ').lower() if plot or debug else 'y'
        except:
            user_resp = input('Ok? ').lower() if plot or debug else 'y'
    else:
        user_resp = 'y'
    if user_resp in ['y', 'yes']:
        hdul[1].data -= s_im
        # hdul[2].data+=err**2
        hdul.insert(3, fits.ImageHDU(s_im.astype(np.float32), name='scatter'))
        hdul[3].header['SVAR'] = err ** 2
        hdul[3].header['THRESH'] = thresh
        hdul[3].header['GLOW'] = str(glow)
        hdul[3].header['GLOWONLY'] = str(glowOnly)
        hdul[3].header['OFFSET'] = str(offset)
    hdul.close()


if __name__ == '__main__':
    import sys

    args = parse_cl()
    file = args.file

    if args.add:
        add_scatter(file, debug=args.debug, plot=not args.silent,
                    thresh=args.threshold, glow=args.glow, offset=args.offset,
                    glowOnly=args.glowonly, noprompt=args.silent)
    else:
        print('Must be called with -a')

#        if len(sys.argv) >2:
#            ofile=sys.argv[2]
#        else:
#            print ' Output file not specified'
#            raise ValueError
#        s_im, err =mkscatter(fits.getdata(file,1),fits.getdata(file,2),
#                             scatter_thresh=args.threshold,
#                             debug=args.debug,
#                             plot=not args.silent, header=fits.getheader(file),
#                             do_glow=args.glow, offset=args.offset,
#                             glowOnly=args.glowonly)
#
#        if ofile==file:
#            ofile+='.scatter.fits.gz'
#
#        hdu=fits.PrimaryHDU(s_im.astype(np.float32))
#        hdu.writeto(ofile)
