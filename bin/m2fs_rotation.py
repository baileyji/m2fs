#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from skimage.feature import peak_local_max
import argparse
from glob import glob
import ipdb
import os.path


def parse_cl():
    parser = argparse.ArgumentParser(description='Rotation Series Analyzer',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    #Finding Settings
    parser.add_argument('-x0', dest='x0', default=25,
                        action='store', required=False, type=int,
                        help='width of x fitting window')
    parser.add_argument('-y0', dest='y0', default=25,
                        action='store', required=False, type=int,
                        help='Width of y fitting window')
    parser.add_argument('-xw', dest='xw', default=4,
                        action='store', required=False, type=float,
                        help='X wid of box')
    parser.add_argument('-yw', dest='yw', default=4,
                        action='store', required=False, type=float,
                        help='y wid of box')

    #Program Settings
    parser.add_argument('-d', dest='dir', default='./',
                        action='store', required=False, type=str,
                        help='Input file directory')
#    parser.add_argument('--plot', dest='plot', default=False,
#                        action='store_true', required=False,
#                        help='Show many plots')
    parser.add_argument('--debug', dest='debug', default=False,
                        action='store_true', required=False,
                        help='Debugging prompts')
                        
    return parser.parse_args()

def find_rotation(files, args):
    
    #Finding settings
    xslice=slice(args.x0,args.x0+args.xw)
    yslice=slice(args.y0,args.y0+args.yw)

    rotation_data=[]
    for f in files:

        print('Processing {}'.format(f))

        #fetch image
        header=fits.getheader(f)
        im=fits.getdata(f)

        #Get the counts
        rotation_data.append(header['ROTATOR'], im[xslice,yslice].sum())

    rotation_data=np.array(rotation_data)

    cx=np.polyfit(rotation_data[:,0], rotation_data[:,1], 2)
    min_x=-cx[1]/cx[0]/2
    print 'Processed {}'.format(files)
    print('Best rotation @ {:.1f}'.format(min_x))

    for r in rotation_data: print('Focus: {}  {}.'.format(*r))

    if args.debug: ipdb.set_trace()

if __name__ =='__main__':
    
    args=parse_cl()
    
    files=(glob(os.path.join(args.dir,'*.fits'))+
           glob(os.path.join(args.dir,'*.fits.gz')))

    print('Running on {}'.format(files))
    

    find_rotation(rfiles,args)

    raw_input('Any key to exit')

