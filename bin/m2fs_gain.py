#This is all kindof garbage, it is predicated on uniform illumination of the CCD
#which we don't have because of 1) solar lines, 2) spatial profile,
#3) fiber throughput variations
# This needs to be redone with a proper flat.



#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#from astropy.io import fits
#from scipy.optimize import curve_fit
#from skimage.feature import peak_local_max
#import argparse
#import ipdb
#import os.path
#import re
#
##
##def parse_cl():
##    parser = argparse.ArgumentParser(description='S/N Check',
##                                     add_help=True,
##                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
##                                 
##    parser.add_argument('file', metavar='FILE', type=str,
##                     help='File to process')
##    
##    #Finding settings
##    parser.add_argument('-d', dest='sep', default=10,
##                        action='store', required=False, type=int,
##                        help='Minimum delta between order peaks')
##    parser.add_argument('-m', dest='min', default=12,
##                        action='store', required=False, type=int,
##                        help='Minimum peak intensity')
##    parser.add_argument('-w', dest='xw', default=25,
##                        action='store', required=False, type=int,
##                        help='Width around peak to use for fitting gaussian')
##
##    #Measurement selection settings
##    parser.add_argument('-s', dest='hwid', default=8.0,
##                        action='store', required=False, type=int,
##                        help='Width to either side of smash')
##
##    parser.add_argument('-c', dest='min_sn', default=0,
##                        action='store', required=False, type=int,
##                        help='Minimum S/N cut')
##
##    parser.add_argument('--tpdat', dest='tpdat', default=False,
##                        action='store_true', required=False,
##                        help='give numbers needed for throughput checking')
##
##    return parser.parse_args()
##
##
#
#def get_bias_corrected_quad(file):
#    hdul=fits.open(file)
#    header=hdul[0].header
#    qdata=hdul[0].data.copy()
#    
#    biassec=header['BIASSEC']
#    bias=[int(s) for s in re.findall(r'\d+',biassec)]
#    
#    trimsec=header['TRIMSEC']
#    crop=[int(s) for s in re.findall(r'\d+',trimsec)]
#    
#    #Compute & subtract median for overscan region rowwise
#    biaslevels=np.median(qdata[crop[2]-1:,bias[0]-1:bias[1]], axis=1)
#    qdata-=biaslevels[:,np.newaxis]
#    
#    #Compute & subtract median bias row
#    qdata-=np.median(qdata[bias[2]:,:], axis=0)
#    
#    #Crop the image
#    qdata=qdata[crop[2]-1:crop[3],crop[0]-1:crop[1]]
#
#    return qdata.astype(float), hdul[0].data.astype(float)
#
#def gather_data(file_pairs, use_mask, average=np.median, clevel=True):
#    var=[]
#    sig=[]
#    for f1,f2 in file_pairs:
#        i1,_=get_bias_corrected_quad(f1)
#        i2,_=get_bias_corrected_quad(f2)
#
#        mean_sig=0.5*(average(i1[use_mask])+average(i2[use_mask]))
#        
#        signal_difference=np.abs(1-average(i1[use_mask])/mean_sig)
#        print 'Mean signal difference is {:.3f}%'.format(signal_difference*100)
#        
#        if clevel:
#            diff=(i1-average(i1[use_mask])*i2/average(i2[use_mask]))[use_mask]
#        else:
#            diff=(i1-i2)[use_mask]
#        
#        mad_variance=np.median(np.abs(diff-np.median(diff)))**2/2.0
#        std_variance=diff.std()**2/2.0
#        var_difference=(std_variance-mad_variance)/mad_variance
#        print 'MAD var and Std var differe by {:.3}%'.format(var_difference*100)
#        var.append(std_variance)
#        sig.append(mean_sig)
#
#    return np.array(sig),np.array(var)
#
#def find_gain(file_pairs, use_mask, average=np.median, clevel=False):
#    gain=[]
#    sig=[]
#    for f1,f2 in file_pairs:
#        i1=get_bias_corrected_quad(f1)
#        i2=get_bias_corrected_quad(f2)
#        
#        mean_sig=0.5*(average(i1[use_mask])+average(i2[use_mask]))
#        
#        signal_difference=np.abs(1-average(i1[use_mask])/mean_sig)
#        print 'Mean signal difference is {:.3f}%'.format(signal_difference*100)
#        
#        if clevel:
#            diff=(i1-average(i1[use_mask])*i2/average(i2[use_mask]))[use_mask]
#        else:
#            diff=(i1-i2)[use_mask]
#        gain.append(mean_sig/diff.std()**2)
#        sig.append(mean_sig)
#    
#    return np.array(sig),np.array(var)
#
#min_cnts=30
#
#dir='/Volumes/mobile/obsdata/feb2015/ut20150305/'
#paris=((1675,1676),
#       (1678,1679),
#       (1680,1681),
#       (1682,1683),
#       (1684,1685),
#       (1686,1687),
#       (1688,1689),
#       (1690,1691),
#       (1692,1693))
##       (1695,1696))
#
#for side in 'rb':
#    for quad in '1234':
#        files=[]
#        for seq_pair in paris:
#            f1='{}{:04}c{}.fits'.format(side,seq_pair[0],quad)
#            f2='{}{:04}c{}.fits'.format(side,seq_pair[1],quad)
#            files.append((os.path.join(dir,f1),os.path.join(dir,f2)))
#        
#        #Determine usable mask
#        use_mask=get_bias_corrected_quad(files[0][0])>min_cnts
#
#        #get data
#        sig,var=gather_data(files, use_mask, np.mean)
#        sig,gain=find_gain(files, use_mask, np.mean)
#
#        #determine gain & readnoise
#        par=np.polyfit(var,sig,1)
#        plt.plot(var,sig,'o')
#        plt.plot(np.arange(0,65000),np.poly1d(par)(np.arange(0,65000)))
#        plt.xlabel('Variance')
#        plt.ylabel('Signal')
#        g,rn=par
##        g,rn=compute_gain_and_rn(sig,var)
#
#        print '{}{}  Gain: {:.3f}  RN:{:.3f}'.format(side,quad,g,rn)
#
#
#
#
#
#ims=np.array([(get_bias_corrected_quad(f1)[0][use_mask].ravel(),
#               get_bias_corrected_quad(f2)[0][use_mask].ravel())
#               for f1,f2 in files])
#
#average=np.median
#mean_sig=average(ims,axis=2)
#sig_corr=mean_sig[:,0]/mean_sig[:,1]
#diff=np.diff(ims, axis=1)[:,0,:]
#good=np.abs(diff) < diff.std(1)[:,np.newaxis]*sigma_lim
#good=good.all(0)
#diff=diff[:,good]
#
#for
#        mean_sig=0.5*(average(i1[use_mask])+average(i2[use_mask]))
#        
#        signal_difference=np.abs(1-average(i1[use_mask])/mean_sig)
#        print 'Mean signal difference is {:.3f}%'.format(signal_difference*100)
#        
#        if clevel:
#            diff=(i1-average(i1[use_mask])*i2/average(i2[use_mask]))[use_mask]
#        else:
#            diff=(i1-i2)[use_mask]
#        
#        mad_variance=np.median(np.abs(diff-np.median(diff)))**2/2.0
#        std_variance=diff.std()**2/2.0
#        var_difference=(std_variance-mad_variance)/mad_variance
#        print 'MAD var and Std var differe by {:.3}%'.format(var_difference*100)
#        var.append(std_variance)
#        sig.append(mean_sig)
#
#    return np.array(sig),np.array(var)
#
#
#
#
#if __name__ =='__main__':
#    
#    import matplotlib as mpl
#    mpl.rcParams['font.size']=16
#    mpl.rcParams['lines.linewidth'] = 2
#    mpl.rcParams['axes.linewidth'] = 2
#    mpl.rcParams['ytick.major.size']=6
#    mpl.rcParams['ytick.major.width']=2
#    mpl.rcParams['xtick.major.size']=6
#    mpl.rcParams['xtick.major.width']=2
#
#    args=parse_cl()
#    
#    compute_sn(args.file,args)
#
#    raw_input('Any key to exit')
