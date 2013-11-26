#!/usr/bin/python
import scipy.interpolate
import numpy as np
from astropy.stats import sigma_clip

def m2fs_scatter(im):

    nrows,ncols=im.shape
    regions=[(0,50), (480,575), (1010,1100), (1525,1600), (2035,2120),
             (2540,2620), (3030,3110), (3525, 3600), (3995 ,4111)]


    #compute the mean scattered light in blocks of xbinwid in each
    # of the scattered light regions between the tetri 
    x=[]
    y=[]
    z=[]
    xbinwid=64
    edgecase=False
    for r in regions:
        for b in xrange(0,4095,xbinwid):
            blockmean=np.median(sigma_clip(im[r[0]:r[1],b:b+xbinwid]))
            z.append(blockmean)
            x.append(0.5*(2*b+xbinwid))
            y.append(0.5*(r[1]+r[0]))
            if b==0:
                edgecase|=True
                edgex=0
            elif b+xbinwid >= ncols:
                edgecase|=True
                edgex=ncols-1
            else:
                edgex=0.5*(2*b+xbinwid)
            if r==regions[0]:
                edgecase|=True
                edgey=0
            elif r==regions[-1]:
                #import ipdb;ipdb.set_trace()
                edgecase|=True
                edgey=nrows-1
            else:
                 edgey=0.5*(r[1]+r[0])
            if edgecase:
                edgecase=False
                z.append(blockmean)
                x.append(edgex)
                y.append(edgey)

#    scatter_data=im.copy()
#    for rs in zip(regions[0:],regions[1:]):
#        scatter_data[rs[0][1]+1:rs[1][0],:]=0

    #Setup for interpolation
    arrxyin=np.zeros((len(x),2),dtype=np.int64)
    arrxyin[:,0]=x #columns
    arrxyin[:,1]=y

    fill_xy=np.meshgrid(np.arange(nrows),np.arange(ncols),indexing='ij')
    arrxyout=np.zeros((len(fill_xy[0].flatten()),2),dtype=np.int64)
    arrxyout[:,0]=fill_xy[1].flatten()
    arrxyout[:,1]=fill_xy[0].flatten()
    
    #Interpolate
    scatter=scipy.interpolate.griddata(arrxyin, np.array(z),
                                    arrxyout, method='cubic')

    #scattered light map
    scatter=scatter.reshape(im.shape)
    
    #imshow(np.flipud(scatter))
    
    return scatter


if __name__ == '__main__':
    import sys
    hdu=fits.read(sys.argv[1])[0]
    header=hdu.header
    im=hdu.data
    im-=m2fs_scatter(im)
    out=PrimaryHDU(im)
    out.header=header
    out.writeto(sys.argv[2])
