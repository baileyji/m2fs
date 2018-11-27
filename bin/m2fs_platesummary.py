#!/usr/bin/python
import numpy as np
from hole_mapper import plate
from hole_mapper import pathconf
import sys

def floatable(x):
    try:
        float(x)
    except ValueError:
        return False
    return True

def summarize_plates():
    """
    Grab all the plates in the plate directory and print out info.
    
    List name, number of fields, number of holes
        list field name, RA/Dec, number of targs, number of skys, minsky,mustkeep, status min/max/mean mags
        
    Grab all the setups from the setup files
    for each setup file
        list setup file name, number of setups
            for each setup in file
                list setup name, setup settings, if a fibermap is present
                
    """
    pnames=plate.get_all_plate_names()
    nlen=min(max(map(len,pnames)), 25)
    namefmt='{:'+'{}'.format(nlen)+'}'
    lines=[]
    for pname in pnames:
        p=plate.get_plate(pname)
        lines.append((namefmt+' {} Field(s), {:4} holes').format(pname, len(p.fields),
                                                        len(p.all_holes)))
        
        
        #field name string length
        fnames=[f.name for f in p.fields]
        nlen=min(max(map(len,fnames)), 25)
        fnamefmt='{:'+'{}'.format(nlen)+'}'
        
        #magnitude keys
        magkey=list(set(k for f in p.fields for t in f.targets for k in t.user.keys()
                        if k in ['v','b','r','g','i'] or 'mag' in k))

        indent='   '
        fmt=(indent+fnamefmt+'  {}  {:4}  {:3}  {:3}  {}'+'  {:4}'*len(magkey))
        lines.append(fmt.format('Field', 'RA Dec', 'N Targ', 'N Sky', 'minsky',
                         'mustkeep', *magkey))
                         


        for f in p.fields:
            mags=[np.array([float(t.user[k]) for t in f.targets
                            if k in t.user and floatable(t.user[k])]) for k in magkey]
             
            for m in mags:m[(m<5) | (m > 30)]=np.nan
            meanmag=[np.nanmean(m) for m in mags]

            meanmag=['{:2.1f}'.format(m) for m in meanmag]
        
            lines.append(fmt.format(f.name, f.info['(ra, dec)'],
                  len(f.targets), len(f.skys), f.info.get('minsky',0),
                  f.info.get('mustkeep', False), *meanmag))


    return lines

#
#def iobserve_file()
#    pnames=plate.get_all_plate_names()
#    nlen=min(max(map(len,pnames)), 25)
#    namefmt='{:'+'{}'.format(nlen)+'}'
#    lines=[]
#    for pname in pnames:
#        p=plate.get_plate(pname)
#        lines.append(fmt.format(p.))


if __name__ == '__main__':
    
    if len(sys.argv) >1 :
        pathconf.ROOT=sys.argv[1]
#    import ipdb;ipdb.set_trace()
    lines=summarize_plates()
    for line in lines:
        print line


    
