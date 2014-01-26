
def info(data):
    """data==str -> fitsfilename"""
    if type(data)==str:
        file=data
    
    import os.path
    ret={}
    if 'r' in os.path.basename(file):
        ret['side']='r'
    else:
        ret['side']='b'
    ret['seqno']=(os.path.basename(file)[1:5])

    return ret

from m2fs.obs.scatter import *
