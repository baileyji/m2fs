#!/usr/bin/env python
import numpy as np
import argparse
import uuid
from subprocess import call
import m2fs.obs


DS9_APP='/Applications/DS9.app'

def parse_cl():
    parser = argparse.ArgumentParser(description='Help undefined',
                                     add_help=True)
    parser.add_argument('-d','--dir', dest='dir',
                        action='store', required=False, type=str,
                        help='',default='./')
    parser.add_argument('-s','--side', dest='side',
                        action='store', required=True, type=str,
                        help='Side to watch',default='r')
    parser.add_argument('--raw', dest='raw',
                        action='store', required=False, type=boolean,
                        help='Do not process data',default=False)
    return parser.parse_args()


if __name__ =='__main__':
    args=parse_cl()

    watchdir=args.dir
    side=args.side.lower()
    
    tmpdir=str(uuid.uuid4())
    
    displayed=''
    
    while True:
        time.sleep(10)

        files=[f for f in glob(watchdir+side+'*.fits')]
        ctimes=[os.path.getctime(f) for f in new_files]
        newest=np.array(self.ctimes).argmin()

        latest=files[newest]

        try:
            seqno=m2fs.obs.info(latest).seqno
        except IOError:
            continue

        merged_name='{}{frameno:04}.fits'.format(side, seqno)
        if displayed == merged name:
            continue

        try:
            m2fs.obs.merge.mergequad(seqno, side=self.side, odir=tmpdir,
                                     )
        except IOError:
            continue

        displayed=merged_name
        call('open {} -mecube {}'.format(DS9_APP, tmpdir+'/'+merged_file)
      
                        

