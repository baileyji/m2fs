from logging import getLogger
import os
from glob import glob
import shutil
from datetime import datetime
import requests
import multiprocessing as mp


def compute_download_paths(file_ids, download_dir, url='https://m2fs.astro.lsa.umich.edu/data'):
    """ files is list of ('utYYYMMDD', {rb}####) """
    # Dates extend before and after runs by a good bit
    pathmap = {'Aug2013': (datetime(2013, 8, 1), datetime(2013, 8, 31)),
               'Nov2013': (datetime(2013, 11, 1), datetime(2013, 12, 2)),
               'Feb2014': (datetime(2014, 2, 1), datetime(2014, 3, 1)),
               'Jun2014': (datetime(2014, 6, 1), datetime(2014, 6, 30)),
               'Sep2014': (datetime(2014, 9, 1), datetime(2014, 9, 30)),
               'Dec2014': (datetime(2014, 12, 1), datetime(2014, 12, 31)),
               'Feb2015': (datetime(2015, 2, 1), datetime(2015, 3, 10)),
               'Apr2015': (datetime(2015, 4, 1), datetime(2015, 4, 28)),
               'Jul2015': (datetime(2015, 7, 1), datetime(2015, 7, 30)),
               'Sep2015': (datetime(2015, 9, 1), datetime(2015, 9, 30)),
               'Nov2015': (datetime(2015, 11, 1), datetime(2015, 11, 30)),
               'Feb2016': (datetime(2016, 2, 1), datetime(2016, 2, 28)),
               'Jun2016': (datetime(2016, 6, 1), datetime(2016, 7, 2)),
               'AugSep2016': (datetime(2016, 8, 1), datetime(2016, 9, 15)),
               'NovDec2016': (datetime(2016, 11, 1), datetime(2016, 12, 15)),
               'FebMar2017': (datetime(2017, 2, 1), datetime(2017, 3, 10)),
               'MayJun2017': (datetime(2017, 5, 1), datetime(2017, 6, 10)),
               'Jul2017': (datetime(2017, 7, 1), datetime(2017, 7, 31)),
               'Sep2017': (datetime(2017, 9, 1), datetime(2017, 9, 30)),
               'Nov2017': (datetime(2017, 11, 1), datetime(2017, 11, 30)),
               'Feb2018': (datetime(2018, 2, 1), datetime(2018, 2, 28)),
               'May2018': (datetime(2018, 5, 1), datetime(2018, 5, 31)),
               'Aug2018': (datetime(2018, 8, 1), datetime(2018, 8, 31)),
               'NovDec2018': (datetime(2018, 11, 1), datetime(2018, 12, 10)),
               'FebMar2019': (datetime(2019, 2, 1), datetime(2019, 3, 10)),
               'MayJun2019': (datetime(2019, 5, 1), datetime(2019, 6, 10)),
               'AugSep2019': (datetime(2019, 8, 1), datetime(2019, 9, 10))}

    to_fetch = []
    for id in file_ids:
        ut, seqstr = id
        run = ''
        dt = datetime.strptime(ut, 'ut%Y%m%d')
        for k, v in pathmap.items():
            if v[0] < dt < v[1]:
                run = k
                break

        if not run:
            getLogger(__name__).warning("Unable to determine download location for {}".format(file_ids))
            continue
        to_fetch.extend([('{}/{}/{}/{}c{}.fits'.format(url, run, ut if run != 'AugSep2016' else 'ALL', seqstr, i),
                          os.path.join(download_dir, ut, '{}c{}.fits'.format(seqstr,i))) for i in range(1, 5)])
    return to_fetch


def download(files, nproc=1):
    getLogger(__name__).info('Downloading files from m2fs.astro.lsa.umich.edu')
    paths = compute_download_paths(files)
    if nproc > 1:
        pool = mp.pool.Pool(nproc)
        pool.map(fetch_url, paths)
        pool.join()
    else:
        for p in paths:
            fetch_url(p)


def fetch_url(x):
    uri, path = x
    r = requests.get(uri, stream=True, auth=('m2fs', 'm2fsdata'))
    if r.status_code != 200:
        getLogger(__name__).error('m2fs.astro.lsa.umich.edu returned {} for {}, skipping'.format(r.status_code, uri))
        return

    total_size = int(r.headers.get('content-length', 0))  # Total size in bytes.
    if os.path.exists(path):
        if os.stat(path).st_size == total_size:
            return
        else:
            getLogger(__name__).info('{} exists with mismatched size, downloading'.format(path))

    getLogger(__name__).debug('Fetching {} to {}'.format(uri, path))

    directory = os.path.dirname(path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(path, 'wb') as f:
        for data in r.iter_content(1024):
            f.write(data)


def retrieve(file_id, local, remote=None, copy=True):
    path_end = os.path.join(file_id[0], file_id[1]+'.fits')
    localfile = os.path.join(local, path_end)
    if os.path.exists(localfile):
        return localfile
    if os.path.exists(localfile+'.gz'):
        return localfile+'.gz'

    if remote:
        remotefile = os.path.join(remote, path_end)
        if copy:
            if not os.path.exists(os.path.dirname(localfile)):
                os.mkdir(os.path.dirname(localfile))
            shutil.copy2(remotefile, os.path.dirname(localfile))
        return localfile if copy else remotefile


def find_seqno(seqno, dir):
    files = [os.path.join(utdir, f)
             for utdir, dirnames, files in os.walk(dir)
             for f in files
             if 'c1.fits' in f and f not in ('rc1.fits', 'bc1.fits') and not f.startswith('.')]
