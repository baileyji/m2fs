from logging import getLogger

import astropy.io.fits
import numpy as np
from jbastro.misc import derangify, rangify
import os.path
from glob import glob
from collections import namedtuple
from astropy.io import fits
from datetime import datetime
from astropy.time import Time, TimeDelta
from astropy.coordinates import Latitude, Longitude
from astropy import units as u
from typing import List




class SeqStr:
    def __init__(self, s: str, default_cr=False, ut=None, run_id=''):
        """[r|b|]#-#,#... {docr|nocr|lamp|arc|led|quartz|bias|dark}
        images: i.e. r0003
        ids: i.e. ut 20220921/r0003
        """
        self.line = s
        words = s.strip().lower().split()
        if words[0][0].isdigit():
            seq = words[0]
            side = 'rb'
        else:
            side = words[0][0]
            assert side in 'rb'
            seq = words[0].strip('rb')

        docr = default_cr
        if 'nocr' in words:
            docr = False
        elif 'docr' in words:
            docr = True

        self.side = side
        self.seq_id = words[0]
        self.seq_str = seq
        self.nums = derangify(seq)
        self.docr = docr
        self.lamp = any(x in words for x in ('lamp', 'arc', 'th', 'xe', 'he', 'ne', 'thxe', 'thar', 'ar', 'benear',
                                             'nearbe', 'arnebe', 'arbene', 'nebear', 'nearbe', 'lihe', 'heli', 'thne',
                                             'tharne', 'thnear', 'li'))
        self.twilight = any(x in words for x in ('twi', 'twilight'))
        self.ledflat = 'led' in words
        self.quartz = 'quartz' in words
        self.bias = 'bias' in words
        self.dark = 'dark' in words
        for w in (w for w in words if w.startswith('ut')):
            try:
                ut = datetime.strptime(w, '#%Y%m%d').strftime('ut%Y%m%d')
            except ValueError:
                pass
        self.ut = ut
        self.run_id = run_id
        self.string = ' '.join(words[1:])
        self.images = tuple(f'{side}{x:04}' for side in self.side for x in self.nums)
        self.rimages = tuple(i for i in self.images if i.startswith('r'))
        self.bimages = tuple(i for i in self.images if i.startswith('b'))

    def __eq__(self, other):
        if self.ut is not None:
            date = self.ut == other.id
        else:
            date = self.run_id == other.run_id
        return self.images == other.images and date

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(str(self.images) + str(self.ut or self.run_id))

    @property
    def ids(self):
        return tuple(os.path.join(self.run_id, x) for x in self.images)

    @property
    def ids_r(self):
        return tuple(os.path.join(self.run_id, x) for x in self.rimages)

    @property
    def ids_b(self):
        return tuple(os.path.join(self.run_id, x) for x in self.bimages)

    def __repr__(self):
        return f'<SeqStr: {self}>'

    def __str__(self):
        return f'{self.side}{self.seq_str}'


class MSpecArmFileset:  # This is a fileset for the quadrants!!!
    """ Set of mspec files since exposure id but may be r b or r and b"""

    def __init__(self, filename, ut=None, run_id='', raw_path='', merged_path=''):
        """
        filename: [r|b]#.fit{s}{.gz}
        Set ut to 'utyyyymmdd' only if files not in folder named by ut, may be ignored even then but exposures
        across runs may collide"""
        self.user_filename = filename
        self.user_file = os.path.basename(filename)
        self.user_path = os.path.dirname(filename)
        self.raw_path = raw_path
        self.merged_path = merged_path

        file, ext = os.path.splitext(self.user_file)
        self.gzipped = ext.lower() in ('.gz',)
        self.gzipped = '.gz' if self.gzipped else ''  # bool evals false if no gzip, concat for correct extension if T
        if self.gzipped:
            file, _ = os.path.splitext(file)

        self.file = file.lower()

        self.user_is_quadrant = self.file.endswith(('c1', 'c2', 'c3', 'c4'))
        if self.user_is_quadrant:
            if not self.raw_path:
                self.raw_path = self.user_path
        else:
            if not self.merged_path:
                self.merged_path = self.user_path

        if self.file.startswith('r'):
            self.side = 'r'
        elif self.file.startswith('b'):
            self.side = 'b'
        else:
            raise ValueError('filename must be sided')

        self.num = self.file.strip('rb')[:-2] if self.user_is_quadrant else self.file.strip('rb')
        self.num = int(self.num)

        self.run_id = run_id
        self.ut = ut
        ut = os.path.split(self.user_path.lower())[-1]
        if ut.startswith('ut'):
            self.ut = ut

        self.id = os.path.join(self.run_id or '', f'{self.side}{self.num:04}')

    def quadrant_files(self, gzipped=None, raw_path=''):
        if gzipped is False:
            gzipped = ''
        elif gzipped is True:
            gzipped = '.gz'
        elif self.user_is_quadrant:
            gzipped = self.gzipped
        else:
            gzipped = ''
        if not raw_path:
            raw_path = self.raw_path

        return tuple(os.path.join(raw_path, f'{self.side}{self.num:04}c{i}.fits{gzipped}') for i in range(1, 5))

    def merged_file(self, gzipped=None, merged_path=''):
        if gzipped is False:
            gzipped = ''
        elif gzipped is True:
            gzipped = '.gz'
        else:
            gzipped = self.gzipped

        if not merged_path:
            merged_path = self.merged_path

        return os.path.join(merged_path, f'{self.side}{self.num:04}.fits{gzipped}')

    def __repr__(self):
        return f'<{self.__class__.__name__}: {self}>'

    def __str__(self):
        return str(self.id)


class MSpecFits(MSpecArmFileset):
    def __init__(self, fileset, load=True, ut=None, run_id='', raw_path='', merged_path='', listfile_record=None):
        if isinstance(fileset, MSpecArmFileset):
            super().__init__(fileset.user_filename, ut=fileset.ut, run_id=fileset.run_id,
                             raw_path=raw_path or fileset.raw_path, merged_path=merged_path or fileset.merged_path)
        else:
            super().__init__(fileset, ut=ut, run_id=run_id, raw_path=raw_path, merged_path=merged_path)
        self._header = None
        self._listfile = None
        self.cr_settings = False
        self.listfile_record = listfile_record
        if load:
            self.load()

    @staticmethod
    def _compute_midpoint(head):
        try:
            time, offset = head['UT-MID'], TimeDelta(0.0, format='sec')
        except KeyError:
            time = head['UT-DATE'] + ' ' + head['UT-TIME']
            offset = TimeDelta(head['EXPTIME'] / 2.0, format='sec')

        return Time(time, format='iso', scale='utc', location=(Longitude(-70.69242, unit=u.degree),
                                                               Latitude(-29.01423, unit=u.degree))) + offset

    def __repr__(self):
        return f'<MSpecFits: {self}>'

    @property
    def dark(self):
        return self.header['EXPTYPE'] == 'dark'

    @property
    def bias(self):
        return self.header['EXPTYPE'] == 'bias'

    @property
    def lamp(self):
        return sum(self.header[k] for k in ('BENEAR', 'LIHE', 'THXE')) > 0

    @property
    def led(self):
        leds = ('LED-UV', 'LED-BL', 'LED-VI', 'LED-NR', 'LED-FR', 'LED-IR')
        return sum(self.header[k] for k in leds) > 0

    @property
    def quartz(self):
        return self.listfile_record.quartz or False

    @property
    def twilight(self):
        return self.listfile_record.twilight or False

    @property
    def temps(self):
        x = [self.header.cards[k] for k in self.header if k.startswith('TEMP')]
        names = [c.comment.lower() for c in x]
        names[names.index('ccd temperature [c]')] = 'ccd'
        temps = namedtuple('{}Temps'.format(self.instrument.upper()), names)
        return temps(*[c.value for c in x])

    @property
    def fitsfile(self):
        return self.merged_file()

    @property
    def seqno(self):
        return self.num

    @property
    def seqid(self):
        return self.id

    def seqno_match(self, seq_id):
        if not isinstance(seq_id, SeqStr):
            seq_id = SeqStr(seq_id)
        return self.id in seq_id.ids

    @property
    def header(self):
        if self._header is None:
            self.load()
        return self._header

    def load(self, reload=False):
        if self._header is not None and not reload:
            return
        try:
            f = self.merged_file()
            self._header = fits.getheader(f)
        except IOError:
            getLogger(__name__).error(f'Unable to load {f} for header')
            self._header = None

    @property
    def instrument(self):
        return self.header['INSTRUME'].lower()

    @property
    def midpoint(self):
        return MSpecFits._compute_midpoint(self.header)

    @property
    def binning(self):
        b = self.header['BINNING']
        return int(b[0]), int(b[-1])

    def merged(self, do_cosmic=None, overwrite=False, pool=4, repair=False,
               raw_path=None, merged_path=None):
        from m2fs.obs.merge import _mergequad
        x = _mergequad(self.quadrant_files(raw_path=raw_path),
                       self.merged_file(merged_path=merged_path),
                       do_cosmic=do_cosmic or self.docr, pool=pool, repair=repair, overwrite=overwrite)
        return x

    @property
    def science_hdu(self) -> astropy.io.fits.ImageHDU:
        return fits.open(self.merged_file())['SCIENCE']

    @property
    def docr(self):
        return self.listfile_record.docr or self.cr_settings


class MSpecFileCollection:
    def __init__(self, files, stacked_path='', gzipped=False):
        self.files = files
        self.r = list(i for i in files if i.side == 'r')
        self.b = list(i for i in files if i.side == 'b')
        self.files.sort(key=lambda x: x.num)
        self.r.sort(key=lambda x: x.num)
        self.b.sort(key=lambda x: x.num)
        self.stacked_path = stacked_path
        self.gzipped = gzipped

    @property
    def sides(self):
        yield 'r', self.r
        yield 'b', self.b
        # raise StopIteration

    def __repr__(self):
        return f'<MSpecFileCollection: r={str([repr(x) for x in self.r])}, b={str([repr(x) for x in self.b])}>'

    def stacked_file(self, side, gzipped=None, stacked_path=''):
        if gzipped is False:
            gzipped = ''
        elif gzipped is True:
            gzipped = '.gz'
        else:
            gzipped = '.gz' if self.gzipped else ''

        if not stacked_path:
            stacked_path = self.stacked_path

        side=side.lower()
        x = self.r if side == 'r' else self.b
        if not x:
            return ''
        seqnos = [a.seqno for a in x]
        return os.path.join(stacked_path, f"{side}{rangify(seqnos, delim='_')}.fits{gzipped}")

    @property
    def flat_or_twilight(self):
        f = self.files[0]
        return f.led or f.quartz or f.twilight

    def __iter__(self):
        return iter(self.files)


def parse_listfile(listfile):
    ret = set()
    datestr = ''
    with open(listfile, 'r') as lf:
        for i, l in enumerate(lf):
            if not l.strip():
                continue
            if l.startswith('#'):
                try:
                    datestr = datetime.strptime(l.strip(), '#%Y%m%d').strftime('ut%Y%m%d')
                except ValueError:
                    pass
                continue
            try:
                ret.add(SeqStr(l, ut=datestr))
            except:
                print(f"Skipping l{i}: '{l}'")
    return list(ret)


def glob_for_fits(source_dir):
    basedirs = glob(source_dir) if '*' in source_dir else (source_dir,)
    return [os.path.join(dirpath, f)
            for basedir in basedirs for dirpath, dirnames, files in os.walk(basedir) for f in files
            if 'c1.fits' in f and f not in ('rc1.fits', 'bc1.fits') and not f.startswith('.')]


def generate_filesets(listfile, source_dir, merged_dir, stacked_dir, gzip_stack=True):
    seq_strs = parse_listfile(listfile)
    files = glob_for_fits(source_dir)
    msf = [MSpecFits(f, merged_path=merged_dir, load=False) for f in files]
    msfd = {x.id: x for x in msf}
    for ss in seq_strs:
        for id in ss.ids:
            msfd[id].listfile_record = ss
    return [MSpecFileCollection([msfd[id] for id in ss.ids], stacked_path=stacked_dir,
                                gzipped=gzip_stack)
            for ss in seq_strs], seq_strs


def find(seq, msf) -> List[MSpecFits]:
    mapping = {im.id: im for pair in msf for im in pair}
    return [mapping.get(id, None) for id in SeqStr(seq).ids]
