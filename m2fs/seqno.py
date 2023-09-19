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
    def __init__(self, s: str, default_cr=False, ut='', run_id=''):
        """[r|b|]#-#,#... {docr|nocr|lamp|arc|led|quartz|bias|dark}
        images: i.e. r0003
        seq_ids: i.e. ut20220921/r0003 or r0003 (will not have a ut)
        """
        self.line = s
        words = s.strip().lower().split()
        seq_id = words[0]
        if words[0].startswith('ut'):
            ut, seq = words[0].split('/')
            words[0] = seq
            words.append(ut)

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
        self.seq_id = seq_id
        self.seq_str = seq  #does not include the ut
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
        if not self.seq_id.startswith('ut'):
            self.seq_id = os.path.join(ut or self.run_id,self.seq_id)
        if ut:
            assert self.ut == os.path.split(self.seq_id)[0]
        if not ut:
            assert self.run_id == os.path.split(self.seq_id)[0]

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

    def to_MSpecFits(self, raw_path='', merged_path=''):
        return [MSpecFits(os.path.join(self.ut, im+'c1.fits'), load=False, run_id=self.run_id,
                          raw_path=os.path.join(raw_path, self.ut) if 'ut' not in raw_path else raw_path,
                          merged_path=merged_path, listfile_record=self) for im in self.images]

    @property
    def sided(self):
        return self.side != 'rb'

    @property
    def seq_ids(self):
        return tuple(os.path.join(self.ut or self.run_id, f'{side}{x}') for side in self.side for x in self.nums)

    @property
    def seq_ids_r(self):
        if 'r' not in self.side:
            return tuple()
        return tuple(os.path.join(self.ut or self.run_id, f'r{x}') for x in self.nums)

    @property
    def seq_ids_b(self):
        if 'b' not in self.side:
            return tuple()
        return tuple(os.path.join(self.ut or self.run_id, f'b{x}') for x in self.nums)

    def __repr__(self):
        return f'<SeqStr: {self}>'

    def __str__(self):
        return self.seq_id


class MSpecCCDFileset:  # This is a fileset for the quadrants!!!
    """ a SIDED Set of mspec quadrant files since"""

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

        self.id = os.path.join(self.ut or self.run_id or '', f'{self.side}{self.num}')

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

    @property
    def seqno(self):
        return self.num

    @property
    def seq_id(self):
        return self.id

    def seqno_match(self, ss: (SeqStr, str)):
        if not isinstance(ss, SeqStr):
            ss = SeqStr(ss)
        return self.id in ss.seq_ids


class MSpecFits(MSpecCCDFileset):
    """A single r or b image"""
    def __init__(self, fileset, load=True, ut=None, run_id='', raw_path='', merged_path='',
                 listfile_record:SeqStr=None):
        if isinstance(fileset, MSpecCCDFileset):
            super(self).__init__(fileset.user_filename, ut=fileset.ut, run_id=fileset.run_id,
                             raw_path=raw_path or fileset.raw_path, merged_path=merged_path or fileset.merged_path)
        else:
            super(self).__init__(fileset, ut=ut, run_id=run_id, raw_path=raw_path, merged_path=merged_path)
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


class MSpecStack:
    def __init__(self, files: list[MSpecFits], stacked_path='', gzipped=False):
        self.files = files
        for f in self.files[1:]:
            assert f.side==self.files[0].side
        for f in self.files[1:]:
            assert f.run_id==self.files[0].run_id
        self.files.sort(key=lambda x: x.num)
        self.stacked_path = stacked_path
        self.gzipped = gzipped

    @property
    def side(self) -> str:
        return self.files[0].side

    @property
    def ut(self) -> str:
        return self.files[0].ut

    @property
    def seq_id(self) -> str:
        ss = rangify([a.seqno for a in self.files], delim='_')
        ut = self.files[0].ut
        run_id = self.files[0].run_id
        return os.path.join(ut or run_id, f'r{ss}')

    def file(self, gzipped=None, path_override='') -> str:
        if len(self.files)==1:
            return self.files[0].merged_file(gzipped=gzipped, merged_path=path_override)
        if gzipped is False:
            gzipped = ''
        elif gzipped is True:
            gzipped = '.gz'
        else:
            gzipped = '.gz' if self.gzipped else ''

        seqnos = [a.seqno for a in self.files]
        return os.path.join(path_override or self.stacked_path,
                            f"{self.side}{rangify(seqnos, delim='_')}.fits{gzipped}")

    @property
    def flat_or_twilight(self):
        f = self.files[0]
        return f.led or f.quartz or f.twilight

    def __iter__(self):
        return iter(self.files)

    @property
    def fits(self):
        return fits.open(self.file())


class MSpecFileCollection:
    """ This is a collection of MSpecFits for r and/or b sides

    The r and b sides are assumed to be on the same night/run but no checking is enforced

    The files are all assumed to be of the same type (e.g. lef, flat, etc as well)

    """
    def __init__(self, files: list[MSpecFits], stacked_path='', gzipped=False):
        bfiles = [i for i in files if i.side == 'b']
        rfiles = [i for i in files if i.side == 'r']
        self.r = None
        self.b = None
        if rfiles:
            self.r = MSpecStack(rfiles, stacked_path=stacked_path, gzipped=gzipped)

        if bfiles:
            self.b = MSpecStack(bfiles, stacked_path=stacked_path, gzipped=gzipped)
        if self.r is not None and self.b is not None:
            assert self.r.ut == self.b.ut
        self.stacked_path = stacked_path
        self.gzipped = gzipped

    @property
    def sides(self):
        yield 'r', self.r
        yield 'b', self.b

    def __repr__(self):
        return f'<MSpecFileCollection: r={str(self.r)}, b={str(self.b)}>'

    @property
    def seq_ids(self) -> tuple[str]:
        ret = []
        if self.r is not None:
            ret.append(self.r.seq_id)
        if self.b is not None:
            ret.append(self.b.seq_id)
        if ret[0].replace('r','b')==ret[1]:
            ret.append(ret[0].replace('r',''))
        return tuple(ret)

    # def stacked_file(self, side, gzipped=None, stacked_path='') -> str:
    #     if self
    #     if gzipped is False:
    #         gzipped = ''
    #     elif gzipped is True:
    #         gzipped = '.gz'
    #     else:
    #         gzipped = '.gz' if self.gzipped else ''
    #
    #     if not stacked_path:
    #         stacked_path = self.stacked_path
    #
    #     side=side.lower()
    #     x = self.r if side == 'r' else self.b
    #     if not x:
    #         return ''
    #     seqnos = [a.seqno for a in x]
    #     return os.path.join(stacked_path, f"{side}{rangify(seqnos, delim='_')}.fits{gzipped}")

    @property
    def flat_or_twilight(self):
        return (self.r or self.b).flat_or_twilight

    @property
    def files(self):
        return ([] if self.r is None else self.r.files) + ([] if self.b is None else self.b.files)

    # def __iter__(self):
    #     return iter(self.r.files)


def parse_listfile(listfile)->list[SeqStr]:
    ret = set()
    datestr = ''
    run_id = ''
    with open(listfile, 'r') as lf:
        for i, l in enumerate(lf):
            if not l.strip():
                continue
            if l.startswith('#'):
                if l.startswith('#rid'):
                    run_id = l[4:].strip().lower()
                else:
                    try:
                        datestr = datetime.strptime(l.strip(), '#ut%Y%m%d').strftime('ut%Y%m%d')
                    except ValueError:
                        pass
                continue
            try:
                ret.add(SeqStr(l, ut=datestr, run_id=run_id))
            except Exception as e:
                getLogger(__name__).warning(f"Skipping l{i}: '{l}' due to {e}")
    return list(ret)


def glob_for_fits(source_dir):
    basedirs = glob(source_dir) if '*' in source_dir else (source_dir,)
    return [os.path.join(dirpath, f)
            for basedir in basedirs for dirpath, dirnames, files in os.walk(basedir) for f in files
            if 'c1.fits' in f and f not in ('rc1.fits', 'bc1.fits') and not f.startswith('.')]


def generate_filesets(listfile, source_dir, merged_dir, stacked_dir, gzip_stack=True) -> tuple[List[MSpecFileCollection], List[SeqStr]]:
    seq_strs = parse_listfile(listfile)
    files = glob_for_fits(source_dir)
    msf = [MSpecFits(f, merged_path=merged_dir, load=False) for f in files]
    msfd = {x.id: x for x in msf}
    try:
        for ss in seq_strs:
            for x in ss.to_MSpecFits(source_dir, merged_dir):
                if x.id not in msfd:
                    getLogger(__name__).warning(f'{x.id} not found in {source_dir}. Run_id support issue?')
                    msfd[x.id]=x

            for id in ss.seq_ids:
                if msfd[id].listfile_record is None:
                    msfd[id].listfile_record = ss
                # try:
                #     msfd[id].listfile_record = ss
                # except KeyError:
                #     getLogger(__name__).warning(f'No matches found for {ss} in {source_dir}')
                #     for x in ss.to_MSpecFits(source_dir,merged_dir):
                #         if try:
                #             assert x.id not in msfd
                #         except AssertionError:
                #             raise
                #         msfd[x.id]=x

    except KeyError:
        raise

    return [MSpecFileCollection([msfd[id] for id in ss.seq_ids], stacked_path=stacked_dir,
                                gzipped=gzip_stack)
            for ss in seq_strs], seq_strs


class M2FSIFUMListfileDataset:
    def __init__(self, lfile, raw_dir, data_dir):
        merged_dir = os.path.join(data_dir,'merged')
        stacked_dir = os.path.join(data_dir, 'stacked')
        if not os.path.exists(merged_dir):
            os.makedirs(merged_dir)
        if not os.path.exists(stacked_dir):
            os.makedirs(stacked_dir)
        msf, seq_strs = generate_filesets(lfile, raw_dir, merged_dir=merged_dir, stacked_dir=stacked_dir)
        self.msf = msf  # A list of MSpecFileCollections
        self.seq_strs = seq_strs  # A list of SeqStr

    def __getitem__(self, item):
        """
        take a string, int or iterable of the same return a single or iterable
        strings are treated as seqstr and will match against MSpecFileCollections
        if r or b is specified it will be pulled from the collection
        initial 0s in strings are ignored

        e.g. 'r430'
        """
        if isinstance(item,(int, str)):
            single=True
            item = [item]
        ss = [SeqStr(i if isinstance(i, str) else f'{i:04}') for i in item]
        matches = []
        for s in ss:
            added = False
            for x in self.msf:
                if s.seq_id in x.seq_ids:
                    matches.append(x)
                    added=True
                    break
                individ_ss = [a.seq_id for a in x.files]
                try:
                    matches.append(x.files[individ_ss.index(s.seq_id)])
                    added = True
                    break
                except ValueError:
                    pass
            if not added:
                matches.append(None)

        for i,s in enumerate(ss):
            if s.sided and matches[i]:
                matches[i]=getattr(matches[i], s.side)

        if single:
            matches = matches[0] if matches else None

        return matches

def find(seq, msf) -> List[MSpecFits]:
    mapping = {im.id: im for pair in msf for im in pair}
    return [mapping.get(id, None) for id in SeqStr(seq).seq_ids]
