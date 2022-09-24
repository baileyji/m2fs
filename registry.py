import atexit
import pickle

DATAFILE = '~/m2fsregistry.pickle'

save = True
_defaults = True
_reg = {'ifum.gaps': {'1x1': (10, 10), '2x1': (10, 10)},
        'ifum.psf': {'1x1': (10, 10), '2x1': (10, 10)}}


def update(k, v):
    global _reg, _defaults
    try:
        if _reg[k] == v:
            return
        else:
            raise KeyError
    except KeyError:
        _reg[k] = v
        _defaults = False


def persist():
    global _defaults
    if save and not _defaults:
        with open(DATAFILE, 'wb') as f:
            pickle.dump(_reg, f)
        _defaults = True


try:
    with open(DATAFILE) as f:
        _reg = pickle.load(f)
except IOError:
    pass


atexit.register(persist)
