from .info import M2FSObsInfo
from m2fs.obs.scatter import *


def info(file, no_load=False):
    return M2FSObsInfo(file, no_load=no_load)
