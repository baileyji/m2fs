import numpy as np

HR = dict(nfib=432, ngroup=8, pad=150, lower_pad=150 / 2, npixels=4096, group_pad=2,
          spatial_binning=2, gap_types=3)
# LSB
# 180/arm 10 groups
LSB = dict(nfib=180, ngroup=10, pad=150, lower_pad=150 / 2, npixels=4096, group_pad=0,
           spatial_binning=2, gap_types=2)
STD = dict(nfib=276, ngroup=6, pad=150, lower_pad=150 / 2, npixels=4096, group_pad=0,
           spatial_binning=2, gap_types=3)

M2FS_MRES = dict(nfib=128, ngroup=8, pad=150, lower_pad=150 / 2, npixels=4096, group_pad=0,
                 spatial_binning=2, gap_types=3)

M2FS_HRES_1ORD = dict(nfib=128, ngroup=8, pad=150, lower_pad=150 / 2, npixels=4096, group_pad=0,
                      spatial_binning=2, gap_types=3)


def nominal_centers_linear(nfib, ngroup, npixels, pad, lower_pad, group_pad, spatial_binning):
    nfib_group = nfib / ngroup
    npix_all = npixels - pad
    npix_group = npix_all / ngroup
    npix_fiber = (npix_group - group_pad) / nfib_group

    group_centers = (np.arange(ngroup) + 0.5) * npix_group  # No padding
    fiber_centers = (np.arange(nfib_group) + 0.5) * npix_fiber

    nominal_centers = (group_centers[:, None] + fiber_centers[None, :]).ravel()
    nominal_centers += lower_pad
    nominal_centers /= spatial_binning
    return nominal_centers


def mspec_fiber_spacing(mode):
    info = getattr(__name__, str.upper(mode))
    nominal_centers = nominal_centers_linear(**info)
    return nominal_centers, np.diff(nominal_centers)
