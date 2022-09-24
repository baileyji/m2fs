from skimage.feature import peak_local_max


def find_peaks(im, min_sep=10, minv=25, maxv=5000):
    """find points to use for psf measurements"""
    points = peak_local_max(im, min_distance=min_sep, threshold_abs=minv, threshold_rel=0)
    if maxv is None:
        return points.squeeze()
    if points.squeeze().ndim > 1:
        return points[im[points[:, 0], points[:, 1]] < maxv]
    else:
        return points[im[points] < maxv]
