import numpy as np
import matplotlib.pyplot as plt


def imshow(im, vmaxp=93):
    plt.imshow(im, origin='lower', vmax=np.percentile(im, vmaxp))
    plt.colorbar()
