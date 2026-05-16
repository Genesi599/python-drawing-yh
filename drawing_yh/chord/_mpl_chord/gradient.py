"""
Create linear color gradients
"""

from matplotlib.colors import ColorConverter, LinearSegmentedColormap
from scipy.ndimage import gaussian_filter

import numpy as np


def linear_gradient(cstart, cend, n=10):
    '''
    Return a gradient list of `n` colors going from `cstart` to `cend`.
    '''
    s = np.array(ColorConverter.to_rgb(cstart))
    f = np.array(ColorConverter.to_rgb(cend))

    rgb_list = [s + (t / (n - 1))*(f - s) for t in range(n)]

    return rgb_list


def gradient(start, end, min_angle, color1, color2, meshgrid, mask, ax,
             alpha):
    '''
    Create a linear gradient from `start` to `end`, which is translationally
    invarient in the orthogonal direction.
    The gradient is then cliped by the mask.
    '''
    xs, ys = start
    xe, ye = end

    X, Y = meshgrid

    # PATCHED: 用 distance ratio 替代 binary mask + gaussian filter,得真正 smooth gradient
    # (原版 (d2end<d2start) 是 binary 0/1,gaussian filter sigma 很小 → Z 几乎 binary,
    #  chord 大段相同色看着像两段。distance ratio = sqrt(d2start) / (sqrt(d2start)+sqrt(d2end))
    #  是平滑距离比,sender 端 Z=0、receiver 端 Z=1、中间真线性过渡。)
    d2start = (X - xs)*(X - xs) + (Y - ys)*(Y - ys)
    d2end   = (X - xe)*(X - xe) + (Y - ye)*(Y - ye)
    ds = np.sqrt(d2start)
    de = np.sqrt(d2end)
    Z = ds / (ds + de + 1e-12)

    # generate the colormap
    n_bin = 100

    color_list = linear_gradient(color1, color2, n_bin)

    cmap = LinearSegmentedColormap.from_list("gradient", color_list, N=n_bin)

    im = ax.imshow(Z, interpolation='bilinear', cmap=cmap,
                   origin='lower', extent=[-1, 1, -1, 1], alpha=alpha)

    im.set_clip_path(mask)
