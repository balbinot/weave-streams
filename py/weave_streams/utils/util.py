
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as mpl_path


def inside_poly(data, vertices):
    return mpl_path(vertices).contains_points(data)


def pc2mu(x):
    return 5 * np.log10(x) - 5


def mu2pc(x):
    return 10 ** ((x + 5) / 5.0)


def plotPreferences():
    font = {"family": "serif", "weight": "normal", "size": 16}
    plt.rcParams["axes.linewidth"] = 2  # set the value globally
    plt.rc("font", **font)
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'


def mkpol(mu, dmu, iblue, ired, offset=0, p=[0.001, 24.5, 1.1, 12.5, 0.9]):
    """ Builds ordered polygon for masking """
    err = lambda x: p[0] + np.exp((x - p[1]) / p[2]) + 1 / np.exp((x - p[3]) / p[4])
    c = iblue - ired + offset
    m = iblue + mu
    mnear = m - dmu / 2.0
    mfar = m + dmu / 2.0
    C = np.r_[c + 2 * err(mfar) + 0.05, c[::-1] - 2 * err(mnear[::-1]) - 0.05]
    M = np.r_[m + 0.5, m[::-1]]

    return np.c_[C, M]
