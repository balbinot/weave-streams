#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as mpl_path
import re
import yaml


path_matcher = re.compile(r'\$\{([^}^{]+)\}')
def path_constructor(loader, node):
  ''' Extract the matched value, expand env variable, and replace the match '''
  value = node.value
  match = path_matcher.match(value)
  env_var = match.group()[2:-1]
  return os.environ.get(env_var) + value[match.end():]

def confLoad(fname):
    yaml.add_implicit_resolver('!path', path_matcher)
    yaml.add_constructor('!path', path_constructor)
    with open(fname, 'r') as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)

    return cfg


def get_wsdb_host():
    home = os.environ['HOME']
    wsdb_file = home + '/.pgpass'
    if os.path.exists(wsdb_file):
        wsdb = open(wsdb_file, 'r').read()
    else:
        try:
            wsdb = os.environ['WSDB_HOST']
        except KeyError:
            print('''you must specify the host name of the database either
with the wsdb_host file in your $HOME directory or throguh a WSDB_HOST
environmental variable
            ''')
            raise
    return wsdb

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
