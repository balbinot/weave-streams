#!/usr/bin/env python

import os
import vaex
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as mpl_path
import re
import yaml
import astropy.coordinates as C
import astropy.units as u
import sqlutilpy
from pathlib import Path




path_matcher = re.compile(r'\$\{([^}^{]+)\}')
def path_constructor(loader, node):
  ''' Extract the matched value, expand env variable, and replace the match '''
  value = node.value
  match = path_matcher.match(value)
  env_var = match.group()[2:-1]
  return os.environ.get(env_var) + value[match.end():]

def fileExists(path):
    my_file = Path(path)
    return my_file.is_file()

def confLoad(fname):
    yaml.add_implicit_resolver('!path', path_matcher)
    yaml.add_constructor('!path', path_constructor)
    with open(fname, 'r') as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)

    return cfg

def deduplicate_dataframe(df, columns=None):
    """
    Implementation of deduplication of a vaex dataframe based on a notebook shared by maartenbreddels
    in https://github.com/vaexio/vaex/issues/746

    :param df: a vaex.DataFrame
    :param columns: the column wrt which we want to deduplicate
    :return: the dataframe with duplicates removed
    """
    initial_columns = df.get_column_names()
    if columns is None:
        column_names = df.get_column_names()
    else:
        column_names = columns
    sets = [df._set(column_name) for column_name in column_names]
    counts = [set.count for set in sets]
    set_names = [df.add_variable(f'set_{column_name}', set, unique=True) for column_name, set in
                 zip(column_names, sets)]

    # create new column 'row_id' that gives each unique 'row' the same number/id
    expression = df[f'_ordinal_values({column_names[0]}, {set_names[0]})'].astype('int64')
    product_count = 1
    for count, set_name, column_name in zip(counts[:-1], set_names[1:], column_names[1:]):
        product_count *= count
        expression = expression + df[f'_ordinal_values({column_name}, {set_name})'].astype('int64') * product_count
    df['row_id'] = expression

    print(df.data_type(expression))

    # this is not 'stable', because it is multithreaded, we may get a different id each time
    index = df._index('row_id')
    unique_row_ids = df.row_id.unique()
    indices = index.map_index(unique_row_ids)

    deduped = df.take(indices)

    return deduped[initial_columns]

def _rename_wsdbcols(df):
    cols = ['gpsfmag', 'gpsfmagerr',
            'rpsfmag', 'rpsfmagerr',
            'ipsfmag', 'ipsfmagerr',
            'zpsfmag', 'zpsfmagerr']
    ncols = ['g_mean_psf_mag', 'g_mean_psf_mag_error',
             'r_mean_psf_mag', 'r_mean_psf_mag_error',
             'i_mean_psf_mag', 'i_mean_psf_mag_error',
             'z_mean_psf_mag', 'z_mean_psf_mag_error']

    for c, nc in zip(cols, ncols):
        df.rename(c, nc)
    return

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

def get_data(sname, phi1min, phi1max, config='config.yaml'):
    """
    Load data from local dataframes in vaex. Could be any generic data generating function.

    Keep vaex file format instead of converting to table, this keeps the process efficient.
    """

    cnames = confLoad(f'{confdir}/columns_gaia_x_ps1.yaml')['cnames']
    cfg = confLoad(config)

    df = vaex.open(f"{cfg['gaiadir']}/gaia-edr3.hdf5")
    df.join(vaex.open(f"{cfg['ps1dir']}/panstarrs1.hdf5"), inplace=True)
    df.add_variable("pi", np.pi)

    ## Do matrix rotation as virtual column for speed
    df.add_virtual_columns_celestial(long_in="ra", lat_in="dec", long_out="nphi1", lat_out="nphi2", _matrix=sname)

    if sname == "orphan":
        # Orphan is special from the start
        absphi2 = np.abs(df.nphi2) < 14
    else:
        absphi2 = np.abs(df.nphi2) < 7

    phi1sel =(df.nphi1 < phi1max) * (df.nphi1 > phi1min)

    ## Make df smaller and workable to add new columns with shorter length
    return df.to_copy(column_names=cnames, selection=phi1sel*absphi2)

def get_field_wsdb(ra, dec, output):


    wsdb = get_wsdb_host()
    wsdb_host = wsdb.split(':')[0]
    wsdb_user = wsdb.split(':')[3]

    wang = 180*u.deg
    gdr2col = "source_id ra dec parallax pmra pmdec phot_g_mean_mag phot_bp_mean_mag phot_rp_mean_mag ebv".split(' ')
    ps1col = "objid ra dec ebv gpsfmag gpsfmagerr rpsfmag rpsfmagerr ipsfmag ipsfmagerr zpsfmag zpsfmagerr".split(' ')
    cols2qry = ', '.join(['g.'+i for i in gdr2col] + ['ps.'+i for i in ps1col])

    ra = C.Angle(ra*u.deg).wrap_at(wang).value

    querystr = """
    select {cols} from
    panstarrs_dr1.stackobjectthin as ps
    FULL OUTER JOIN gaia_edr3_aux.panstarrs1bestneighbour as gps
     ON ps.objid = gps.original_ext_source_id
    FULL OUTER JOIN gaia_edr3.gaia_source as g
        ON g.source_id = gps.source_id
    WHERE
    q3c_radial_query(ps.ra, ps.dec, {rac:.4f}, {decc:.4f}, 1) AND
    (ps.ginfoflag3&panstarrs_dr1.detectionflags3('STACK_PRIMARY'))>0 AND
    (ps.rpsfmag-ps.rkronmag)<0.05 AND
    (ps.gpsfmag < 23);
    """.format(cols=cols2qry, rac=ra, decc=dec)

    data = sqlutilpy.get(querystr, db='wsdb', host=wsdb_host, user=wsdb_user)
    ddict = {}
    for k,c in enumerate(gdr2col + ps1col):
        ddict[c] = data[k]
    dftmp = vaex.from_dict(ddict)

    path = 'pointed_data/'
    os.makedirs(path, exist_ok=True)
    dftmp.export_hdf5(path+output, progress=True)

    return




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
