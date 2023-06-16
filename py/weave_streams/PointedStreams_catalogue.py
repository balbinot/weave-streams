#!/usr/bin/env python
# System stuff
import os

import numpy as np
import vaex
# Local imports
from weave_streams.utils.util import (angular_separation, confLoad,
                                      get_wsdb_host)

# Get configuration and data directories
dirpref = os.path.dirname(__file__)
datadir = dirpref + "/data/"
confdir = dirpref + "/conf/"

wsdb = get_wsdb_host()
wsdb_host = wsdb.split(":")[0]
wsdb_user = wsdb.split(":")[3]


def filterExport(sname, wfile, ps1file, config="config.yaml"):
    """
    INPUT:
        sname: stream name (gd1, pal5, orphan, or tripsc)
        wfile: filtered catalogue from the streamswide selection
        ps1file: file containing the outter join between Gaia eDR3/DR3 and PS1
                 (see WideStreams_catalogue wsdb data fetch)
        config: yaml file with global config
    """

    c = confLoad(f"{confdir}/{sname}.yaml")  # stream specific configuration
    cfg = confLoad(config)  # global configuration
    c["label"]

    ra, dec = np.loadtxt(f'{datadir}/fields/{c["label"]}_WEAVE_tiles.dat', unpack=True)

    dfw = vaex.open(wfile)
    dfps1 = vaex.open(ps1file)

    # Steps:
    # [x] Remove all Gaia stars from dfps1
    # [x] Apply ingi filter to dfps1
    # [x] Concatenate with dfw
    # [x] Remove duplicates
    # [x] Filter faint magnitude; set in conf/yaml config files.
    # [x] Export FITS; done via exportFits function on CLI

    # TODO
    # [ ] Check if priority is being properly exported
    # [ ]

    notGaia = dfps1.source_id == -9999
    df = dfps1[(notGaia) * dfps1.ingi]
    pddf = df.to_pandas_df()
    pddf = pddf.drop_duplicates(subset=["objid"])
    df = vaex.from_pandas(pddf)

    dfuni = dfw.concat(df)

    maglim = c["pointedfaint"]
    dmag = cfg["maglim"][1] - cfg["maglim"][0]

    magcut = (dfuni.g_mean_psf_mag < maglim) * (dfuni.g_mean_psf_mag > maglim - dmag)
    _magcut = magcut.values

    for rac, decc in zip(ra, dec):
        sep = angular_separation(rac, decc, dfuni.ra.values, dfuni.dec.values) < 1
        print(rac, decc, len(sep[sep * _magcut]))
        try:
            J = J | sep
        except:
            J = sep

    dfuni["inpointed"] = J

    expdf = dfuni[magcut * dfuni.inpointed].to_copy()

    expdf.export_hdf5(f"{sname}_dataframe_edr3_POINTED_filtered.hdf5", progress=True)

    return dfuni[magcut * dfuni.inpointed].to_copy()

    # for rac, decc in zip(ra, dec):
    #    cen = C.SkyCoord(ra=RA*u.deg, dec=DEC*u.deg)
    #    sep = cen.separation(coo)
    #    i = magcut*(sep < 1*u.deg)
    #    print(k, n, RA, DEC, len(coo[i]))
    #    try:
    #        II = II|i
    #    except:
    #        II = i


# %%


# %%
# dfuni['inpointed'] = II


if __name__ == "__main__":
    from sys import argv

    from pylab import *

    df = filterExport(argv[1], argv[2], argv[3])
    J = df.source_id.isna()
    df[~J].plot("ra", "dec", colormap="gray_r")
    df[~J].scatter("ra", "dec", length_check=False, s=1, c="r")

    print(df[~J].column_names)

    figure()
    df.scatter("g_mean_psf_mag - i_mean_psf_mag", "g_mean_psf_mag", length_check=False)
    ylim()[::-1]

    plt.show()
