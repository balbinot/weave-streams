#/usr/bin/env python

import vaex
from astropy.io import fits
from astropy.table import Table
from weave_streams import crossmatcher
from weave_streams.utils.util import get_wsdb_host

wsdb = get_wsdb_host()
host = wsdb.split(':')[0]
user = wsdb.split(':')[3]


## Format Cetus
# missing source_id edr3

for f in ['data/Cetus.fits', 'data/TriPsc.fits']:

    print(f)
    tbl = Table.read(f)
    ra, dec = tbl['ra'], tbl['dec']
    source_id, gmag= crossmatcher.doit('gaia_edr3.gaia_source', ra,dec, 'source_id, phot_g_mean_mag', rad=2., host=host, db='wsdb', user=user)

    tbl['source_id'] = source_id
    tbl['ext_G_mag'] = gmag
    df = vaex.from_astropy_table(tbl)
    df.export_hdf5(f.replace('.fits', '_w_ids.hdf5'), progress=True)


