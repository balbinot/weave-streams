# System stuff
import os
import vaex
import numpy as np
import astropy.units as u
from astropy import coordinates as C
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d

## Local imports
from weave_streams.coords import gd1, orphan, pal5, tripsc
from weave_streams.utils.util import confLoad

# Utility dict to link stream name to coordinate class.
# There is probably a better way of doing this
sclassdict = {'gd1': gd1.GD1, 'orphan': orphan.Orphan, 'pal5': pal5.Pal5, 'tripsc': tripsc.TriPsc}

## Get configuration and data directories
dirpref = os.path.dirname(__file__)
datadir = dirpref+'/data/'
confdir = dirpref+'/conf/'

wang = 360*u.deg  # Wrapping angle for RA

def filter_data(sname, catalogue, config='config.yaml'):

    sclass = sclassdict[sname]
    c = confLoad(f'{confdir}/{sname}.yaml')  # stream specific configuration
    cfg = confLoad(config)                   # global configuration

    lab = c['label']

    tf = datadir+c['trackpos'][0]

    ## Not very smart, but do not want to change the YAML convention
    #cols = [(0,1), (0, 2), (0,1), (0,1)]
    cols = np.array(c['trackpos'][1])
    colnames = np.array(c['trackpos'][2])
    _cidx = np.where((colnames == 'phi1')|(colnames == 'phi2'))[0]
    print(colnames, _cidx, cols)
    cols=cols[_cidx]

    #obrange = [(-100, 40), (-30, 30), (60, 120), (-7, 7)]
    obrange = c['phi1range']

    #frame = [gd1.GD1, pal5.Pal5, orphan.Orphan, tripsc.TriPsc]
    #tdist = [2.5, 2, 3, 1]
    tdist = c['tdist']

    Gfaint  = cfg['maglim'][1]
    Gbright = cfg['maglim'][0]

    rfaint  = cfg['psmaglim'][1]
    rbright = cfg['psmaglim'][0]

    phi1, phi2 = np.loadtxt(tf, usecols=cols, unpack=True)
    jj = (phi1 < obrange[1])*(phi1 > obrange[0])
    phi1 = phi1[jj]
    phi2 = phi2[jj]
    fp2 = interp1d(phi1, phi2, bounds_error=False) ## interpolated track

    df = vaex.open(catalogue)

    # Generic cuts that apply to all streams
    incmd = ((df.ingi*df.inG)|df.inHB)
    inpos = (np.abs(df.tphi2) < tdist)
    insurvey = (df.r_mean_psf_mag < rfaint)*(df.r_mean_psf_mag > rbright)

    # Get position of streams-pointed fields
    rphi1, rphi2 = np.loadtxt(f"{datadir}/{c['pointedfields']}", unpack=True)

    ## Sect in PM, using single polygon cut OR along the track cut
    inpm = df[c['pmsel'][0]]
    if len(c['pmsel']) > 1:
        for _s in c['pmsel'][1:]:
            inpm *= df[_s]

    j = incmd * inpm * inpos * insurvey

    df['_pmcut'] = inpm.values

    print('Exporting filtered catalogue as hdf5')
    df[j].export_hdf5(catalogue.replace('.hdf5', '_filtered.hdf5'), progress=True)

    return

    #coo = C.SkyCoord(phi1, phi2, unit='deg', frame=sclass)
    #rcoo = C.SkyCoord(rphi1, rphi2, unit='deg', frame=sclass)
    #c2 = coo.transform_to(C.ICRS)
    #rc2 = rcoo.transform_to(C.ICRS)
    #ra, dec = c2.ra.wrap_at(wang).deg, c2.dec.deg
    #rra, rdec = rc2.ra.wrap_at(wang).deg, rc2.dec.deg

    #alltarg = 0
    ##ra, dec = np.loadtxt(fname, usecols=(1,2), unpack=True)
    #df = vaex.open(catalogue)
    #ra = df.ra.values[j.values]
    #dec = df.dec.values[j.values]

    #try:
    #    RA = np.r_[RA, ra]
    #    DEC = np.r_[DEC, dec]
    #except:
    #    RA = ra
    #    DEC = dec

    #coo = C.SkyCoord(ra, dec, unit='deg', frame=C.ICRS)
    #for x, y in zip(rra, rdec):
    #    ell = Ellipse(xy=[x,y],  width=2./np.cos(np.deg2rad(y)), height=2., facecolor='None', edgecolor='k', lw=2, zorder=100)
    #    p.gca().add_patch(ell)
    #    coocen = C.SkyCoord(x, y, unit='deg', frame=C.ICRS)
    #    sep = coo.separation(coocen)
    #    jjj = sep < 1*u.deg
    #    plot(ra[jjj], dec[jjj], 'k.', ms=1)
    #    print(lab, len(ra[jjj]))
    #    alltarg += len(ra[jjj])

    #print(lab, alltarg, 'all')
    #xlim(np.max(rra)+3, np.min(rra)-3)
    #ylim(np.min(rdec)-3, np.max(rdec)+3)

    #sra, sdec = rc2.ra.deg, rc2.dec.deg
    #np.savetxt(lab+'_WEAVE_tiles.dat', np.array([sra, sdec]).T, fmt=['%.5f', '%.5f'])
    ##p.plot(sra, sdec, 'o')


    return


def formatFits(filteredCatalogues, external=['sag', 'cetus', 'tripsc'], config='config.yaml'):
    """
    Get filtered hdf5 dataframes and format them to FITS, ready to submit as target lists.

    Add a list of external catalogues. TODO

    """

    cfg = confLoad(config)                   # global configuration

    df = vaex.open_many(filteredCatalogues)
    dfnew = df

    ## Extract data from DF
    ra = dfnew.ra.values
    dec = dfnew.dec.values
    source_id = dfnew.source_id.values
    #obj_id = dfnew.objid.values ## without underscore it is the WSDB version
    obj_id = dfnew.obj_id.values ## without underscore it is the WSDB version
    revid = np.zeros_like(ra).astype(np.int64) + 3
    Gmag = dfnew.phot_g_mean_mag.values
    rmag = dfnew.r_mean_psf_mag.values
    priority = np.zeros_like(ra) + 10

    for ename in external:
        # Concatenate with external catalogues
        f = parseExternalCatalogues(ename, cfg=cfg)
        #return (ra, dec, Gmag, source_id, rmag, revid, ps1id)

        ra = np.r_[ra, f[0]]
        dec = np.r_[dec, f[1]]
        Gmag = np.r_[Gmag, f[2]]
        source_id = np.r_[source_id, f[3]]
        rmag = np.r_[rmag, f[4]]
        revid = np.r_[revid, f[5]]
        obj_id = np.r_[obj_id, f[6]]
        priority = np.r_[priority, f[7]]


    ## Fix datatype
    ra = ra.astype(np.float64)
    dec =  dec.astype(np.float64)
    source_id = source_id.astype(np.int64)
    revid = revid.astype(np.int64)
    obj_id = obj_id.astype(np.int64)
    Gmag = Gmag.astype(np.float64)
    rmag =  rmag.astype(np.float64)

    print("all, unique = ", len(source_id), len(np.unique(source_id)))

    ## Find and remove them
    kkk = np.setdiff1d(np.arange(len(source_id)), np.unique(source_id, return_index=True)[1])
    source_id[kkk]
    fill = np.zeros_like(source_id)     #  array_length = 10
    fill[kkk] = 1                 #  indexes = [2,5,6]
    fill = fill.astype(np.bool)

    print('Duplicate IDs:', source_id[fill])
    ra        = ra[~fill]
    dec       = dec[~fill]
    source_id = source_id[~fill]
    revid     = revid[~fill]
    obj_id    = obj_id[~fill]
    rmag      = rmag[~fill]
    Gmag      = Gmag[~fill]
    priority  = priority[~fill]

    # Confirm removal
    print('all, unique =', len(source_id), len(np.unique(source_id)), 'shoud be equal now')

    # %%
    obj_id[obj_id==999999] = 0 ### This is needed if using Gaia crossmatch; WSDB sets it to zero

    # %%
    obj_id = obj_id.astype(np.int64)


#RA -- in degrees (double precision)
#DEC -- in degrees (double precision)
#SOURCE_ID -- gaia source_id (64 bit integer)
#GAIA_REV_ID -- gaia revision id ( integer number), It can be 2 for GDR2, 3 for GDR3 or 0 if there is no Gaia Source.
#PS1_ID -- PS1 id number if known/exists (0 otherwise) (64 bit integer)

    # %%
    tbl = Table()
    tbl.add_column(source_id, name='SOURCE_ID')
    tbl.add_column(obj_id, name='PS1_ID')
    tbl.add_column(revid, name='GAIA_REV_ID')
    tbl.add_column(ra, name='RA')
    tbl.add_column(dec, name='DEC')
#    tbl.add_column(Gmag, name='PHOT_G_MEAN_MAG')
#    tbl.add_column(rmag, name='R_MEAN_PSF_MAG')
    # %%

    catversion = '220622'

    ofile = f'STREAMSWIDE_{catversion}.fits'

    # %%
    tbl.write(ofile, overwrite=True)
    fits.setval(ofile, 'VERSION', value=f'{catversion}')
    fits.setval(ofile, 'CATALOG', value='STREAMSWIDE')

    ### TODO: save dataframe to use in pointed script
    #df3.export_hdf5(ofile.replace('.fits', '.DF.hdf5'))
    #dfnew.export_hdf5(ofile.replace('.fits', '.DFnew.hdf5'))

def parseExternalCatalogues(feature="sag", cfg=None):

    """
    feature: str (one of: sag, cetus)
    """

    c = confLoad(f'{datadir}/external/{feature}.yaml')  # stream specific configuration

    print(f"Parsing external catalog for {feature} using {c['filename']} as input")

    edr3 = vaex.open(f"{cfg['gaiadir']}/gaia-edr3.hdf5")
    edr3.join(vaex.open(f"{cfg['ps1dir']}/dr2_best_neighbour.hdf5"), inplace=True)
    edr3.join(vaex.open(f"{cfg['ps1dir']}/panstarrs1.hdf5"), inplace=True)

    if '.fits' in c['filename']:
        data = fits.open(f"{datadir}/{c['filename']}")[1]

        # Names of columns (required: RA, DEC, Gmag, source_id)
        _ra = c['cnames'][0]
        _dec = c['cnames'][1]
        _G =  c['cnames'][2]
        _sourceid = c['cnames'][3]

        ra        = data.data[_ra][0]
        dec       = data.data[_dec][0]
        Gmag      = data.data[_G][0]
        source_id = data.data[_sourceid][0]
        priority  = np.zeros_like(ra) + 10
    elif '.hdf5' in c['filename']:
        data = vaex.open(f"{datadir}/{c['filename']}")
        ra        = data.ra.values
        dec       = data.dec.values
        Gmag      = data.phot_g_mean_mag.values
        source_id = data.source_id.values
        priority  = data.priority.values

    Gfaint  = cfg['maglim'][1]
    Gbright = cfg['maglim'][0]
    rfaint  = cfg['psmaglim'][1]
    rbright = cfg['psmaglim'][0]
    decmin  = cfg['decmin']

    #jsag = (sagG < Gfaint)*(sagG > Gbright)*(sagdec>-10)
    #j = (rmag < rfaint)*(rmag > rbright)*(dec> decmin)
    #j = (Gmag < Gfaint)*(Gmag > Gbright)*(dec >  decmin)

    # Apply filter
    #ra        = ra[j]
    #dec       = dec[j]
    #Gmag      = Gmag[j]
    #source_id = source_id[j]
    revid     = np.zeros_like(ra).astype(np.int64) + np.int(c['gaiarev'])
    ps1id     = np.zeros_like(ra).astype(np.int64)

    dfext = vaex.from_arrays(ra=ra, dec=dec, G=Gmag, source_id=source_id.astype(np.int64), priority=priority)

    if c['gaiarev'] > 2:
        edr3joined = edr3.join(dfext, left_on='source_id', right_on='source_id', rsuffix='_ext')
        jj = ~edr3joined.source_id_ext.isna()
        dfext_fixid = edr3joined[jj]

    # If GDR2 do xmatch to update source_ids
    elif c['gaiarev'] <= 2:

        edr3joined = edr3.join(dfext, left_on='dr2_source_id', right_on='source_id', rsuffix='_ext')
        jj = ~edr3joined.source_id_ext.isna()
        dfext_fixid = edr3joined[jj]

        ## Check if sources with different source_id have similar magnitudes
        #newid = dfext_fixid.source_id.values
        #oldid = dfext_fixid.source_id_ext.values
        #newG = dfext_fixid.phot_g_mean_mag.values
        #oldG = dfext_fixid.G.values
        #ttt = newG - oldG
        #print(mean(ttt), std(ttt))

    # Get stuff again, but from official Gaia eDR3
    ra        = dfext_fixid.ra.values
    dec       = dfext_fixid.dec.values
    source_id = dfext_fixid.source_id.values ##This eDR3 surce_id; external catalogue has _ext suffix
    Gmag      = dfext_fixid.phot_g_mean_mag.values
    rmag      = dfext_fixid.r_mean_psf_mag.values
    priority  = dfext_fixid.priority.values
    j = (rmag < rfaint)*(rmag > rbright)*(dec> decmin)
    #j = (Gmag < Gfaint)*(Gmag > Gbright)*(dec >  decmin)

    # Apply filter (inplace)
    ra        = ra[j]
    dec       = dec[j]
    Gmag      = Gmag[j]
    source_id = source_id[j]
    rmag      = rmag[j]
    priority  = priority[j]
    revid = np.zeros_like(ra).astype(np.int64) + 3
    ps1id = np.zeros_like(ra).astype(np.int64)

    return (ra, dec, Gmag, source_id, rmag, revid, ps1id, priority)

