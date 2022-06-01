# System stuff
import os
from collections import OrderedDict
import vaex
import numpy as np
import astropy.coordinates as coord
import astropy.units as u
from scipy.interpolate import interp1d

## Local imports
from weave_streams.utils.ebv import get_SFD_dust, coef
from weave_streams.utils.util import mkpol, pc2mu, inside_poly, confLoad
from weave_streams.coords import gd1, orphan, pal5, tripsc

# Utility dict to link stream name to coordinate class.
# There is probably a better way of doing this
sclassdict = {'gd1': gd1.GD1, 'orphan': orphan.Orphan, 'pal5': pal5.Pal5, 'tripsc': tripsc.TriPsc}

from vaex.astro import transformations  # # Needed for vaex>=3.0
vaex.astro.transformations.comat["gd1"] = gd1.R
vaex.astro.transformations.comat["pal5"] = pal5.R
vaex.astro.transformations.comat["orphan"] = orphan.R
vaex.astro.transformations.comat["tripsc"] = tripsc.R

## Get configuration and data directories
dirpref = os.path.dirname(__file__)
datadir = dirpref+'/data/'
confdir = dirpref+'/conf/'

def get_data(sname, phi1min, phi1max, config='config.yaml'):
    """
    Load data from local dataframes in vaex. Could be any generic data generating function.

    Keep vaex file format instead of converting to table, this keeps the process efficient.
    """

    cnames = confLoad(f'{confdir}/columns_gaia_x_ps1.yaml')['cnames']
    cfg = confLoad(config)

    df = vaex.open(f"{cfg['gaiadir']}/gaia-edr3.hdf5")
    df.join(vaex.open(f"{cfg['ps1dir']}/panstarrs1.hdf5"), inplace=True)
    df.add_variable("pi", math.pi)

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

def makecat(sname, output='default', config='config.yaml'):

    sclass = sclassdict[sname]
    c = confLoad(f'{confdir}/{sname}.yaml')  # stream specific configuration
    cfg = confLoad(config)

    age = c["age"]
    feh = c["feh"]
    lage = np.round(np.log10(age * 1e9), 2)

    try:
        pmpoly = np.loadtxt(datadir+c["pmpoly"])
        haspmpoly = True
    except:
        haspmpoly = False
        print("No PM polygon")

    f = np.loadtxt(datadir+c["trackpos"][0], usecols=c["trackpos"][1], unpack=False)
    t = OrderedDict()

    for i, n in enumerate(c["trackpos"][2]):
        t[n] = f[:, i]

    if "trackpm" in c.keys():
        f = np.loadtxt(datadir+c["trackpm"][0], usecols=c["trackpm"][1], unpack=False)
        t2 = OrderedDict()
        for i, n in enumerate(c["trackpm"][2]):
            t2[n] = f[:, i]
    else:
        t2 = {}
        t2["phi1"] = t["phi1"]

    #Produce interpolated tracks
    #fphi2 = interp1d( t["phi1"], t["phi2"], bounds_error=False, fill_value=(t["phi2"][0], t["phi2"][-1]),)
    fphi2 = interp1d(t["phi1"], t["phi2"], bounds_error=False,
                     fill_value="extrapolate")

    # This is not used
    if "sphi2" in t.keys():
        # if only width is given
        fW = interp1d(t["phi1"], t["sphi2"], bounds_error=False,
                      fill_value=(t["sphi2"][0], t["sphi2"][-1]))
    elif "phi2up" in t.keys():
        # if up and lower bound are given
        fWup = interp1d(t["phi1"], t["phi2up"], bounds_error=False,
                        fill_value=(t["phi2up"][0], t["phi2up"][-1]),)
        fWdown = interp1d( t["phi1"], t["phi2down"], bounds_error=False,
                          fill_value=(t["phi2down"][0], t["phi2down"][-1]),)
    else:
        # if nothing is given, assume dphi2
        fW = interp1d(t["phi1"], np.zeros_like(t["phi1"]) + c["dphi2"],
                      bounds_error=False)

    if "pmphi1" in t2.keys():
        fvx = interp1d(t2["phi1"], t2["pmphi1"], bounds_error=False,
                       fill_value=(t2["pmphi1"][0], t2["pmphi1"][-1]),)
    else:
        fvx = interp1d(t2["phi1"], np.zeros_like(t2["phi1"]),
                       bounds_error=False)

    if "pmphi2" in t2.keys():
        fvy = interp1d(t2["phi1"], t2["pmphi2"], bounds_error=False,
                       fill_value=(t2["pmphi2"][0], t2["pmphi2"][-1]),)
    else:
        fvy = interp1d(t2["phi1"], np.zeros_like(t2["phi1"]),
                       bounds_error=False)

    # Get isochrone MIST both in Gaia and PS1 passband
    ifeh, iage, label, iG, iBP, iRP, phase = np.loadtxt(f"{cfg['isocdir']}/Gaia/MIST_-1.5to-2.1.cmd",
                                                        usecols=(7, 1, 0, 22, 23, 24, -1),
                                                        unpack=True)
    iifeh, iiage, llabel, ig, ir, ii, iz, pphase = np.loadtxt(f"{cfg['isocdir']}/PS1/MIST_-1.5to-2.1.cmd",
                                                              usecols=(7, 1, 0, 9, 10, 11, 12, -1),
                                                              unpack=True)

    j = (ifeh == feh) & (iage == lage) & (phase <= 2)
    jj = (iifeh == feh) & (iiage == lage) & (pphase <= 2)
    iiG = iG[j]
    iiRP = iRP[j]
    iiBP = iBP[j]
    iig = ig[jj]
    iir = ir[jj]
    iii = ii[jj]
    iiz = iz[jj]

    # Read catalogue
    df = get_data(sname, c["phi1range"][0], c["phi1range"][1], config=config)
    df.add_variable('pi', np.pi)

    if "dist" in t2.keys():
        cc = np.polyfit(t2["phi1"], t2["dist"], 2)
        df["distance"] = np.polyval(cc, df["nphi1"])
    else:
        df["distance"] = vaex.vrange(0, len(df))
        df["distance"] = (df["distance"] * 0 + c["dist"])
        print("No distance info, using fixed")

    radial_velocity = np.zeros(len(df))
    df['radial_velocity'] = radial_velocity

    coo = coord.SkyCoord( ra=df.ra.values * u.deg, dec=df.dec.values * u.deg,
                         pm_ra_cosdec=df.pmra.values * u.mas / u.yr,
                         pm_dec=df.pmdec.values * u.mas / u.yr,
                         distance=df.distance.values * u.kpc,
                         radial_velocity=radial_velocity * u.km / u.s,)

    v_sun = coord.Galactocentric().galcen_v_sun
    observed = coo.transform_to(coord.Galactic)
    rep = observed.cartesian.without_differentials()
    rep = rep.with_differentials(observed.cartesian.differentials["s"] + v_sun)
    c2 = coord.Galactic(rep).transform_to(sclass)

    df["rpmphi1"] = c2.pm_phi1_cosphi2.value
    df["rpmphi2"] = c2.pm_phi2.value

    c3 = coo.transform_to(sclass)
    df["pmphi1"] = c3.pm_phi1_cosphi2.value
    df["pmphi2"] = c3.pm_phi2.value
    df["phi1"] = c3.phi1.value
    df["phi2"] = c3.phi2.value

    ebv = get_SFD_dust(df.l.values, df.b.values, bdir=cfg['ebvmaps'])

    aV = 3.1 * np.array(ebv)
    aG = coef["G"] * aV
    aRP = coef["RP"] * aV
    aBP = coef["BP"] * aV
    ag = coef["g"] * aV
    ar = coef["r"] * aV
    ai = coef["i"] * aV
    az = coef["z"] * aV

    ## Apply extinction correction and normalize distance to mean distance
    fmu = pc2mu(df.distance.values * 1000)
    if "dist" in t2.keys():
        mdmod = np.mean(pc2mu(t2["dist"] * 1000))
    else:
        mdmod = pc2mu(c["dist"] * 1000)

    ddmod = fmu - mdmod

    ## Extinction corrected
    df["G0"] = df.phot_g_mean_mag.values - aG
    df["RP0"] = df.phot_rp_mean_mag.values - aRP
    df["g0"] = df.g_mean_psf_mag.values - ag
    df["r0"] = df.r_mean_psf_mag.values - ar
    df["i0"] = df.i_mean_psf_mag.values - ai
    df["z0"] = df.z_mean_psf_mag.values - az

    # Distance homogenized
    df["dG0"] = df.phot_g_mean_mag.values - ddmod - aG
    df["dRP0"] = df.phot_rp_mean_mag.values - ddmod - aRP
    df["dg0"] = df.g_mean_psf_mag.values - ddmod - ag
    df["dr0"] = df.r_mean_psf_mag.values - ddmod - ar
    df["di0"] = df.i_mean_psf_mag.values - ddmod - ai
    df["dz0"] = df.z_mean_psf_mag.values - ddmod - az

    # %%
    polG  = mkpol(mdmod, 0.5, iiG, iiRP, offset=c["GRPoff"], p=[0.005, 24.6, 1.5, 12.0, 0.9])
    polgi = mkpol(mdmod, 0.5, iig, iii, offset=c["gioff"],   p=[0.005, 24.6, 1.5, 13.5, 0.9])
    poliz = mkpol(mdmod, 0.5, iii, iiz, offset=c["gioff"],   p=[0.005, 24.6, 1.5, 13.5, 0.9])

    df.select_lasso("dG0 - dRP0", "dG0", polG[:, 0], polG[:, 1], name="inG")
    df.select((df.dG0 > c["HBG"][0]) * (df.dG0 < c["HBG"][1]) * (df.dG0 - df.dRP0 < 0.6),
              name="inHB")

    df.select_lasso("dg0 - di0", "dg0", polgi[:, 0], polgi[:, 1], name="ingi")
    df.select_lasso("di0 - dz0", "di0", poliz[:, 0], poliz[:, 1], name="iniz")

    df["inG"] = inside_poly(np.c_[df.dG0.values - df.dRP0.values, df.dG0.values], polG)
    df["ingi"] = inside_poly(np.c_[df.dg0.values - df.di0.values, df.dg0.values], polgi)
    df["iniz"] = inside_poly(np.c_[df.di0.values - df.dz0.values, df.di0.values], poliz)
    df["inHB"] = (df.dG0 > c["HBG"][0]) * (df.dG0 < c["HBG"][1]) * (df.dG0 - df.dRP0 < 0.6)

    dpm = c["dpm"] ## Delta PM to select around stream track

    df["fvx"] = fvx(df.phi1.values)
    df["fvy"] = fvy(df.phi1.values)

    if sname == "orphan":
        ## Orphan uses reflex corrected tracks
        df["invx"] = (df.rpmphi1 < df.fvx + dpm) * (df.rpmphi1 > df.fvx - dpm)
        df["invy"] = (df.rpmphi2 < df.fvy + dpm) * (df.rpmphi2 > df.fvy - dpm)
        df.select((df.rpmphi1 < df.fvx + dpm) * (df.rpmphi1 > df.fvx - dpm), name="invx")
        df.select((df.rpmphi2 < df.fvy + dpm) * (df.rpmphi2 > df.fvy - dpm), name="invy")
    else:
        df["invx"] = (df.pmphi1 < df.fvx + dpm) * (df.pmphi1 > df.fvx - dpm)
        df["invy"] = (df.pmphi2 < df.fvy + dpm) * (df.pmphi2 > df.fvy - dpm)
        df.select((df.pmphi1 < df.fvx + dpm) * (df.pmphi1 > df.fvx - dpm), name="invx")
        df.select((df.pmphi2 < df.fvy + dpm) * (df.pmphi2 > df.fvy - dpm), name="invy")

    if haspmpoly != False:
        print("Doing the polygon PM cut")
        df.select_lasso("rpmphi1", "rpmphi2", pmpoly[:, 0], pmpoly[:, 1], name="pmpoly")
        df["pmcut"] = inside_poly(np.c_[df.rpmphi1.values, df.rpmphi2.values], pmpoly)

    df["tphi2"] = df.phi2.values - fphi2(df.phi1.values)
    if output=='default':
        # default name
        df.export_hdf5("{0}_dataframe_edr3.hdf5".format(sname), progress=True)
        return df
    elif output:
        # user provided name
        df.export_hdf5(output, progress=True)
        return df
    else:
        # dont export
        return df