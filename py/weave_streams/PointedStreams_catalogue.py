#!/usr/bin/env python
# System stuff
import os
import vaex
import numpy as np
import astropy.coordinates as C
import astropy.units as u
import sqlutilpy

## Local imports
from weave_streams.utils.ebv import get_SFD_dust, coef
from weave_streams.utils.util import angular_separation, mkpol, pc2mu, inside_poly, confLoad, get_wsdb_host
from weave_streams.utils.util import get_wsdb_host
from weave_streams.coords import gd1, orphan, pal5, tripsc
from weave_streams.WideStreams_catalogue import makecat

## Get configuration and data directories
dirpref = os.path.dirname(__file__)
datadir = dirpref+'/data/'
confdir = dirpref+'/conf/'

wsdb = get_wsdb_host()
wsdb_host = wsdb.split(':')[0]
wsdb_user = wsdb.split(':')[3]


def filterExport(sname, wfile, ps1file, config='config.yaml'):

    c = confLoad(f'{confdir}/{sname}.yaml')  # stream specific configuration
    cfg = confLoad(config)                   # global configuration
    lab = c['label']

    ra, dec = np.loadtxt(f'{datadir}/fields/{c["label"]}_WEAVE_tiles.dat', unpack=True)

    dfw = vaex.open(wfile)
    dfps1 = vaex.open(ps1file)


    ## Steps:
    # [x] Remove all Gaia stars from dfps1
    # [x] Apply ingi filter to dfps1
    # [x] Concatenate with dfw
    # [x] Remove duplicates
    # [ ] Filter faint magnitude
    # [ ] Export FITS

    notGaia = dfps1.source_id == -9999
    df = dfps1[(notGaia)*dfps1.ingi]
    pddf = df.to_pandas_df()
    pddf = pddf.drop_duplicates(subset=['objid'])
    df = vaex.from_pandas(pddf)

    dfuni = dfw.concat(df)

    maglim = c['pointedfaint']
    dmag = cfg['maglim'][1] - cfg['maglim'][0]

    magcut = (dfuni.g_mean_psf_mag < maglim)*(dfuni.g_mean_psf_mag > maglim - dmag)
    _magcut = magcut.values

    for rac, decc in zip(ra, dec):
        sep = angular_separation(rac, decc, dfuni.ra.values, dfuni.dec.values) < 1
        print(rac, decc, len(sep[sep*_magcut]))
        try:
            J = J|sep
        except:
            J = sep

    dfuni['inpointed'] = J
    return dfuni[magcut*dfuni.inpointed].to_copy()

    #for rac, decc in zip(ra, dec):
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
    #dfuni['inpointed'] = II



if __name__=='__main__':
    from sys import argv
    from pylab import *
    df = filterExport(argv[1], argv[2], argv[3])
    J = df.source_id.isna()
    df[~J].plot('ra', 'dec', colormap='gray_r')
    df[~J].scatter('ra', 'dec', length_check=False, s=1, c='r')
    plt.show()



else:

    # %%
    config = """
    streams:
        gd1:
            trackpos:
                - tracks/track.txt
                - [0, 1, 3]
                - [phi1, phi2, sphi2]
            trackpm:
                - 'tracks/pms.txt'
                - [0, 3, 4, 5, 7] #file, colums(phi1, distkpc, edistkpc, pmphi1, pmphi2)
                - [phi1, dist, edist, pmphi1, pmphi2]
            pmpoly: 'tracks/gd1_pmpoly.dat'
            age: 12
            feh: -2.1
            gioff: 0.09
            GRPoff: 0.09
            phi1range: [-100, 40]
            dpm: 1.0
            HBG: [15, 16]

        pal5:
            trackpos:
                - tracks/pal5_track_final.txt
                - [0, 1, 2, 3]
                - [phi1, phi2down, phi2, phi2up]
            trackpm:
                - tracks/stream_obs_pal5_10.txt
                - [0, 1, 2, 4, 5]
                - [phi1, phi2, dist, pmphi1, pmphi2]
            pmpoly: None
            age: 10
            feh: -2.0
            gioff: 0.09
            GRPoff: 0.09
            phi1range: [-30, 30]
            dpm: 1.0
            HBG: [16.5, 17.5]

        tripsc:
            trackpos:
                - tracks/tripsc_track.dat
                - [0, 1]
                - [phi1, phi2]
            pmpoly: 'tracks/tripsc_pmpoly.dat'
            age: 12
            feh: -2.0
            dist: 26
            gioff: 0.09
            GRPoff: 0.09
            phi1range: [-8, 8]
            dpm: 1.0
            dphi2: 1.0
            HBG: [16.5, 17.5]
        orphan:
            trackpos:
                - tracks/orphan_homo_phi2.dat
                - [0, 1]
                - [phi1, phi2]
            trackpm:
                - tracks/orphan_homo_pm1_dist.dat
                - [0, 1, 2]
                - [phi1, pmphi1, dist]
            pmpoly: None
            age: 12
            feh: -2.1
            gioff: 0.09
            GRPoff: 0.09
            phi1range: [20, 140]
            dpm: 1.0
            dphi2: 1
            HBG: [17., 18.]


    """

    # %%
    # Cardelli 89 for Gaia + Pan-STARRS photometric system
    coef = {'G': 0.8620,
            'BP': 1.07198,
            'RP': 0.6465,
            'g': 1.16529,
            'r': 0.86813,
            'i': 0.67659,
            'z': 0.51743}

    # %%


    # %%
    #sname = 'gd1'
    #sclass = gd1.GD1
    #fl = glob('pointed_data/GD*.hdf5')
    #df = vaex.open_many(fl)

    #sname = 'pal5'
    #sclass = pal5.Pal5
    #fl = glob('pointed_data/Palomar5*.hdf5')
    #df = vaex.open_many(fl)

    sname = 'orphan'
    sclass = orphan.Orphan
    fl = glob('pointed_data/Orphan*.hdf5')
    df = vaex.open_many(fl)

    #sname = 'tripsc'
    #sclass = tripsc.TriPsc
    #fl = glob('pointed_data/Tri-Pis*.hdf5')
    #df = vaex.open_many(fl)


    cfg = safe_load(config)
    c = cfg['streams'][sname]

    age =   c['age']
    feh =   c['feh']
    lage = np.round(np.log10(age*1e9),2)

    try:
        pmpoly = np.loadtxt(c['pmpoly'])
        haspmpoly = True
    except:
        haspmpoly = False
        print("No PM polygon")

    f = np.loadtxt(c['trackpos'][0], usecols=c['trackpos'][1], unpack=False)
    t = OrderedDict()

    for i, n in enumerate(c['trackpos'][2]):
        t[n] = f[:,i]

    if 'trackpm' in c.keys():
        f = np.loadtxt(c['trackpm'][0], usecols=c['trackpm'][1], unpack=False)
        t2 = OrderedDict()
        for i, n in enumerate(c['trackpm'][2]):
            t2[n] = f[:, i]
    else:
        t2 = {}
        t2['phi1'] = t['phi1']

    ## Produce interpolated tracks
    fphi2 = interp1d(t['phi1'], t['phi2'], bounds_error=False, fill_value=(t['phi2'][0], t['phi2'][-1]), )
    fphi2 = interp1d(t['phi1'], t['phi2'], bounds_error=False, fill_value='extrapolate')
    if 'sphi2' in t.keys():
        # if only width is given
        fW =     interp1d(t['phi1'], t['sphi2'], bounds_error=False, fill_value=(t['sphi2'][0], t['sphi2'][-1]))
    elif 'phi2up' in t.keys():
        # if up and lower bound are given
        fWup    = interp1d(t['phi1'], t['phi2up'], bounds_error=False, fill_value=(t['phi2up'][0], t['phi2up'][-1]))
        fWdown  = interp1d(t['phi1'], t['phi2down'], bounds_error=False, fill_value=(t['phi2down'][0], t['phi2down'][-1]))
    else:
        # if nothing is given, assume dphi2
        fW =     interp1d(t['phi1'], np.zeros_like(t['phi1']) + c['dphi2'], bounds_error=False)

    if 'pmphi1' in t2.keys():
        fvx = interp1d(t2['phi1'], t2['pmphi1'], bounds_error=False, fill_value=(t2['pmphi1'][0], t2['pmphi1'][-1]))
    else:
        fvx = interp1d(t2['phi1'], np.zeros_like(t2['phi1']), bounds_error=False)

    if 'pmphi2' in t2.keys():
        fvy = interp1d(t2['phi1'], t2['pmphi2'], bounds_error=False, fill_value=(t2['pmphi2'][0], t2['pmphi2'][-1]))
    else:
        fvy = interp1d(t2['phi1'], np.zeros_like(t2['phi1']), bounds_error=False)

    XX = np.arange(c['phi1range'][0], c['phi1range'][1], 0.01)

    # Get isochrone PADOVA; deprecated
    #ifeh, iage, label, iG, iBP, iRP = np.loadtxt('isoc_gaia.dat', usecols=(1,2,9, 25,26,28), unpack=True)
    #j = (ifeh==feh)&(iage==lage)&(label<4)
    #iifeh, iiage, llabel, ig, ir, ii, iz = np.loadtxt('isoc_ps1.dat', usecols=(1,2,9, 25,26,27,28), unpack=True)
    #jj = (iifeh==feh)&(iiage==lage)&(llabel<4)

    # Get isochrone MIST both in Gaia and PS1 passband
    ifeh, iage, label, iG, iBP, iRP, phase = np.loadtxt('isocs/Gaia/MIST_-1.5to-2.1.cmd', usecols=(7, 1, 0, 22, 23, 24, -1),
               unpack=True)
    j = (ifeh==feh)&(iage==lage)&(phase<=2)

    iifeh, iiage, llabel, ig, ir, ii, iz, pphase = np.loadtxt('isocs/PS1/MIST_-1.5to-2.1.cmd', usecols=(7, 1, 0, 9, 10, 11, 12,
                                                         -1), unpack=True)
    jj = (iifeh==feh)&(iiage==lage)&(pphase<=2)

    iiG  = iG[j]
    iiRP = iRP[j]
    iiBP = iBP[j]
    iig  = ig[jj]
    iir  = ir[jj]
    iii  = ii[jj]
    iiz  = iz[jj]

    # %%


    # %%
    df.plot('nphi1', 'nphi2', f='log1p')

    # %%
    ## Reflex correction block. Do not use for WEAVE.
    # Compute reflex correction using interpolated distance at each position
    if 'dist' in t2.keys():
        cc = np.polyfit(t2['phi1'], t2['dist'], 2)
        df['distance'] = np.polyval(cc, df['nphi1'])
    else:
        df['distance'] = vaex.vrange(0, len(df))  # add a 'virtual range' column, which takes no memory
        df['distance'] = df['distance']* 0 + c['dist']  # multiply by zero, and add the initial value
        print("No distance info, using fixed")

    ## This may be needed
    #distances[distances<0] = 100*u.kpc

    # %%
    ebv       = df.ebv.values

    # %%
    aV        = 3.1        * np.array(ebv)
    aG        = coef['G']  * aV
    aRP       = coef['RP'] * aV
    aBP       = coef['BP'] * aV
    ag        = coef['g']  * aV
    ar        = coef['r']  * aV
    ai        = coef['i']  * aV
    az        = coef['z']  * aV

    # %%
    ## Apply extinction correction and normalize distance to mean distance
    fmu   = pc2mu(df.distance.values*1000)
    if 'dist' in t2.keys():
        mdmod = np.mean(pc2mu(t2['dist']*1000)) ## Using the original track, as this does not bias the distance to places with more stars
    else:
        mdmod = pc2mu(c['dist']*1000)

    # %%
    ## This brings the stream to the same average distance
    ddmod = fmu - mdmod

    # %%
    df

    # %%
    ddmod = ddmod.astype(np.float32)
    df.execute()

    df['aG'] = aG
    df['aRP'] = aRP
    df['ag'] = ag
    df['ar'] = ar
    df['ai'] = ai

    df['ddmod'] = ddmod

    ## Extinction corrected
    df['G0'] = df.phot_g_mean_mag - df.aG
    df['RP0'] = df.phot_rp_mean_mag - df.aRP
    df['g0'] = df.gpsfmag - df.ag
    df['r0'] = df.rpsfmag - df.ar
    df['i0'] = df.ipsfmag - df.ai

    # Distance homogenized
    df['dG0'] = df.phot_g_mean_mag    -    df.ddmod - df.aG
    df['dRP0'] = df.phot_rp_mean_mag  -    df.ddmod - df.aRP
    df['dg0'] = df.gpsfmag            -    df.ddmod - df.ag
    df['dr0'] = df.rpsfmag            -    df.ddmod - df.ar
    df['di0'] = df.ipsfmag            -    df.ddmod - df.ai

    # %%
    df.dg0.values

    # %%
    polG  = mkpol(mdmod, 0.5, iiG, iiRP, offset = c['GRPoff'], p = [0.005, 24.6, 1.5, 12.0, 0.9])
    polgi = mkpol(mdmod, 0.5, iig, iii,  offset = c['gioff'],  p = [0.005, 24.6, 1.5, 13.5, 0.9])
    poliz = mkpol(mdmod, 0.5, iii, iiz,  offset = c['gioff'],  p = [0.005, 24.6, 1.5, 13.5, 0.9])

    df['ingi'] = inside_poly(np.c_[df.dg0.values - df.di0.values , df.dg0.values], polgi)
    df.select_lasso('dg0 - di0', 'dg0', polgi[:,0], polgi[:,1], name='ingi')

    # %%


    # %%


    # %%
    dpm = c['dpm']
    df['fvx'] = fvx(df.phi1.values)
    df['fvy'] = fvy(df.phi1.values)

    # %%
    if sname=='gd1':
        df.select('(ingi)', name='inCMD')
    elif sname=='pal5':
        df.select('(ingi)', name='inCMD')
    elif sname=='orphan':
        df.select('(ingi)', name='inCMD')
    elif sname=='tripsc':
        df.select('(ingi)', name='inCMD')


    df.select('(source_id < 0)|(phot_g_mean_mag>20.5)', name='notgaia')
    notgaia = (df.source_id.values < 0)|(df.phot_g_mean_mag > 20.5)

    # %%
    df.scatter('phi1', 'phi2', selection='inCMD*notgaia', s=1, c='k', length_check=False, alpha=0.1)

    # %%
    df.plot('dg0 - di0', 'dg0', selection='inCMD*notgaia', limits=[[-2, 3], [23, 13]])

    # %%
    df.export_hdf5(sname+'_PS1only.hdf5', selection='inCMD*notgaia')

    # %%


    # %%
    fig = plt.figure(figsize=(16,10))
    ax1 = plt.subplot2grid((22, 4), (0, 0),  colspan=3, rowspan=7)
    ax2 = plt.subplot2grid((22, 4), (7, 0),  colspan=3, rowspan=7, sharex=ax1)
    ax3 = plt.subplot2grid((22, 4), (14, 0), colspan=3, rowspan=7, sharex=ax1)
    ax4 = plt.subplot2grid((22, 4), (0, 3), rowspan=6)
    ax5 = plt.subplot2grid((22, 4), (7, 3), rowspan=6)


    statsstr="""Statistics"""
    plt.text(0.75, 0.33, statsstr, ha='left', va='top',
              transform=fig.transFigure)
    plt.subplots_adjust(hspace=1.0,  left=0.065, bottom=0.05, wspace=0.3, top=0.945)

    sca(ax1)
    df.scatter('phi1', 'phi2', selection='pmcut*inCMD*(abs(tphi2)<2)', alpha=0.5, c='k', s=8, lw=0)
    plot(XX, fphi2(XX), 'C0-', lw=4, alpha=0.8)

    fov = plt.Circle((c['phi1range'][0]+10, 4), 1, color='C3', fill=False, lw=2)
    gca().add_artist(fov)
    gca().set_ylim(-10,10)


    sca(ax2)
    if sname=='orphan':
        df.scatter('phi1', 'rpmphi1', selection='pmcut*inCMD*intrack', alpha=0.5, c='C0', s=3, lw=0)
        df.scatter('phi1', 'rpmphi2', selection='pmcut*inCMD*intrack', alpha=0.5, c='C1', s=3, lw=0)
        df.scatter('phi1', 'rpmphi1', selection='inrossi', alpha=0.5, c='r', s=8, lw=0)
        df.scatter('phi1', 'rpmphi2', selection='inrossi', alpha=0.5, c='g', s=8, lw=0)
        ylim(-3, 6)
    else:
        df.scatter('phi1', 'pmphi1', selection='pmcut*inCMD*intrack', alpha=0.5, c='C0', s=3, lw=0)
        df.scatter('phi1', 'pmphi2', selection='pmcut*inCMD*intrack', alpha=0.5, c='C1', s=3, lw=0)
        df.scatter('phi1', 'pmphi1', selection='inrossi', alpha=0.5, c='r', s=8, lw=0)
        df.scatter('phi1', 'pmphi2', selection='inrossi', alpha=0.5, c='g', s=8, lw=0)
    gca().plot(XX, fvx(XX), 'k-', lw=2, zorder=99)
    gca().plot(XX, fvy(XX), 'k-', lw=2, zorder=99)
    gca().plot(XX, fvx(XX)+dpm, 'C0--', zorder=99)
    gca().plot(XX, fvx(XX)-dpm, 'C0--', zorder=99)
    gca().plot(XX, fvy(XX)+dpm, 'C1--', zorder=99)
    gca().plot(XX, fvy(XX)-dpm, 'C1--', zorder=99)


    dphi = c['phi1range'][1] - c['phi1range'][0]
    sca(ax3)
    df.plot1d(df.phi1, shape=int(dphi), limits=[c['phi1range'][0], c['phi1range'][1]], selection='pmcut*inCMD*intrack')



    sca(ax4)
    df.scatter(df.dg0 - df.di0, df.dg0, c='k', s=3, alpha=0.2, selection='pmcut*intrack*inCMD', lw=0)
    df.scatter(df.dg0 - df.di0, df.dg0, c='r', s=8, alpha=1, selection="rossi", lw=0)
    gca().plot(polgi[:,0], polgi[:,1], 'w-', lw=3)
    gca().plot(polgi[:,0], polgi[:,1], 'k-', lw=1)
    gca().set_xlim(-0.5, 2)
    gca().set_ylim(21, 14)

    sca(ax5)
    df.scatter(df.dG0 - df.dRP0, df.dG0, c='k', s=3, alpha=0.2, selection='pmcut*intrack*inCMD', lw=0)
    df.scatter(df.dG0 - df.dRP0, df.dG0, c='r', s=8, alpha=1, selection='rossi', lw=0)
    polG  = mkpol(mdmod, 0.5, iiG, iiRP, offset = c['GRPoff'], p = [0.005, 24.6, 1.5, 12., 0.9])
    gca().plot(polG[:,0], polG[:,1], 'w-', lw=3)
    gca().plot(polG[:,0], polG[:,1], 'k-', lw=1)
    gca().set_xlim(-0.5, 2)
    gca().set_ylim(21, 14)



    ax3.set_xlabel(r'$\phi_1$ [deg]')
    ax1.set_ylabel(r'$\phi_1$ [deg]')
    ax2.set_ylabel(r'$\mu_{\phi_1}$ or $\mu_{\phi_2}$ [mas/yr]')
    ax3.set_ylabel(r'$\Sigma$ [not deg$^{-2}$]')
    ax4.set_ylabel(r'$g_0$')
    ax4.set_xlabel(r'$(g-i)_0$')
    ax5.set_ylabel(r'$G_0$')
    ax5.set_xlabel(r'$(G-RP)_0$')
    ax1.tick_params(labelbottom=False)
    ax2.tick_params(labelbottom=False)
    ax1.set_xlabel('')
    ax2.set_xlabel('')
    tight_layout()

    # %%
    df.export_hdf5('{0}_dataframe.hdf5'.format(sname))

    # %%
    df2 = vaex.open('{0}_dataframe.hdf5'.format(sname))

    # %%
    df2.plot('phi1', 'phi2', selection='((ingi*inG))*invx*invy', f='log1p', shape=(128//2, 128//2), smooth_pre=0.5,  vmin=1.5, vmax=3)
    #df.scatter('phi1', 'phi2', selection='(inHB|(ingi*inG))*invx*(inrossi)', s=3, c='w', alpha=0.2)
    #df.scatter('phi1', 'phi2', selection='inrossi', s=3, c='k', alpha=0.2)

    # %%
