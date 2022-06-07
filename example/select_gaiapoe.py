## This code is not supposed to work.
## This is the selection of MS stars for Sergey.

SelCongifg.parallax_zeropoint = 17./1000. #in mas

def gaiapoe_sel(**kw):
    parallax = kw['parallax']
    POE = parallax/kw['parallax_error']
    dist = 1/(parallax + SelConfig.parallax_zeropoint)
    # Absolute magnitude needs extinction in the G-band "AG"
    MG = kw['phot_g_mean_mag'] - 5*np.log10(1000*dist) - 5 - kw['AG']
    Vtan = dist * np.sqrt(kw['pmra']**2 + kw['pmdec']**2) * 4.74047
    return (dist < 3.0) & (MG > 3.2) & (POE >= 5.0) & betw(Vtan, 200, 800)
