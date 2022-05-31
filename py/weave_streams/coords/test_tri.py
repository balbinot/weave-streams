#!/usr/bin/env python

from pylab import *
from astropy import coordinates as C
from astropy import units as u

from tripsc import TriPsc

ra, dec = np.loadtxt('../tracks/galstreams.footprint.Tri-Pis.dat',
                     usecols=(1,2), unpack=True)


coo = C.SkyCoord(ra, dec, unit='deg')
c2 = coo.transform_to(TriPsc)
phi1, phi2 = c2.phi1.deg, c2.phi2.deg
c3 = C.SkyCoord(phi1, phi2, unit='deg', frame=TriPsc)
c4 = coo.transform_to(C.ICRS)
ra2, dec2 = c4.ra.deg, c4.dec.deg

figure()
plot(phi1, phi2, 'o')

figure()
plot(ra, dec, 'o', zorder=99, ms=3)
plot(ra2, dec2, 'o', label='back', ms=10, zorder=98)
legend()


for x, y in zip(phi1,phi2):
    print(x,y)

show()
