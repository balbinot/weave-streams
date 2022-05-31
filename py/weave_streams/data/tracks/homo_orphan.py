from pylab import *
from scipy.interpolate import interp1d

x, y = np.loadtxt('orphan_track.txt', unpack=True)
x1, pm1 = np.loadtxt('orphan_pm1.txt', unpack=True)
x2, dist = np.loadtxt('orphan_dist.txt', unpack=True)


fy = interp1d(x, y, bounds_error=False, kind='cubic')
fpm1 = interp1d(x1, pm1, bounds_error=False, kind='cubic')
fdist = interp1d(x2, dist, bounds_error=False, kind='quadratic')


X = np.arange(-100, 150, 1)

figure()
plot(X, fy(X))
plot(x, y, 'o')
ylabel('y')

figure()
plot(X, fpm1(X))
plot(x1, pm1, 'o')
ylabel('pm1')

figure()
plot(X, fdist(X))
plot(x2, dist, 'o')
ylabel('dist')



for xx in X:
    print(xx, fpm1(xx), fdist(xx))



show()
