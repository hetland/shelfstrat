
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import netCDF4

import octant.tools

nc = netCDF4.Dataset('simulations/shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04/shelfstrat_his.nc')

yidx = 46  # 50 m isobath

dx = 1.0/nc.variables['pm'][1, 1]
h = nc.variables['h'][yidx:yidx+1, :].mean()

u = nc.variables['u'][:, -1, yidx:yidx+1, :].mean(axis=-2)
v = nc.variables['v'][:, -1, yidx, :]
v = 0.5*(v[:, 1:] + v[:, :-1])

salt = nc.variables['salt'][:, -1, yidx:yidx+1, :].mean(axis=-2)
temp = nc.variables['temp'][:, -1, yidx:yidx+1, :].mean(axis=-2)

t = nc.variables['ocean_time'][:]/86400.0

k = 0.5 * (u**2 + v**2)
k = k[:, :-1]  # trim duplicate column from periodic BC


f = 1e-4
N2 = 1e-4

Rd = np.sqrt(N2) * h / f

NWindow = 64
freq, Pxx = signal.welch(k, fs=2*np.pi/dx, nperseg=NWindow, noverlap=NWindow/2)



slope = []
for n in np.arange(len(t)):
    Pxx_n = Pxx[n:n+24, :].mean(axis=0)
    plt.loglog(freq*Rd, Pxx_n, color='k', alpha=0.2, lw=0.25)
    slope.append(np.polyfit(np.log10(freq[5:15]*Rd), np.log10(Pxx_n[5:15]), 1)[0])

for n in np.arange(len(t))[24:72:8]:
    Pxx_n = Pxx[n:n+24, :].mean(axis=0)
    plt.loglog(freq*Rd, Pxx_n, color='k', lw=0.25)

for n in np.arange(len(t))[[24, 72]]:
    Pxx_n = Pxx[n:n+24, :].mean(axis=0)
    plt.loglog(freq*Rd, Pxx_n, color='r')


k = np.logspace(-5, -2)
y2 = (k*Rd/3.0)**-1.666   # 5/3
y3 = (k*Rd/3.0)**-3

plt.plot(k*Rd, y3, '-r', lw=2)
plt.plot(k*Rd, y2, '-m', lw=2)

plt.show()