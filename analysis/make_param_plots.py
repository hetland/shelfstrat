
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from datetime import datetime, timedelta

import octant
import octant.roms
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='NetCDF filename to plot')
args = parser.parse_args()

xidx = 100
tidx = 14*8

figsize = (10, 6)

nc = netCDF4.Dataset(args.filename)

y = nc.variables['y_rho'][:, xidx] / 1000.0
h = nc.variables['h'][:, xidx]
zr = octant.roms.nc_depths(nc, grid='rho')[tidx][:, :, xidx]
temp = nc.variables['temp'][tidx, :, :, xidx]
salt = nc.variables['salt'][tidx, :, :, xidx]
rho = nc.variables['rho'][tidx, :, :, xidx]
u = nc.variables['u'][tidx, :, :, xidx]
t = octant.cf.time(nc, 'ocean_time')

days = np.floor(t[tidx] / 86400.0)
hours = np.floor( (t[tidx]-days*86400) / 24.0)
minutes = 0
timestr = 'day %02d %02d:%02d' % (days, hours, minutes)


fig_dens = plt.figure(figsize=figsize)
ax_dens = fig_dens.add_subplot(111)
crho = ax_dens.contour(y * np.ones_like(rho), zr, rho, np.linspace(rho.min(), rho.max(), 20), colors='k', linestyles='-')
plt.clabel(crho, fontsize=8, fmt='%3.2f', inline_spacing=1)
ax_dens.plot(y, -h, '-k', lw=2)
ax_dens.fill( np.hstack((y, y[-1], y[0])), np.hstack((-h, -h.max(), -h.max())), '0.8', lw=0)
umax = np.floor((100.0*np.abs(u).max()))/100.0
cf_u = ax_dens.contourf(y * np.ones_like(u), zr, u, np.linspace(-umax, umax, 13), extend='both', cmap=plt.cm.RdBu_r)
cax = plt.axes([0.17, 0.3, 0.35, 0.02])
plt.colorbar(cf_u, cax=cax, orientation='horizontal').set_label(r'Along-shore velocity [m s$^{-1}$]')
ax_dens.set_ylim(-h.max(), 0)
ax_dens.set_xlim(y[1], y[-2])
ax_dens.set_xlabel('Distance offshore [km]')
ax_dens.set_ylabel('Depth [m]')
ax_dens.set_title('Density - %s' % timestr)
plt.savefig('dens_%d.pdf' % tidx)



salt_cdict = {'red':  ((0.00, 0.4,  0.4),
                       (0.35, 0.3,  0.3),
                       (0.66, 1.0,  1.0),
                       (0.85, 0.9,  0.9),
                       (0.93, 0.75,  0.75),
                       (1.00, 0.83, 0.83)),
             'green': ((0.00,  0.4, 0.4),
                       (0.125, 0.3, 0.3),
                       (0.375, 1.0, 1.0),
                       (0.64,  1.0, 1.0),
                       (0.75,  0.5, 0.5),
                       (0.93,  0.5, 0.5),
                       (1.00,  0.8, 0.8)),
             'blue':  ((0.00, 0.7, 0.7),
                       (0.11, 1.0, 1.0),
                       (0.34, 1.0, 1.0),
                       (0.65, 0.0, 0.0),
                       (0.85,  0.6, 0.6),
                       (1.00, 0.8, 0.8))}

salt_cmap = plt.matplotlib.colors.LinearSegmentedColormap('salt_cmap', salt_cdict, 256)

xx = nc.variables['x_rho'][:, :130] / 1000.0
yy = nc.variables['y_rho'][:, :130] / 1000.0
sss = nc.variables['salt'][tidx, -1, :, :130]

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111)
# cf_salt = ax.contourf(xx, yy, sss, np.linspace(sss.min(), 35, 20), cmap=salt_cmap)
cf_salt = ax.contourf(xx, yy, sss, 20, cmap=salt_cmap)
# ax.contour(xx, yy, sss, np.arange(28, 35.25, 0.25), linewidths=0.2, colors='k')
 # ax.plot(xx[:, xidx], yy[:, xidx], '--k', lw=2)
# ax.set_title('Surface salinity - %s' % timestr)
ax.set_aspect(1.0)
plt.colorbar(cf_salt, ticks=[sss.min(), 35], format='%3.2f')
plt.savefig('sss_%d.pdf' % tidx)

plt.show()
