
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from datetime import datetime, timedelta
import argparse
import os

import octant
import octant.roms

parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='NetCDF filename to plot')
parser.add_argument('--xidx', type=int, default=30, help='x-index slice to plot')
parser.add_argument('--tidx', type=int, default=14*8, help='x-index slice to plot')
args = parser.parse_args()

xidx = args.xidx
tidx = args.tidx

figsize = (6.5, 4)

nc = netCDF4.Dataset(args.filename)
casedir = os.path.split(os.path.split(args.filename)[0])[-1]
case = casedir.split('_')[1:]
keys = case[::2]
vals = [float(val) for val in case[1::2]]
params = dict(zip(keys, vals))
params['Ri'] = params['N2'] * params['f']**2 / params['M2']**2
params['phi'] = params['M2'] * 1e-3 * np.sqrt(params['Ri']) / params['f']**2

y = nc.variables['y_rho'][:, xidx] / 1000.0
h = nc.variables['h'][:, xidx]
zr = octant.roms.nc_depths(nc, grid='rho')[tidx][:, :, xidx]
temp = nc.variables['temp'][tidx, :, :, xidx]
salt = nc.variables['salt'][tidx, :, :, xidx]
rho = nc.variables['rho'][tidx, :, :, xidx]
u = nc.variables['u'][tidx, :, :, xidx]
t = nc.variables['ocean_time'][:]

days = np.floor(t[tidx] / 86400.0)
hours = np.floor( (t[tidx]-days*86400) / 24.0)
minutes = 0
timestr = 'day %02d %02d:%02d' % (days, hours, minutes)


####################################
# TEMPERATURE CONTOURS

plt.figure(figsize=figsize)
temp_ticks = np.arange(np.floor(temp.min()), np.ceil(temp.max()))
ctemp = plt.contour(y * np.ones_like(temp), zr, temp, temp_ticks, colors='k')
for cnt in temp_ticks:
    plt.contour(y * np.ones_like(temp), zr, temp, np.arange(cnt+0.25, cnt+1, 0.25), colors='0.8', linestyles='solid')

plt.clabel(ctemp, fontsize=8, fmt='%3.1f', inline_spacing=1)
plt.plot(y, -h, '-k', lw=2)
plt.fill( np.hstack((y, y[-1], y[0])), np.hstack((-h, -h.max(), -h.max())), '0.8', lw=0)
plt.ylim(-h.max(), 0)
plt.xlim(y[1], y[-2])
plt.xlabel('Distance offshore [km]')
plt.ylabel('Depth [m]')
plt.title('Temperature - %s' % timestr)
plt.savefig('temp_%s_%d.pdf' % (casedir, tidx))

####################################
# SALINITY CONTOURS

plt.figure(figsize=figsize)
salt_ticks = np.arange(np.floor(salt.min()), np.ceil(salt.max()))
csalt = plt.contour(y * np.ones_like(salt), zr, salt, salt_ticks, colors='k')
for cnt in salt_ticks:
    plt.contour(y * np.ones_like(salt), zr, salt, np.arange(cnt+0.25, cnt+1, 0.25), colors='0.8', linestyles='solid')

plt.clabel(csalt, fontsize=8, fmt='%3.1f', inline_spacing=1)
plt.plot(y, -h, '-k', lw=2)
plt.fill( np.hstack((y, y[-1], y[0])), np.hstack((-h, -h.max(), -h.max())), '0.8', lw=0)
plt.ylim(-h.max(), 0)
plt.xlim(y[1], y[-2])
plt.xlabel('Distance offshore [km]')
plt.ylabel('Depth [m]')
plt.title('Salinity - %s' % timestr)
plt.savefig('salt_%s_%d.pdf' % (casedir, tidx))

####################################
# DENSITY CONTOURS

fig_dens = plt.figure(figsize=figsize)
ax_dens = fig_dens.add_subplot(111)
dens_ticks = np.arange(np.floor(rho.min()), np.ceil(rho.max()))
crho = ax_dens.contour(y * np.ones_like(rho), zr, rho, dens_ticks, colors='k')
for cnt in dens_ticks:
    ax_dens.contour(y * np.ones_like(rho), zr, rho, np.arange(cnt+0.25, cnt+1, 0.25), colors='0.8', linestyles='solid')

plt.clabel(crho, fontsize=8, fmt='%3.1f', inline_spacing=1)
ax_dens.plot(y, -h, '-k', lw=2)
ax_dens.fill( np.hstack((y, y[-1], y[0])), np.hstack((-h, -h.max(), -h.max())), '0.8', lw=0)
if tidx == 0:
    u_ticks = np.linspace(0, 0.5, 11)
    cf_u = ax_dens.contourf(y * np.ones_like(u), zr, u, u_ticks, extend='both', cmap=plt.cm.Reds)
    cax = plt.axes([0.2, 0.3, 0.3, 0.02])
    ucb = plt.colorbar(cf_u, cax=cax, orientation='horizontal', ticks=[0, 0.5])
    ucb.ax.set_xticklabels([str(np.round(0, 2)), str(np.round(0.5, 2))])
    ucb.set_label(r'Along-shore velocity [m s$^{-1}$]')
else:
    umax = np.abs(u).max()
    u_ticks = np.linspace(-umax, umax, 21)
    cf_u = ax_dens.contourf(y * np.ones_like(u), zr, u, u_ticks, extend='both', cmap=plt.cm.RdBu_r)
    cax = plt.axes([0.2, 0.3, 0.3, 0.02])
    ucb = plt.colorbar(cf_u, cax=cax, orientation='horizontal', ticks=[-umax, 0, umax])
    ucb.ax.set_xticklabels([str(np.round(-umax, 2)), '0', str(np.round(umax, 2))])
    ucb.set_label(r'Along-shore velocity [m s$^{-1}$]')
    
ax_dens.set_ylim(-h.max(), 0)
ax_dens.set_xlim(y[1], y[-2])
ax_dens.set_xlabel('Distance offshore [km]')
ax_dens.set_ylabel('Depth [m]')
ax_dens.set_title('Density - %s' % timestr)
plt.savefig('dens_%s_%d.pdf' % (casedir, tidx))


xx = nc.variables['x_rho'][:] / 1000.0
yy = nc.variables['y_rho'][:] / 1000.0
ssd = nc.variables['rho'][tidx, -1]

####################################
# DENSITY CONTOURS

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)

dens_ticks = np.arange(np.floor(ssd.min()), np.ceil(ssd.max()))
crho_plan = ax.contour(xx, yy, ssd, dens_ticks, colors='k')
for cnt in dens_ticks:
    ax.contour(xx, yy, ssd, np.arange(cnt+0.25, cnt+1, 0.25), colors='0.8', linestyles='solid')

plt.clabel(crho_plan, fontsize=8, fmt='%3.1f', inline_spacing=1)

ax.set_title('Surface density - %s' % timestr)
ax.set_aspect(1.0)

# ax.text(0.05, 0.9,r'''N$^2$ = %4.2e 
# M$^2$ = %4.2e
# f = %4.2e''' % (params['N2'], params['M2'], params['f']),
#              horizontalalignment='left', verticalalignment='center',
#              transform = ax.transAxes)
# 
# ax.text(0.50, 0.9,r'''Ri = %4.2f
# $\phi$ = %4.2f''' % (params['Ri'], params['phi']),
#              horizontalalignment='left', verticalalignment='center',
#              transform = ax.transAxes)

# cb = plt.colorbar(cf_salt, ticks=[sss_ticks.min(), sss_ticks.max()])
# cb.ax.set_yticklabels([str(np.round(sss_ticks.min(), 2)), str(np.round(sss_ticks.max(), 2))])


plt.savefig('ssd_%s_%d.pdf' % (casedir, tidx))


plt.show()










