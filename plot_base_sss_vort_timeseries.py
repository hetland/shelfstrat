
import numpy as np
import matplotlib.pyplot as plt

import octant
import netCDF4
import cmocean

filename = 'simulations/shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04/shelfstrat_his.nc'
y_offset = 0.03


nc = netCDF4.Dataset(filename)

xr = nc.variables['x_rho'][:]/1000.
yr = nc.variables['y_rho'][:]/1000.

xp = nc.variables['x_psi'][:]/1000.
yp = nc.variables['y_psi'][:]/1000.

time = nc.variables['ocean_time'][:]/86400.

fig, axs = plt.subplots(4, 3, sharex=False, sharey=False, figsize=(9, 6))

vort_axs = axs[[0, 2]]
sss_axs = axs[[1, 3]]

days = [0, 2, 4, 6, 8, 10]
tidx = [0, 16, 32, 48, 64, 80]

for sax, vax, day, tidx in zip(sss_axs.flat, vort_axs.flat, days, tidx):
    sax.text(20, 100, 'day %d' % day)
    
    sss = nc.variables['salt'][tidx, -1]
    
    u = nc.variables['u'][tidx, -1]
    v = nc.variables['v'][tidx, -1]
    vort = (np.diff(v, axis=-1) - np.diff(u, axis=-2))/1000.0
    vort /= 1e-4
    
    cfs = sax.contourf(xr, yr, sss, [28, 29, 30, 31, 32, 33, 34, 35], cmap=cmocean.cm.salinity, extend='min')
    cfv = vax.contourf(xp, yp, vort, np.linspace(-3, 3, 13), cmap=cmocean.cm.vorticity, extend='max')
    
    vax.set_xticks(np.arange(0, 300, 50))
    vax.set_yticks(np.arange(0, 150, 50))
    
    sax.set_xticks(np.arange(0, 300, 50))
    sax.set_yticks(np.arange(0, 150, 50))
    
    if day!=6:
        sax.set_xticklabels([])
        sax.set_yticklabels([])
    
    vax.set_xticklabels([])
    vax.set_yticklabels([])
    
    sax.tick_params(labelsize=10)
    
    bbox = vax.get_position()
    bbox.y0 -= y_offset
    bbox.y1 -= y_offset
    vax.set_position(bbox)
    
    vax.set_aspect(1.0)
    sax.set_aspect(1.0)
    
    vax.set_xlim(0, 256)
    vax.set_ylim(0, 128)

    sax.set_xlim(0, 256)
    sax.set_ylim(0, 128)

caxv = fig.add_axes([0.915, 0.68, 0.015, 0.22])
cbv = plt.colorbar(cfv, cax=caxv)
cbv.set_ticks(np.linspace(-3, 3, 7))
cbv.set_label(r'$\zeta/f$ [s$^{-1}$]', fontsize=10)
caxv.tick_params(labelsize=10)

caxs = fig.add_axes([0.915, 0.08, 0.015, 0.22])
cbs = plt.colorbar(cfs, cax=caxs)
cbs.set_ticks(np.arange(28, 36))
cbs.set_label('Salinity [g kg$^{-1}$]', fontsize=10)
caxs.tick_params(labelsize=10)

sss_axs[1, 0].set_xlabel('Cross-shore distance [km]', fontsize=10)
sss_axs[1, 0].set_ylabel('Along-shore distance [km]', fontsize=10)

plt.show()
