
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import octant
import octant.roms

ncfile = 'simulations/shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04/shelfstrat_his.nc'

nc = netCDF4.Dataset(ncfile)
time = nc.variables['ocean_time'][:]/86400.0

delta = 0.1

HMF = np.zeros((len(time), 129), 'd')
VMF = np.zeros((len(time), 129), 'd')
HBF = np.zeros((len(time), 129), 'd')
VBF = np.zeros((len(time), 130), 'd')

for tidx in range(len(time)):
    zw = octant.roms.nc_depths(nc, 'w')[tidx]
    dz = np.diff(zw, axis=0)
    pm = nc.variables['pm'][:]
    pn = nc.variables['pn'][:]
    dV = dz/(pm*pn)

    rho0 = nc.variables['rho0'][:]
    rho = nc.variables['rho'][tidx]
    b = 9.8*(rho0-(1000.0+rho))/rho0
    bm = b.mean(axis=-1)[:, :, np.newaxis]
    bp = b - bm

    u = nc.variables['u'][tidx]
    um = u.mean(axis=-1)[:, :, np.newaxis]
    up = u - um

    vm = 0.0
    vp = nc.variables['v'][tidx]
    
    wp = nc.variables['w'][tidx, 1:-1]
    
    up, vp = octant.tools.shrink(up, vp)
    Uy = np.diff(um, axis=-2)/1000.0
    dVu = octant.tools.shrink(dV, up.shape)
    
    # <u' v'> U_y
    HMF[tidx, :] = np.sum( np.sum(dVu*up*vp, axis=-1) * Uy[:, :, 0], axis=0) 

    zr = octant.roms.nc_depths(nc, 'rho')[tidx]
    Uz = np.diff(um, axis=0)/np.diff(zr, axis=0)
    
    wp_b, bp_b = octant.tools.shrink(wp, bp)
    dV_b = octant.tools.shrink(dV, wp_b.shape)
    
    # <b' w'>
    VBF[tidx, :] = np.sum( np.sum(dV_b*bp_b*wp_b, axis=-1) , axis=0) 
    
    up_w, wp_w = octant.tools.shrink(up, wp)
    Uz_w = octant.tools.shrink(Uz, up_w.shape)
    dV_w = octant.tools.shrink(dV, up_w.shape)

    # <u' w'> U_z
    VMF[tidx, :] = np.sum( np.sum(dV_w*up_w*wp_w, axis=-1) * Uz_w[:, :, 0], axis=0) 

    bp = octant.tools.shrink(bp, vp.shape)
    
    # <b' v'> (B_x / B_z)
    HBF[tidx, :] = np.sum( np.sum(dVu*bp*vp, axis=-1) * (1e-3/delta), axis=0) 

#
np.save('HMF', HMF.sum(axis=-1)/dV.sum())
np.save('VMF', VMF.sum(axis=-1)/dV_w.sum())

np.save('HBF', HBF.sum(axis=-1)/dVu.sum())
np.save('VBF', VBF.sum(axis=-1)/dV.sum())

#
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)

ax.plot(time, VBF.sum(axis=-1)/dV.sum(), '-r', lw=2)
ax.plot(time, HMF.sum(axis=-1)/dV.sum(), '-g', lw=2)
ax.plot(time, VMF.sum(axis=-1)/dV_w.sum(), '-b', lw=2)
ax.grid(True)
ax.set_xlabel('time [days]')
ax.text(9, 2e-8, r"$<w' b'>$", color='r')

# ax.set_ylabel('Mean energy flux [m$^2$ s$^{-3}$]')
# ax.text(15.5, 10.5, 'Horizontal eddy buoyancy flux', color='r')
# ax.text(15.5, 0.5, 'Horizontal Reynolds stress', color='k')
#
# plt.savefig('energy_conversion.pdf')
