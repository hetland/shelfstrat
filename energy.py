
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import netCDF4

import octant
import octant.roms


ncfilename = 'simulations/shelfstrat_M2_4.47e-07_N2_1.00e-04_f_1.00e-04/shelfstrat_his.nc'

nc = netCDF4.Dataset(ncfilename)

Nskip = 2

t = nc.variables['ocean_time'][::Nskip]/86400.0

eke = []
mke = []
tke = []

ekep = []
mkep = []
tkep = []

dx = 1.0 / nc.variables['pm'][:]
dy = 1.0 / nc.variables['pn'][:]

salt0 = nc.variables['salt'][0]
ds = salt0.ptp()
s0 = salt0.max()
s_crit = s0 - 0.05*ds

for tidx in range(len(nc.variables['ocean_time']))[::Nskip]:
    
    salt = nc.variables['salt'][tidx]
    salt_bar = salt[:, :, 1:-1].mean(axis=-1)[:, :, None]
    
    u = nc.variables['u'][tidx]
    u_bar = u[:, :, :-1].mean(axis=-1)[:, :, None]
    up = u - u_bar
    vp = nc.variables['v'][tidx]
    up = octant.tools.shrink(up, (30, 128, 256))
    vp = octant.tools.shrink(vp, (30, 128, 256))
    
    u = octant.tools.shrink(u, (30, 128, 256))
    
    # zr = octant.roms.nc_depths(nc, grid='rho')[tidx]
    zw = octant.roms.nc_depths(nc, grid='w')[tidx]
    dz = np.diff(zw, axis=0)
    dV = (dx*dy*dz)[:, 1:-1, 1:-1]
    dV_bar = dV[:, :, :-1].mean(axis=-1)[:, :, None]
    
    eke.append( ((up**2 + vp**2)*dV).sum() / dV.sum() )
    tke.append( ((u**2 + vp**2)*dV).sum() / dV.sum() )
    mke.append( (u_bar[:, 1:-1, :]**2 *dV_bar).sum() / dV_bar.sum() )
    
    idx = salt[:, 1:-1, 1:-1] < s_crit
    ekep.append( ((up**2 + vp**2)*dV)[idx].sum() / dV[idx].sum() )
    tkep.append( ((u**2 + vp**2)*dV)[idx].sum() / dV[idx].sum() )
    
    idx = salt_bar[:, 1:-1, :] < s_crit
    mkep.append( (u_bar[:, 1:-1, :]**2 *dV_bar)[idx].sum() / dV_bar[idx].sum() )


eke = np.asarray(eke)
mke = np.asarray(mke)
tke = np.asarray(tke)

ekep = np.asarray(ekep)
mkep = np.asarray(mkep)
tkep = np.asarray(tkep)
