
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os

import octant
import octant.roms



import argparse
parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str,
                   help='NetCDF filename to plot')
args = parser.parse_args()

rootdir = os.path.dirname(args.filename)

case = rootdir.split('_')[1:]
keys = case[::2]
vals = [float(val) for val in case[1::2]]
params = dict(zip(keys, vals))

alpha = 1e-3

params['Ri'] = params['N2'] * params['f']**2 / params['M2']**2
params['phi'] = params['M2'] * alpha * np.sqrt(params['Ri']) / params['f']**2

bins_M2 = np.linspace(-8., -6., 301)
bins_M2_c = 0.5 * (bins_M2[1:] + bins_M2[:-1])

bins_N2 = np.linspace(-6, 0, 201)
bins_N2_c = 0.5 * (bins_N2[1:] + bins_N2[:-1])

bins_Ri = np.linspace(-3, 3, 201)
bins_Ri_c = 0.5 * (bins_Ri[1:] + bins_Ri[:-1])

bins_phi = np.linspace(-2, 2, 201)
bins_phi_c = 0.5 * (bins_phi[1:] + bins_phi[:-1])

nc = netCDF4.Dataset(args.filename)
time = nc.variables['ocean_time'][:]/86400.0

pm = nc.variables['pm'][:]
pn = nc.variables['pn'][:]

f = nc.variables['f'][1, 1]

M2hist = []
N2hist = []
Rihist = []
phihist = []

step = 20

for tidx in range(len(time))[::step]:
    zr = octant.roms.nc_depths(nc, grid='rho')[tidx]
    
    rho = nc.variables['rho'][tidx]
    M2x, M2y = octant.tools.hgrad(rho, zr, pm, pn)
    N2 = octant.tools.N2(rho, zr, rho_0=nc.variables['rho0'][0])
    
    shp = (N2.shape[0], N2.shape[1]-2, N2.shape[2]-2)
    M2x = octant.tools.shrink(M2x, shp)
    M2y = octant.tools.shrink(M2y, shp)
    N2 = octant.tools.shrink(N2, shp)
    
    M2 = 9.8 * np.sqrt(M2x**2 + M2y**2) / 1025.0
    
    Ri = N2 * f**2 / M2**2
    phi = M2 * alpha * np.sqrt(Ri) / f**2
    
    dz = np.diff(zr, axis=0)
    dV = dz / (pm * pn)
    dV = octant.tools.shrink(dV, shp)
    
    M2hist.append( np.histogram(np.log10(M2[:, :46, :]), bins=bins_M2, weights=dV[:, :46, :])[0] )
    N2hist.append( np.histogram(np.log10(N2[:, :46, :]), bins=bins_N2, weights=dV[:, :46, :])[0] )
    Rihist.append( np.histogram(np.log10(Ri[:, :46, :]), bins=bins_Ri, weights=dV[:, :46, :])[0] )
    phihist.append( np.histogram(np.log10(phi[:, :46, :]), bins=bins_phi, weights=dV[:, :46, :])[0] )
    
    print ' timestep %03d/%d' % (tidx+1, len(time)/step)

M2hist = np.asarray(M2hist)
N2hist = np.asarray(N2hist)    
Rihist = np.asarray(Rihist)    
phihist = np.asarray(phihist)    

###################################
# Plotting

M2hist = np.ma.masked_where(M2hist ==0, M2hist)
N2hist = np.ma.masked_where(N2hist ==0, N2hist)
Rihist = np.ma.masked_where(Rihist ==0, Rihist)
phihist = np.ma.masked_where(phihist ==0, phihist)

fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)

t = time[::step]

pcm1 = ax1.contourf(t, bins_M2_c, M2hist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax1.set_ylabel(r'log$_{10}$(M$^2$)')
ax1.set_xticklabels([])
ax1.set_title(r'M$^2$ = %4.2e    N$^2$ = %4.2e    Ri = %4.2e      $\phi$ = %4.2e' % (params['M2'], params['N2'], params['Ri'], params['phi']))
ax1.plot(t, np.ones_like(t)*np.log10(params['M2']), '-b', lw=3, alpha=0.3)
ax1.grid(True)

pcm2 = ax2.contourf(t, bins_N2_c, N2hist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax2.set_ylabel(r'log$_{10}$(N$^2$)')
ax2.set_xticklabels([])
ax2.plot(t, np.ones_like(t)*np.log10(params['N2']), '-b', lw=3, alpha=0.3)
ax2.grid(True)

pcm3 = ax3.contourf(t, bins_Ri_c, Rihist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax3.set_ylabel(r'log$_{10}$(Ri)')
ax3.set_xticklabels([])
ax3.plot(t, np.ones_like(t)*np.log10(params['Ri']), '-b', lw=3, alpha=0.3)
ax3.grid(True)


pcm4 = ax4.contourf(t, bins_phi_c, phihist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax4.set_ylabel(r'log$_{10}$($\phi$)')
ax4.set_xlabel('Time [days]')
ax4.plot(t, np.ones_like(t)*np.log10(params['phi']), '-b', lw=3, alpha=0.3)
ax4.grid(True)
       
# cax = fig.add_axes([0.93, 0.35, 0.02, 0.3])
# cb = colorbar(pcm1, cax=cax)


case = args.filename.split('/')[1]
plt.savefig('volhist_%s.png' % case, dpi=300)
# plt.show()
    

