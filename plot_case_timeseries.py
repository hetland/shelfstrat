
import os, sys
import argparse

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import matplotlib.ticker


import octant
import octant.roms

from case_dictionary import case_dictionary

# supress warnings that occur in division in the polyfit_omega function
import warnings
warnings.filterwarnings("ignore")


def get_case_params(casedir):
    case = casedir.split('_')[1:]  # remove leading description 'shelfstrat'
    keys = case[::2]
    vals = [float(val) for val in case[1::2]]
    params = dict(zip(keys, vals))
    params['Ri'] = params['N2'] * params['f']**2 / params['M2']**2
    params['S'] = params['M2'] * 1e-3 * np.sqrt(params['Ri']) / params['f']**2
    params['delta'] = np.sqrt(params['Ri']) * params['S']
    return params

use_saved_histograms = False


parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory with top-level cases listed')
args = parser.parse_args()

# cases = case_dictionary(args.directory)
# file = cases.find(Ri=10.0, S=0.3, N2=1e-4)
# # file = cases.find(Ri=1.0, S=0.1, N2=1e-4)

# params = cases[file[0]]
# hisfilename = os.path.join(args.directory, file[0], 'shelfstrat_his.nc')

file = args.directory.split('/')[-1]
params = get_case_params(file)
hisfilename = os.path.join(args.directory, 'shelfstrat_his.nc')
print params

###########################################################################
nc = netCDF4.Dataset(hisfilename)

x = nc.variables['x_rho'][:]/1000.0
y = nc.variables['y_rho'][:]/1000.0
t = nc.variables['ocean_time'][:]/86400.0

# fig, axs = plt.subplots(3, 3, figsize=(10, 4.8))
# # axs = np.flipud(axs)
# ax_crnr = axs[2, 0]
#
# for ax, tidx in zip(axs.flat, range(0, 18*8, 16)):
#     sss = nc.variables['salt'][tidx, -1, :, :]
#     cnt = ax.contourf(x, y, sss, 10, cmap=plt.cm.YlGnBu_r)
#     ax.set_aspect(1.0)
#     plt.text(0.05, 0.95, '%5.2f days'%t[tidx],
#              horizontalalignment='left',
#              verticalalignment='top',
#              transform=ax.transAxes)
#
#     if ax == ax_crnr:
#         ax.set_xlabel('Along-shore distance [km]')
#         ax.set_ylabel('Cross-shore distance [km]')
#     else:
#         ax.set_xticklabels([])
#         ax.set_yticklabels([])
#
# cax = fig.add_axes((0.92, 0.25, 0.01, 0.5))
# plt.colorbar(cnt, cax=cax).set_label(r'Surface Salinity [g kg$^{-1}$]')
# plt.subplots_adjust(left=0.07, bottom=0.095, right=0.90, top=0.995,
#                   wspace=0.05, hspace=0.05)
# # plt.show()

###########################################################################
# property histograms.

pm = nc.variables['pm'][:]
pn = nc.variables['pn'][:]

f = nc.variables['f'][1, 1]
alpha = 1e-3

bins_M2 = np.linspace(-8.5, -2.5, 301)
bins_M2_c = 0.5 * (bins_M2[1:] + bins_M2[:-1])

bins_N2 = np.linspace(-6, 0, 201)
bins_N2_c = 0.5 * (bins_N2[1:] + bins_N2[:-1])

bins_Ri = np.linspace(-3, 3, 201)
bins_Ri_c = 0.5 * (bins_Ri[1:] + bins_Ri[:-1])

bins_S = np.linspace(-2, 1, 201)
bins_S_c = 0.5 * (bins_S[1:] + bins_S[:-1])

bins_delta = np.linspace(-2, 2, 201)
bins_delta_c = 0.5 * (bins_delta[1:] + bins_delta[:-1])

M2hist = []
N2hist = []
Rihist = []
Shist = []
deltahist = []

skip = int(np.floor( float(len(t))/80.0 ))
if skip == 0:
    skip = 1

skip=1

timesteps = len(t[::skip])
print timesteps

if use_saved_histograms:
    M2hist = np.load('M2hist.npy')
    N2hist = np.load('N2hist.npy')    
    Rihist = np.load('Rihist.npy')    
    Shist = np.load('Shist.npy')    
    deltahist = np.load('deltahist.npy')
else:
    for tidx in range(len(t))[::skip]:
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
        S = np.sqrt(N2) * alpha / f
        delta = S * np.sqrt(Ri)
    
        dz = np.diff(zr, axis=0)
        dV = dz / (pm * pn)
        dV = octant.tools.shrink(dV, shp)
    
        M2hist.append( np.histogram(np.log10(M2[:, :46, :]), bins=bins_M2, weights=dV[:, :46, :])[0] )
        N2hist.append( np.histogram(np.log10(N2[:, :46, :]), bins=bins_N2, weights=dV[:, :46, :])[0] )
        Rihist.append( np.histogram(np.log10(Ri[:, :46, :]), bins=bins_Ri, weights=dV[:, :46, :])[0] )
        deltahist.append( np.histogram(np.log10(delta[:, :46, :]), bins=bins_delta, weights=dV[:, :46, :])[0] )
        Shist.append( np.histogram(np.log10(S[:, :46, :]), bins=bins_S, weights=dV[:, :46, :])[0] )
    
        print ' timestep %03d/%d' % (tidx/skip + 1, timesteps)

    M2hist = np.asarray(M2hist)
    N2hist = np.asarray(N2hist)    
    Rihist = np.asarray(Rihist)    
    Shist = np.asarray(Shist)    
    deltahist = np.asarray(deltahist)    
    
    np.save('M2hist', M2hist)
    np.save('N2hist', N2hist)
    np.save('Rihist', Rihist)
    np.save('Shist', Shist)
    np.save('deltahist', deltahist)

###################################
# Plotting
t = t[::skip]

majorLocator1   = matplotlib.ticker.MultipleLocator(1)
majorLocator2   = matplotlib.ticker.MultipleLocator(1)
majorLocator3   = matplotlib.ticker.MultipleLocator(1)
majorLocator4   = matplotlib.ticker.MultipleLocator(1)
majorLocator5   = matplotlib.ticker.MultipleLocator(1)
majorFormatter = matplotlib.ticker.FormatStrFormatter('$10^{%d}$')

M2hist = np.ma.masked_where(M2hist ==0, M2hist)
N2hist = np.ma.masked_where(N2hist ==0, N2hist)
Rihist = np.ma.masked_where(Rihist ==0, Rihist)
deltahist = np.ma.masked_where(deltahist ==0, deltahist)
Shist = np.ma.masked_where(Shist ==0, Shist)

fig = plt.figure(figsize=(8,12))
ax1 = fig.add_subplot(611)
ax2 = fig.add_subplot(612)
ax3 = fig.add_subplot(613)
ax4 = fig.add_subplot(614)
ax5 = fig.add_subplot(615)
ax_energy = fig.add_subplot(616)

pcm1 = ax1.contourf(t, bins_M2_c, M2hist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax1.set_ylabel(r'$M^{\,2}$')
ax1.set_xticklabels([])
ax1.set_title(r'M$^2$ = %4.2e    N$^2$ = %4.2e    Ri = %4.2e      $\delta$ = %4.2e' % (params['M2'], params['N2'], params['Ri'], params['delta']))
ax1.plot(t, np.ones_like(t)*np.log10(params['M2']), '-b', lw=3, alpha=0.3)
ax1.yaxis.set_major_locator(majorLocator1)
ax1.yaxis.set_major_formatter(majorFormatter)
ax1.grid(True)

pcm2 = ax2.contourf(t, bins_N2_c, N2hist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax2.set_ylabel(r'$N^{\,2}$')
ax2.set_xticklabels([])
ax2.plot(t, np.ones_like(t)*np.log10(params['N2']), '-b', lw=3, alpha=0.3)
ax2.yaxis.set_major_locator(majorLocator2)
ax2.yaxis.set_major_formatter(majorFormatter)
ax2.grid(True)

pcm3 = ax3.contourf(t, bins_Ri_c, Rihist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax3.set_ylabel(r'$Ri$')
ax3.set_xticklabels([])
ax3.plot(t, np.ones_like(t)*np.log10(params['Ri']), '-b', lw=3, alpha=0.3)
ax3.yaxis.set_major_locator(majorLocator3)
ax3.yaxis.set_major_formatter(majorFormatter)
ax3.grid(True)


pcm4 = ax4.contourf(t, bins_S_c, Shist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax4.set_ylabel(r'$S$')
ax4.plot(t, np.ones_like(t)*np.log10(params['S']), '-b', lw=3, alpha=0.3)
ax4.set_xticklabels([])
ax4.yaxis.set_major_locator(majorLocator4)
ax4.yaxis.set_major_formatter(majorFormatter)
ax4.grid(True)
       
pcm5 = ax5.contourf(t, bins_delta_c, deltahist.T, np.linspace(0, 1e10, 20), cmap=plt.cm.Reds, extend='max')
ax5.set_ylabel(r'$\delta$')
ax5.plot(t, np.ones_like(t)*np.log10(params['delta']), '-b', lw=3, alpha=0.3)
ax5.set_xticklabels([])
ax5.yaxis.set_major_locator(majorLocator5)
ax5.yaxis.set_major_formatter(majorFormatter)
ax5.grid(True)

# Energy stuff.

u = nc.variables['u'][::skip, -1, :, :]
v = nc.variables['v'][::skip, -1, :, :]

u, v = octant.tools.shrink(u, v)

umean = u.mean(axis=-1)[:, :, None]
up = u - umean
vp = v

tke = 0.5*(u**2 + v**2)
eke = 0.5*(up**2 + v**2)
mke = 0.5*(umean**2)

norm = mke.mean(axis=-1).mean(axis=-1)[0]

ax_energy.plot(t, tke.mean(axis=-1).mean(axis=-1)/norm, '-k', lw=2)
ax_energy.plot(t, eke.mean(axis=-1).mean(axis=-1)/norm, '-r', lw=2)
ax_energy.plot(t, mke.mean(axis=-1).mean(axis=-1)/norm, '-b', lw=2)
ax_energy.set_xlabel('Time [days]')
ax_energy.grid(True)

plt.savefig(file + '_histogram.png', dpi=300)
plt.show()
