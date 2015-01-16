
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import netCDF4

import octant

from case_dictionary import case_dictionary

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory with top-level cases listed')
args = parser.parse_args()

cases = case_dictionary(args.directory)

fig, axs = plt.subplots(4, 5, figsize=(13.5, 6))
axs = np.flipud(axs)

ax_crnr = axs[0, 0]

# files = [cases.find(Ri=Ri, delta=delta, N2=1e-4)
#                 for delta in [0.1, 0.2, 0.3, 0.5]
#                 for Ri in [1, 2, 3, 5, 10]]

files = [cases.find(Ri=Ri, S=S, N2=1e-4) 
            for S in [0.1, 0.2, 0.3, 0.5]
            for Ri in [1, 2, 3, 5, 10]]


delta_colors = {'0.1': (1.0, 0.0, 0.0), 
                '0.2': (0.8, 0.1, 0.2),
                '0.3': (0.2, 0.1, 0.8),
                '0.5': (0.0, 0.0, 1.0),}

delta_linestyles = {'0.1': '-', 
                    '0.2': '--',
                    '0.3': '-.',
                    '0.5': '-',}


def polyfit_omega(n=6):
    'fit an order-n polynomial to the maximum growth rate as a function of delta, the slope parameter.'
    delta, mu = np.mgrid[-1.2:2.2:1001j, 0:4.2:1001j]
    
    tmu = np.tanh(mu)
    omega = np.sqrt( (1.0+delta)*(mu - tmu)/tmu 
                     -0.25*(delta/tmu + mu)**2 ).real
    omega2 = (1.0+delta)*(mu - tmu)/tmu - 0.25*(delta/tmu + mu)**2
    
    omega = np.ma.masked_where(np.isnan(omega), omega)
    
    omega_max = omega.max(axis=1)
    idx = np.where(~omega_max.mask)
    
    omega_max = omega_max[idx]
    delta = delta[:, 0][idx]
    
    p = np.polyfit(delta, omega_max, n)
    
    return p

omega_poly = polyfit_omega()

Ris = []
deltas = []
ekes = []

fig_all = plt.figure()
ax_all_normed = fig_all.add_subplot(211)
ax_all_unnormed = fig_all.add_subplot(212)

for ax, file in zip(axs.flat, files):
    
    hisfilename = os.path.join(args.directory, file[0], 'shelfstrat_his.nc')
    params = cases[file[0]]
    print hisfilename
    
    omega = np.polyval(omega_poly, params['delta'])
    omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri'])
    
    # timescale = 1.0 * omega_dim_ref / omega_dim
    
    nc = netCDF4.Dataset(hisfilename)
    time = nc.variables['ocean_time'][:] / 86400.0
    
    # x = nc.variables['x_rho'][:]/1000.0
    # y = nc.variables['y_rho'][:]/1000.0
    
    u = nc.variables['u'][:, -1, :, :]
    v = nc.variables['v'][:, -1, :, :]
    
    u, v = octant.tools.shrink(u, v)
    
    umean = u.mean(axis=-1)[:, :, None]
    up = u - umean
    vp = v
    
    tke = 0.5*(u**2 + v**2)
    eke = 0.5*(up**2 + v**2)
    mke = 0.5*(umean**2)
    
    norm = mke.mean(axis=-1).mean(axis=-1)[0]
    
    ax.plot(time*omega_dim, tke.mean(axis=-1).mean(axis=-1)/norm, '-k')
    ax.plot(time*omega_dim, eke.mean(axis=-1).mean(axis=-1)/norm, '-r')
    ax.plot(time*omega_dim, mke.mean(axis=-1).mean(axis=-1)/norm, '-b')
    
    if params['delta'] == 0.1:
        ax_crnr.plot(time*omega_dim, eke.mean(axis=-1).mean(axis=-1)/norm/np.sqrt(params['Ri']), '-r', lw=0.25)
    
    delta_str = '%0.1f' % params['delta']
    ax_all_normed.plot(time*omega_dim, eke.mean(axis=-1).mean(axis=-1)/norm/np.sqrt(params['Ri']), 
                    linestyle=delta_linestyles[delta_str], color=delta_colors[delta_str])
    
    ax_all_unnormed.plot(time, eke.mean(axis=-1).mean(axis=-1)/norm/np.sqrt(params['Ri']), 
                    linestyle=delta_linestyles[delta_str], color=delta_colors[delta_str])
    
    ax.plot([1, 2], [1, 1], 'k-', alpha=0.5, lw=4)
    
    ax.set_ylim(0, 2)
    ax.set_xlim(0, 40.0)
    
    Ris.append(params['Ri'])
    deltas.append(params['delta'])
    ekes.append( (eke.mean(axis=-1).mean(axis=-1)/norm).max() )
    
    if ax == ax_crnr:
        ax.set_xlabel('Time [days]')
        ax.set_ylabel('Normalized energy')
    else:
        # ax.set_xticklabels([])
        ax.set_yticklabels([])

ref_timescale = 10.0 # days
ref_delta = 0.1
ref_Ri = 1.0
ref_f = 1e-4
omega_ref = np.polyval(omega_poly, ref_delta) # non-dim
omega_dim_ref = 86400.0 * omega_ref * ref_f / np.sqrt(ref_Ri)  # rad/days
ax_all_normed.plot([ref_timescale*omega_dim_ref, ref_timescale*omega_dim_ref], [0, 0.8], '--k', lw=3)

fig.subplots_adjust(left=0.05, bottom=0.095, right=0.95, top=0.995,
                  wspace=0.05, hspace=0.05)

Ri = np.array(Ris)
delta = np.array(deltas)
eke = np.array(ekes)
S = delta/np.sqrt(Ri)
p = np.polyfit(np.log(S), eke, 1)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(S, eke, 'ko')
logS = np.linspace(-5, 0, 50)
ax.plot(np.exp(logS), p[0]*logS + p[1], 'k--')
ax.set_xlim(1e-2, 1)
ax.set_ylim(0, 2)
r = np.corrcoef(np.log(S), eke)[0, 1]

plt.show()
