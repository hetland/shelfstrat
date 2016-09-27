
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import netCDF4

import octant

from case_dictionary import case_dictionary

# supress warnings that occur in division in the polyfit_omega function
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory with top-level cases listed')
args = parser.parse_args()

cases = case_dictionary(args.directory)


# files = [cases.find(Ri=Ri, delta=delta, N2=1e-4)
#             for delta in [0.1, 0.2, 0.3, 0.5]
#             for Ri in [1, 2, 3, 5, 10]]
# fig, axs = plt.subplots(4, 5, figsize=(13.5, 6))
# normalize_time = False

files = [cases.find(Ri=Ri, S=S, N2=1e-4)
            for S in [0.1, 0.2, 0.3, 0.5]
            for Ri in [1, 2, 3, 5, 10]]
fig, axs = plt.subplots(4, 5, figsize=(13.5, 6))
normalize_time = True

axs = np.flipud(axs)
ax_crnr = axs[0, 0]


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

# ref_timescale = 7.0 # days
# ref_delta = 0.1
# ref_Ri = 1.0
# ref_f = 1e-4
# omega = np.polyval(omega_poly, ref_delta) # non-dim
# omega_dim = 86400.0 * omega * ref_f / np.sqrt(ref_Ri)  # rad/days
# timescale_factor = ref_timescale * omega_dim

for ax, file in zip(axs.flat, files):
    
    hisfilename = os.path.join(args.directory, file[0], 'shelfstrat_his.nc')
    params = cases[file[0]]
    
    omega = np.polyval(omega_poly, params['delta'])
    omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri']) # rad/days
    
    if normalize_time:
        # timescale = timescale_factor / omega_dim
        omega = np.polyval(omega_poly, params['delta'])
        omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri']) # rad/days
        timescale = 50.0 * np.sqrt(params['S']) / omega_dim
    else:
        timescale = ref_timescale
    
    print(hisfilename)
    
    nc = netCDF4.Dataset(hisfilename)
    time = nc.variables['ocean_time'][:] / 86400.0
    
    tidx = np.where( time >= timescale )[0]
    if len(tidx) == 0:
        tidx = len(time) - 1
    else:
        tidx = tidx.min()
    
    print('   {0:d}/{1:d} -- {2:f}'.format(tidx, len(time), time[tidx]))
    
    x = nc.variables['x_psi'][:]/1000.0
    y = nc.variables['y_psi'][:]/1000.0
    
    u = nc.variables['u'][tidx, -1, :, :]
    v = nc.variables['v'][tidx, -1, :, :]
    
    vort = (np.diff(v, axis=-1) - np.diff(u, axis=-2))/1000.0
    
    ax.contourf(x, y, vort/params['f'], 11, cmap=plt.cm.RdBu_r, vmin=-3, vmax=3)
    ax.set_aspect(1.0)
    
    Rd = np.sqrt(params['N2']) * 50.0 / params['f']
    Uscale = np.sqrt(params['N2']/params['Ri'])*50.0
    Ladv = Uscale / params['f']
    
    ax.text(0.05, 0.9, 
         # '$R_d$=%5.2f km\n$L_{adv}$=%5.2f km\n$T$=%5.2f days\nTo=%5.2f days'% (Rd/1000.0, Ladv/1000.0, time[tidx],timescale),
         # '$R_d$=%5.2f km\n$L_{adv}$=%5.2f km\n$T$=%5.2f days'% (Rd/1000.0, Ladv/1000.0, time[tidx]),
         '$R_d$=%5.2f km\n$T$=%6.2f days\n$T_0$=%6.2f'% (Rd/1000.0, time[tidx], timescale),
         horizontalalignment='left', 
         verticalalignment='top',
         transform=ax.transAxes,
         fontsize=8)
    
    def expsplit(qlist):
        res = ()
        for q in qlist:
            qstr = '%e' % q
            res += tuple(map(float, qstr.split('e')))
        return res
    
    paramstrs = expsplit([params['M2'], params['f']]) + (params['delta'],)
    ax.text(0.5, 0.9, 
         '$M^2$ = %5.2fx10$^{%d}$ s$^{-2}$\n$f$ = %5.2fx10$^{%d}$ s$^{-1}$\n$\delta$ = %5.2f' % paramstrs,
         horizontalalignment='left', 
         verticalalignment='top',
         transform=ax.transAxes,
         fontsize=8)


    if ax == ax_crnr:
        ax.set_xlabel('Along-shore distance [km]')
        ax.set_ylabel('Cross-shore distance [km]')
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    


plt.subplots_adjust(left=0.05, bottom=0.095, right=0.95, top=0.995,
                  wspace=0.05, hspace=0.05)
plt.savefig('Ri_S_sss.png', dpi=300)
# plt.savefig('Ri_S_sss.pdf', dpi=300)
plt.show()
