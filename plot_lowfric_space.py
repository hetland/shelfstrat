import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import netCDF4

from case_dictionary import case_dictionary

lowfric_cases = case_dictionary('sims_lowfric')
vlowfric_cases = case_dictionary('sims_vlfric')

#tidx = 112   # 14 days
tidx = -1

fig, axs = plt.subplots(6, 2, figsize=(6, 10))
axs = np.flipud(axs)

ax_crnr = axs[0, 0]

lowfric_files = [os.path.join('sims_lowfric', lowfric_cases.find(Ri=Ri, delta=delta, N2=1e-4)[0])
                 for delta in [0.1, 0.2, 0.3, 0.5, 1.0]
                 for Ri in [1, 10]]

vlowfric_files = [os.path.join('sims_vlfric', vlowfric_cases.find(Ri=Ri, delta=delta, N2=1e-4)[0])
                 for delta in [1.0]
                 for Ri in [1, 10]]

files = lowfric_files + vlowfric_files

for ax, file in zip(axs.flat, files):
    
    hisfilename = os.path.join(file, 'shelfstrat_his.nc')
    
    if 'sims_lowfric' in file:
        params = lowfric_cases[os.path.split(file)[-1]]
    else:
        params = vlowfric_cases[os.path.split(file)[-1]]
    
    print hisfilename
    
    nc = netCDF4.Dataset(hisfilename)
    time = nc.variables['ocean_time'][:]
    
    x = nc.variables['x_rho'][:]/1000.0
    y = nc.variables['y_rho'][:]/1000.0
    sss = nc.variables['salt'][tidx, -1, :, :]
    
    ax.contourf(x, y, sss, 10, cmap=plt.cm.YlGnBu_r)
    ax.set_aspect(1.0)
    
    Rd = np.sqrt(params['N2']) * 50.0 / params['f']
    Uscale = np.sqrt(params['N2']/params['Ri'])*50.0
    Ladv = Uscale / params['f']
    
    plt.text(.05, 0.95, r''' $Ri$=%4.2f
$\delta$=%4.2f
$S$=%4.2f
tidx=%d
''' % (params['Ri'], params['delta'], params['S'], len(time)),
          horizontalalignment='left',
          verticalalignment='top',
          fontsize=10,
          transform = ax.transAxes)
#
#  
#     plt.text(0.05, 0.95, r''' $R_d$=%4.2f km
# $L_{adv}$=%4.2f km
# ''' % (Rd/1000.0, Ladv/1000.0),
#          horizontalalignment='left',
#          verticalalignment='top',
#          fontsize=10,
#          transform = ax.transAxes)

#     plt.text(0.05, 0.95, r''' f=%4.2e
# M$^2$=%4.2e
# ''' % (params['f'], params['M2']),
#          horizontalalignment='left',
#          verticalalignment='top',
#          fontsize=8,
#          transform = ax.transAxes)
    
    if ax == ax_crnr:
        ax.set_xlabel('Along-shore distance [km]')
        ax.set_ylabel('Offshore distance [km]')
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        

plt.subplots_adjust(left=0.15, bottom=0.095, right=0.95, top=0.995,
                  wspace=0.05, hspace=0.05)
plt.show()
