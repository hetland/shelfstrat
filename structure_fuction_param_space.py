
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from scipy.stats import binned_statistic

from case_dictionary import case_dictionary

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory with top-level cases listed')
args = parser.parse_args()

cases = case_dictionary(args.directory)

tidx = 112   # 14 days

fig, axs = plt.subplots(4, 5, figsize=(13.5, 6))
axs = np.flipud(axs)

ax_crnr = axs[0, 0]

files = [cases.find(Ri=Ri, delta=delta, N2=1e-4) 
            for delta in [0.1, 0.2, 0.3, 0.5]
            for Ri in [1, 2, 3, 5, 10]]


def plot_slopes(power, N=10, **kwargs):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    yints = np.logspace(np.log10(ymin)-10, np.log10(ymax), N+10)
    for yint in yints:
        x = np.logspace(np.log10(xmin), np.log10(xmax), 10)
        y = yint*(x/x[0])**power
        ax.loglog(x, y, **kwargs)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


for ax, file in zip(axs.flat, files):
    
    hisfilename = os.path.join(args.directory, file[0], 'shelfstrat_his.nc')
    params = cases[file[0]]
    
    print hisfilename
    
    nc = netCDF4.Dataset(hisfilename)
    time = nc.variables['ocean_time'][:]
    
    xu = nc.variables['x_u'][50, :]
    u = nc.variables['u'][tidx, -1, 50, :]
    Uscale = np.sqrt(params['N2']/params['Ri'])*50.0
    
    x = xu
    var = u / Uscale
    N = len(u)
    
    var = (u[:, None] - u[None, :])**2
    dist = np.abs(x[:, None] - x[None, :])
    mask = np.triu(np.ones((N, N), 'd'))
    dist = dist.flat[mask.flat==1]
    var = var.flat[mask.flat==1]
    v, bins, bn = binned_statistic(dist, var, bins=np.arange(0, 1e5, 1e3))
    
    Rd = np.sqrt(params['N2'])*50.0/params['f']
    
    ax.loglog(bins[1:], v, color='k', lw=3)
    plot_slopes(4, ls='-', color='b', lw=0.5)
    plot_slopes(2.0/3.0, ls='-', color='r', lw=0.5)
    
    ax.plot((Rd, Rd), (1e-7, 1.0))
    
    plt.text(0.05, 0.95, r''' f=%4.2e
M$^2$=%4.2e
''' % (params['f'], params['M2']),
         horizontalalignment='left',
         verticalalignment='top',
         fontsize=8,
         transform = ax.transAxes)
    
    if ax == ax_crnr:
        ax.set_xlabel('Separation [m]')
        ax.set_ylabel(r'<$\delta u(r)>')
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
    ax.set_xlim(1e3, 1e5)
    ax.set_ylim(1e-7, 0.1)

plt.subplots_adjust(left=0.05, bottom=0.095, right=0.95, top=0.995,
                  wspace=0.05, hspace=0.05)
plt.show()
