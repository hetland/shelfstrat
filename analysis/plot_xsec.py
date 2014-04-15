
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

import octant
import octant.roms


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='NetCDF filename to plot')
parser.add_argument('--prop', type=str, default='rho', help='property to plot')
parser.add_argument('--xidx', type=int, default=1, help='x-index slice to plot')
parser.add_argument('--tidx', type=int, default=-1, help='time index to plot')
parser.add_argument('--outfile', type=str, default='xsec.png', help='time index to plot')
parser.add_argument('--vmin', type=float, default=0, help='min property')
parser.add_argument('--vmax', type=float, default=32.0, help='max property')
args = parser.parse_args()

nc = netCDF4.Dataset(args.filename)

y = nc.variables['y_rho'][:, args.xidx]

zr = octant.roms.nc_depths(nc, grid='rho')[args.tidx, :, :, args.xidx]
p = nc.variables[args.prop][args.tidx, :, :, args.xidx]

y = y * np.ones_like(p)
plt.contourf(y, zr, p, 20, extend='both', vmin=args.vmin, vmax=args.vmax)
plt.savefig(args.outfile)

