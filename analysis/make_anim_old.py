
import matplotlib 
matplotlib.rc('xtick', labelsize=8) 
matplotlib.rc('ytick', labelsize=8) 

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import argparse
import os
from matplotlib.ticker import FormatStrFormatter

import octant
import octant.roms


parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='NetCDF filename to plot')
parser.add_argument('--xidx', type=int, default=100, help='x-index slice to plot')
args = parser.parse_args()

salt_cdict = {'red':  ((0.00, 0.4,  0.4),
                       (0.35, 0.3,  0.3),
                       (0.66, 1.0,  1.0),
                       (0.85, 0.9,  0.9),
                       (0.93, 0.75,  0.75),
                       (1.00, 0.83, 0.83)),
             'green': ((0.00,  0.4, 0.4),
                       (0.125, 0.3, 0.3),
                       (0.375, 1.0, 1.0),
                       (0.64,  1.0, 1.0),
                       (0.75,  0.5, 0.5),
                       (0.93,  0.5, 0.5),
                       (1.00,  0.8, 0.8)),
             'blue':  ((0.00, 0.7, 0.7),
                       (0.11, 1.0, 1.0),
                       (0.34, 1.0, 1.0),
                       (0.65, 0.0, 0.0),
                       (0.85,  0.6, 0.6),
                       (1.00, 0.8, 0.8))}
                       
salt_cmap = plt.matplotlib.colors.LinearSegmentedColormap('salt_cmap', salt_cdict, 256)

nc = netCDF4.Dataset(args.filename)

x = nc.variables['x_rho'][:]/1000.0
y = nc.variables['y_rho'][:]/1000.0
h = nc.variables['h'][:, args.xidx]
t = nc.variables['ocean_time'][:]

fig = plt.figure(figsize=(10.0,3.25))

rootdir = os.path.dirname(args.filename)

case = rootdir.split('_')[1:]
keys = case[::2]
vals = [float(val) for val in case[1::2]]
params = dict(zip(keys, vals))

Ri = params['N2'] * params['f']**2 / params['M2']**2

phi = params['M2'] * 1e-3 * np.sqrt(Ri) / params['f']**2

frames_dir = os.path.join(rootdir, 'frames')
if not os.path.exists(frames_dir):
        os.makedirs(frames_dir)

print '### Saving frames to %s' % frames_dir

for tidx in range(len(t)):
    
    days = np.floor(t[tidx] / 86400.0)
    hours = np.floor( (t[tidx]-days*86400) / 3600.0)
    minutes = 0
    timestr = 'day %02d %02d:%02d' % (days, hours, minutes)
    
    fig.clf()
    ax_plan = fig.add_axes([0.05, 0.1, 0.5, 0.8])
    
    rho = nc.variables['rho'][tidx, -1]
    
    if tidx == 0:
        
        vmin = nc.variables['rho'][0].min()
        vmax = nc.variables['rho'][0].max()
        print 'vmin, vmax = ', vmin, vmax
        cticks = np.linspace(vmin, vmax, 20)
        print cticks
    
    ###########
    # Plan plot
    ax_plan.contourf(x, y, rho, cticks, extend='both', vmin=vmin, vmax=vmax, cmap=salt_cmap)
    ax_plan.set_xlabel('Along-shore distance [km]', fontsize=10)
    ax_plan.set_ylabel('Cross-shore distance [km]', fontsize=10)
    ax_plan.set_title(timestr)
    ax_plan.set_aspect(1.0)
    
    ####################
    # Cross-section plot
    
    ax_xsec = fig.add_axes([0.65, 0.05, 0.33, 0.88])
    
    zrx = octant.roms.nc_depths(nc, grid='rho')[tidx, :, :, args.xidx]
    rhox = nc.variables['rho'][tidx, :, :, args.xidx]
    
    yx = nc.variables['y_rho'][:, args.xidx]/1000.0
    xx = nc.variables['x_rho'][:, args.xidx]/1000.0
    ax_plan.plot(xx, yx, '--k')
    
    yxx = yx * np.ones_like(rhox)
    print 'rhox.min, rhox.max = ', rhox.min(), rhox.max()
    cnt = ax_xsec.contourf(yxx, zrx, rhox, cticks, extend='both', vmin=vmin, vmax=vmax, cmap=salt_cmap)
    ax_xsec.set_xlabel('Along-shore distance [km]', fontsize=10)
    ax_xsec.set_ylabel('Depth [m]', fontsize=10)
    cb = plt.colorbar(cnt, orientation='horizontal', ticks=[vmin, vmax])
    cb.set_label(r'Density [kg m$^{-3}$]', fontsize=10)
    cb.ax.set_xticklabels([str(np.round(vmin, 2)), str(np.round(vmax, 2))])
    ax_xsec.fill( np.hstack((yx, yx[-1], yx[0])), np.hstack((-h, -h.max(), -h.max())), '0.8', lw=0)
    ax_xsec.plot(yx, -h, '-k', lw=2)
    ax_xsec.text(0.05, 0.4,r'''Ri = %5.3f
$\phi$  = %5.3f

f  = %4.2e
N$^2$ = %4.2e
M$^2$ = %4.2e''' % (Ri, phi, params['f'], params['N2'], params['M2']),
             horizontalalignment='left', verticalalignment='center',
             transform = ax_xsec.transAxes, fontsize=8)
    
    fig.savefig(os.path.join(frames_dir, 'frame_%04d.png' % tidx))
    print ' ... saved frame_%04d.png at %s' % (tidx, timestr)

frames = os.path.join(frames_dir, 'frame_%04d.png')
outfile = os.path.join(rootdir, 'rho_anim.mp4')
print 'ffmpeg -y -sameq -r 10 -i %s %s' % (frames, outfile)
os.system('ffmpeg -y -sameq -r 10 -i %s %s' % (frames, outfile))

