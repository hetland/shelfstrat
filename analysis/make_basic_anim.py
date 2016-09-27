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
parser.add_argument('--prop', type=str, default='rho', help='property to plot')
args = parser.parse_args()

prop = args.prop
nc = netCDF4.Dataset(args.filename)

x = nc.variables['x_rho'][:]/1000.0
y = nc.variables['y_rho'][:]/1000.0
h = nc.variables['h'][:, args.xidx]
t = nc.variables['ocean_time'][:]

fig = plt.figure(figsize=(10.0, 3.24), dpi=100)

rootdir = os.path.dirname(args.filename)

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
    
    if prop == 'salt':
        cmap = plt.cm.YlGnBu_r
        clabel_str = r'Salinity [g kg$^{-1}$]'
    elif prop == 'temp':
        cmap = plt.cm.OrRd
        clabel_str = r'Temperature [$^\circ$C]'
    else: # prop == 'rho'
        cmap = plt.cm.BuPu_r
        clabel_str = r'Density [kg m$^{-3}$]'
    
    p = nc.variables[prop][tidx, -1]
    if tidx == 0:
        vmin = nc.variables[prop][0].min()
        vmax = nc.variables[prop][0].max()
        print 'vmin, vmax = ', vmin, vmax
        cticks = np.linspace(vmin, vmax, 20)
        print cticks
    
    ###########
    # Plan plot
    ax_plan.contourf(x, y, p, cticks, extend='both', vmin=vmin, vmax=vmax, cmap=cmap)
    ax_plan.set_xlabel('Along-shore distance [km]', fontsize=10)
    ax_plan.set_ylabel('Cross-shore distance [km]', fontsize=10)
    ax_plan.set_title(timestr)
    ax_plan.set_aspect(1.0)
    
    ####################
    # Cross-section plot
    
    ax_xsec = fig.add_axes([0.65, 0.05, 0.33, 0.88])
    
    zrx = octant.roms.nc_depths(nc, grid='rho')[tidx, :, :, args.xidx]
    rhox = nc.variables[prop][tidx, :, :, args.xidx]
    
    yx = nc.variables['y_rho'][:, args.xidx]/1000.0
    xx = nc.variables['x_rho'][:, args.xidx]/1000.0
    ax_plan.plot(xx, yx, '--k')
    
    yxx = yx * np.ones_like(rhox)
    # print 'rhox.min, rhox.max = ', rhox.min(), rhox.max()
    cnt = ax_xsec.contourf(yxx, zrx, rhox, cticks, extend='both', vmin=vmin, vmax=vmax, cmap=cmap)
    ax_xsec.set_xlabel('Along-shore distance [km]', fontsize=10)
    ax_xsec.set_ylabel('Depth [m]', fontsize=10)
    cb = plt.colorbar(cnt, orientation='horizontal', ticks=[vmin, vmax])
    cb.set_label(clabel_str, fontsize=10)
    cb.ax.set_xticklabels([str(np.round(vmin, 2)), str(np.round(vmax, 2))])
    ax_xsec.fill( np.hstack((yx, yx[-1], yx[0])), np.hstack((-h, -h.max(), -h.max())), '0.8', lw=0)
    ax_xsec.plot(yx, -h, '-k', lw=2)
    
    fig.savefig(os.path.join(frames_dir, '%s_frame_%04d.png' % (prop, tidx)), dpi=100)
    print ' ... saved %s_frame_%04d.png at %s' % (prop, tidx, timestr)

frames = os.path.join(frames_dir, prop + '_frame_%04d.png')
outfile = os.path.join(rootdir, prop + '_anim.mp4')
print 'ffmpeg -y -r 10 -i %s -c:v libx264 -crf 15 %s' % (frames, outfile)
os.system('/usr/local/bin/ffmpeg -y -r 10 -i %s -c:v libx264 -pix_fmt yuv420p -crf 15 %s' % (frames, outfile))

