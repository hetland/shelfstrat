
import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

def get_suffix(**kwargs):
    keys = kwargs.keys()
    keys.sort()
    suffix=''
    for key in keys:
        suffix += '_%s_%4.2e' % (str(key), kwargs[key])
    
    return suffix
    

fig, axs = plt.subplots(3, 3, figsize=(10, 6))

def plot_frame(tidx=-1, nframe=0):
    n=0
    for M2 in np.logspace(-6, -8, 3):
        for N2 in np.logspace(-4, -6, 3):
            axs.flat[n].cla()
            suffix = get_suffix(M2=M2, N2=N2)
            rootdir = './simulations/shelfstrat%s' % suffix
            nc = netCDF4.Dataset(os.path.join(rootdir, 'shelfstrat_his.nc'))
            x = nc.variables['x_rho'][:]/1000.0
            y = nc.variables['y_rho'][:]/1000.0
            try:
                ds = nc.variables['salt'][tidx, 0] - nc.variables['salt'][tidx, -1]
            except:
                axs.flat[n].set_xticklabels([])
                axs.flat[n].set_yticklabels([])
                n += 1
                continue
            if M2 == 1e-6:
                Smax = 5
            elif M2 == 1e-7:
                Smax = 0.5
            else:
                Smax = 0.05
            axs.flat[n].contourf(x, y, ds, np.linspace(0, Smax, 21))
            axs.flat[n].set_aspect(1.0)
            axs.flat[n].set_title(suffix.replace('_', ' '), fontsize=8)
            axs.flat[n].set_xlim(0, x.max())
            axs.flat[n].set_ylim(0, y.max())
            if n != 6:
                axs.flat[n].set_xticklabels([])
                axs.flat[n].set_yticklabels([])
            print rootdir
            n += 1
            nc.close()
    
    fig.savefig('sss_frames_%04d.png' % nframe)


nframe = 0
for tidx in range(321)[::10]:
    plot_frame(tidx, nframe=nframe)
    print 'tidx = ', tidx
    nframe += 1