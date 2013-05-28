""" This script creates a grid and forcing file for an idealized near-field plume case
"""



import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt
import netCDF4

import octant
import octant.roms



def make_grd(rootdir = '../project/',
             Hmin = 5.0, alpha = 0.001, 
             f = 1e-4,
             dx = 1e3, dy = 1e3,
             shp = (131, 259)):
    
    if not os.path.exists(rootdir):
        print ' ### Making directory %s' % rootdir
        os.mkdir(rootdir)
    
    M, L = shp
    x = np.arange(L, dtype='d') * dx
    y = np.arange(M, dtype='d') * dy
    
    xvert, yvert = np.meshgrid(x, y)
    
    grd = octant.grid.CGrid(xvert, yvert)
    
    # custom set coreolis paremeter (otherwise determined from grd.lat)
    grd.f = f
    
    # create an exponential depth profile, with decay radius Rh, and a value
    # of Hmin and Hmax at the river and ocean ends, respectively.
    cff1 = (grd.y_rho - grd.y_rho[1])*alpha + Hmin
    cff1 += 0.01 * np.random.randn(*grd.y_rho.shape) * cff1
    cff1[0] = cff1[1]
    cff2 = Hmin
    grd.h = np.maximum(cff1, cff2)
    
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(grd.x_psi/1000.0, grd.y_psi/1000.0, '-k', linewidth=0.1)
    ax.plot(grd.x_psi.T/1000.0, grd.y_psi.T/1000.0, '-k', linewidth=0.1)
    hm = np.ma.masked_where(grd.mask==0, grd.h)
    pc = ax.pcolormesh(xvert/1000.0, yvert/1000.0, hm)
    plt.colorbar(pc, orientation='horizontal')
    ax.set_xlim(0, xvert.max()/1000.0)
    ax.set_ylim(0, yvert.max()/1000.0)
    ax.set_aspect(1.0)
    plt.savefig(rootdir + '/grd.png', dpi=300)
    plt.close(fig)
    
    print 'pickling grid..'
    f = open(os.path.join(rootdir, 'grd.pickle'), 'w')
    cPickle.dump(grd, f)
    f.close()
    
    grd.proj = None  # define non-georeference grid
    print 'Writing netcdf GRD file..'
    octant.roms.write_grd(grd, os.path.join(rootdir, 'shelfstrat_grd.nc'), verbose=True)

if __name__ == '__main__':
    from os.path import expanduser
    home = expanduser("~")
    
    # make_grd(rootdir='../project/test_grid')
    make_grd(rootdir=os.path.join(home, 'Projects/shelfstrat/project/test'), shp=(131, 259))