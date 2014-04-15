
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

import octant
import octant.roms


files = [
    # 'simulations/shelfstrat_M2_1.00e-06_N2_5.00e-01_f_1.00e-04/shelfstrat_his.nc',
    # 'simulations/shelfstrat_M2_1.00e-06_N2_1.00e-01_f_1.00e-04/shelfstrat_his.nc',
    # 'simulations/shelfstrat_M2_1.00e-06_N2_5.00e-02_f_1.00e-04/shelfstrat_his.nc',
    # 'simulations/shelfstrat_M2_1.00e-06_N2_1.00e-02_f_1.00e-04/shelfstrat_his.nc',
    # 'simulations/shelfstrat_M2_1.00e-06_N2_5.00e-03_f_1.00e-04/shelfstrat_his.nc',
    'simulations/shelfstrat_M2_1.00e-06_N2_1.00e-03_f_1.00e-04/shelfstrat_his.nc',
    'simulations/shelfstrat_M2_1.00e-06_N2_1.50e-04_f_1.00e-04/shelfstrat_his.nc',
    'simulations/shelfstrat_M2_1.00e-06_N2_1.00e-04_f_1.00e-04/shelfstrat_his.nc',
    'simulations/shelfstrat_M2_1.00e-06_N2_5.00e-05_f_1.00e-04/shelfstrat_his.nc', 
    'simulations/shelfstrat_M2_1.00e-06_N2_1.00e-07_f_1.00e-04/shelfstrat_his.nc',
    'simulations/shelfstrat_M2_1.00e-06_N2_1.00e-08_f_1.00e-04/shelfstrat_his.nc',
    'simulations/shelfstrat_M2_1.00e-06_N2_0.00e+00_f_1.00e-04/shelfstrat_his.nc',
]


tidx = 8 * 7

bins_M2 = np.linspace(-8, -3, 101)
bins_M2_c = 0.5 * (bins_M2[1:] + bins_M2[:-1])

bins_N2 = np.linspace(-8, -1, 301)
bins_N2_c = 0.5 * (bins_N2[1:] + bins_N2[:-1])

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

frame = 0

for tidx in range(56)[::8]:

    ax1.cla()
    ax2.cla()

    for file in files:
        nc = netCDF4.Dataset(file)
    
        pm = nc.variables['pm'][:]
        pn = nc.variables['pn'][:]
    
        zr = octant.roms.nc_depths(nc, grid='rho')[tidx]
        
        rho = nc.variables['rho'][tidx]
        M2x, M2y = octant.tools.hgrad(rho, zr, pm, pn)
        N2 = octant.tools.N2(rho, zr, rho_0=nc.variables['rho0'][0])
    
        shp = (N2.shape[0], N2.shape[1]-2, N2.shape[2]-2)
        M2x = octant.tools.shrink(M2x, shp)
        M2y = octant.tools.shrink(M2y, shp)
        N2 = octant.tools.shrink(N2, shp)
    
        M2 = np.sqrt(M2x**2 + M2y**2)
    
        dz = np.diff(zr, axis=0)
        dV = dz / (pm * pn)
        dV = octant.tools.shrink(dV, shp)
    
        M2hist = np.histogram(np.log10(M2), bins=bins_M2, weights=dV)[0]
        N2hist = np.histogram(np.log10(N2), bins=bins_N2, weights=dV)[0]
    
        color=plt.cm.jet( float(files.index(file))/float(len(files)) )
    
        ax1.plot(bins_M2_c, M2hist, color=color)
        ax1.set_ylim(0, 1e11)
        
        ax2.plot(bins_N2_c, N2hist, color=color)
        ax1.set_ylim(0, 5e11)
        
        
    plt.savefig('foo_%04d.png' % frame)
    frame += 1


plt.show()
    

