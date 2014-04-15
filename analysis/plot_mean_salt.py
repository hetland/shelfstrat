
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

plt.clf()

tidx = 7 * 8

for file in files:
    try:
        nc = netCDF4.Dataset(file)
    except:
        continue

    y = nc.variables['y_rho'][:, 0]/1000.0

    zw = octant.roms.nc_depths(nc, grid='w')[-1]
    dz = np.diff(zw, axis=0)
    
    t = nc.variables['ocean_time'][tidx] / 86400.0
    print file[38:46], ' time = ', t

    salt = nc.variables['salt'][tidx-8:tidx].mean(axis=0)

    sbar = np.sum(salt*dz, axis=0) / np.sum(dz, axis=0)

    sbarmean = sbar.mean(axis=1)
    sini = nc.variables['salt'][0, -1, :, 0]
    
    plt.plot(y, sbarmean - sini, label=file[38:46], color=plt.cm.jet(float(files.index(file))/float(len(files))))


plt.legend()
plt.show()

