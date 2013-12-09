
import numpy as np
import netCDF4
from datetime import datetime
import os

import octant

def make_ini(rootdir='../project/', 
             s_rho=30, theta_s = 3.0, theta_b = 0.4, hc = 5.0,
             Vtransform=1, Vstretching=1,
             R0=1027.0, T0=25.0, S0=35.0, TCOEF=1.7e-4, SCOEF=7.6e-4,
             M20=1e-7, M2_yo=50e3, M2_r=5e3,
             N20=1e-4, N2_zo=50.0, N2_r=50.0):
    '''
    Create an initialization file.
    
    Horizontal stratification is controlled by salinity (only)
    Vertical stratification is controlled by temperature (only)
    
    Stratification properties are conservative through 
    a linear equation of state:
    
       R0 == 1027.0d0                   ! kg/m3
       T0 == 25.0d0                      ! Celsius
       S0 == 35.0d0                     ! PSU
    TCOEF == 1.7d-4                     ! 1/Celsius
    SCOEF == 7.6d-4                     ! 1/PSU
    
    
    '''
    
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)    
    grd = np.load(os.path.join(rootdir, 'grd.pickle'))
    
    alpha = TCOEF
    beta = SCOEF
     
    g = 9.8
    
    dy = grd.dy[0,0]   # assume constant grid spacing.
    
    M2 = M20 * np.exp( (M2_yo - grd.y_rho) / M2_r )
    M2[grd.y_rho < M2_yo] = M20
    s = np.cumsum(M2 * dy / (g * beta), axis=0)
    s -= s[-1] - S0
    s = s * np.ones((s_rho, 1, 1), 'd')
    print 'Coastal salinity: ', s.min(), ' for M2=',M20
    
    z = octant.depths.get_zrho(Vtransform, Vstretching, s_rho, theta_s, theta_b, grd.h, hc)
    Hz = octant.depths.get_Hz(Vtransform, Vstretching, s_rho, theta_s, theta_b, grd.h, hc)
    
    N2 = N20 * np.exp( -(N2_zo - z) / N2_r )
    N2[z < N2_zo] = N20
    
    t = np.zeros_like(s)
    for n in range(s_rho):
        t[n] = T0 - np.trapz(N2[n:] / (g * alpha), x=z[n:], axis=0)
    
    
    #########################################
    # Create NetCDF file
    
    nc = netCDF4.Dataset(os.path.join(rootdir, 'shelfstrat_ini.nc'), 'w', format='NETCDF3_CLASSIC')
    nc.Description = 'Initial conditions for ideal shelf'
    nc.Author = 'Rob Hetland'
    nc.Created = datetime.now().isoformat()
    nc.type = 'ROMS FRC file'
    
    shp = grd.x_rho.shape
    nc.createDimension('xi_rho', shp[1])
    nc.createDimension('eta_rho', shp[0])
    nc.createDimension('s_rho', s_rho)
    nc.createDimension('xi_u', shp[1]-1)
    nc.createDimension('eta_u', shp[0])
    nc.createDimension('xi_v', shp[1])
    nc.createDimension('eta_v', shp[0]-1)
    
    nc.createDimension('ocean_time', 1)
    nc.createVariable('ocean_time', 'd', ('ocean_time',))
    nc.variables['ocean_time'][0] = 0.0
    nc.variables['ocean_time'].units = 'days'
    
    def write_nc_var(var, name, dimensions, units=None):
        nc.createVariable(name, 'f8', dimensions)
        if units is not None:
            nc.variables[name].units = units
        nc.variables[name][:] = var
    
    # variable_list = ('salt', 'temp', 'u', 'v', 'zeta', 'ubar', 'vbar')
    
    write_nc_var(s, 'salt', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), 'PSU')
    write_nc_var(t, 'temp', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), 'degC')
    
    write_nc_var(0.0, 'u', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), 'm s-1')
    write_nc_var(0.0, 'v', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), 'm s-1')
    
    write_nc_var(0.0, 'zeta', ('ocean_time', 'eta_rho', 'xi_rho'), 'm')
    write_nc_var(0.0, 'ubar', ('ocean_time', 'eta_u', 'xi_u'), 'm s-1')
    write_nc_var(0.0, 'vbar', ('ocean_time', 'eta_v', 'xi_v'), 'm s-1')
    
    nc.close()


if __name__ == '__main__':
    from os.path import expanduser
    home = expanduser("~")

    make_ini(rootdir=os.path.join(home, 'Projects/shelfstrat/project/test'))
    