
import numpy as np
import netCDF4
from datetime import datetime
import os

def make_frc(rootdir='../project/', 
             s_rho=30, 
             sustr0=0.0,
             svstr0=0.0,
             Tramp=3.0):
    
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)
    
    grd = np.load(os.path.join(rootdir, 'grd.pickle'))
    
    # one year of forcing, just to be sure.
    t = np.linspace(0, 365, 365*24)
    
    ### write river file
    print 'Writing netcdf FRC file..'
    
    shp = grd.x_rho.shape
    
    nc = netCDF4.Dataset(os.path.join(rootdir, 'shelfstrat_frc.nc'), 'w', format='NETCDF3_CLASSIC')
    nc.Description = 'Surface wind forcing for ideal shelf'
    nc.Author = 'Rob Hetland'
    nc.Created = datetime.now().isoformat()
    nc.type = 'ROMS FRC file'
    
    #########################################
    # River forcing
    
    nc.createDimension('xi_rho', shp[1])
    nc.createDimension('eta_rho', shp[0])
    nc.createDimension('s_rho', s_rho)
    
    def write_nc_var(var, name, dimensions, units=None):
        nc.createVariable(name, 'f8', dimensions)
        if units is not None:
            nc.variables[name].units = units
        nc.variables[name][:] = var
    
    ramp = 1.0 - np.exp(t/Tramp)
    sustr = sustr0 * ramp
    svstr = svstr0 * ramp
    
    nc.createDimension('sms_time', len(t))
    
    write_nc_var(t, 'sms_time', ('sms_time', ), 'days')
    write_nc_var(sustr, 'svstr', ('sms_time', ), 'days')
    write_nc_var(svstr, 'sustr', ('sms_time', ), 'days')
    
    
    nc.close()


if __name__ == '__main__':
    from os.path import expanduser
    home = expanduser("~")

    make_frc(rootdir=os.path.join(home, 'Projects/shelfstrat/project/test'))
