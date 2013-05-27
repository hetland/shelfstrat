
import numpy as np
import netCDF4
from datetime import datetime
import os

def make_frc(rootdir='../project/', 
             s_rho=30, 
             Qbar_river=10000.0):
    
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)
    
    grd = np.load(os.path.join(rootdir, 'grd.pickle'))
    
    Nsrc = 5
    river_temp = 25.0
    river_Eposition = 1
    
    # one year of forcing, just to be sure.
    t = np.linspace(0, 365, 365*24)
    Qbar = Qbar_river*(1 - np.exp(-t/0.3))
    
    ### write river file
    print 'Writing netcdf FRC file..'
    
    shp = grd.x_rho.shape
    
    nc = netCDF4.Dataset(os.path.join(rootdir, 'shelfplume_frc.nc'), 'w', format='NETCDF3_CLASSIC')
    nc.Description = 'River forcing for ideal plume'
    nc.Author = 'Rob Hetland'
    nc.Created = datetime.now().isoformat()
    nc.type = 'ROMS FRC file'
    
    #########################################
    # River forcing
    
    nc.createDimension('river', Nsrc)
    nc.createDimension('river_time', len(t))
    nc.createDimension('xi_rho', shp[1])
    nc.createDimension('eta_rho', shp[0])
    nc.createDimension('s_rho', s_rho)
    
    def write_nc_var(var, name, dimensions, units=None):
        nc.createVariable(name, 'f8', dimensions)
        if units is not None:
            nc.variables[name].units = units
        nc.variables[name][:] = var
    
    write_nc_var(np.arange(Nsrc)+1, 'river', ('river', ))
    write_nc_var(np.arange(Nsrc)+30, 'river_Xposition', ('river', ))
    write_nc_var(np.ones(Nsrc), 'river_Eposition', ('river', ))
    write_nc_var(np.ones(Nsrc), 'river_direction', ('river', ))
    write_nc_var(np.ones((s_rho, Nsrc))/float(s_rho), 'river_Vshape', ('s_rho', 'river'))
    write_nc_var(t, 'river_time', ('river_time', ), 'days')
    write_nc_var(Qbar[:, np.newaxis]*np.ones((len(t), Nsrc))/float(Nsrc), 'river_transport',
                 ('river_time', 'river'), 'm3 s-1')
    write_nc_var(river_temp*np.ones((len(t), s_rho, Nsrc)), 'river_temp',
                 ('river_time', 's_rho', 'river'))
    write_nc_var(np.zeros((len(t), s_rho, Nsrc)), 'river_salt', ('river_time', 's_rho', 'river'))
    write_nc_var(3.0*np.ones(Nsrc), 'river_flag', ('river', ))
    
    
    ###########################################
    # wind forcing
    
    t = np.linspace(0, 365, 365*24)
    sustr = 5.0e-4 * np.ones_like(t)
    svstr = -5.0e-4 * np.ones_like(t)
    
    nc.createDimension('sms_time', len(t))
    
    write_nc_var(t, 'sms_time', ('sms_time', ), 'days')
    write_nc_var(sustr, 'svstr', ('sms_time', ), 'days')
    write_nc_var(svstr, 'sustr', ('sms_time', ), 'days')
    
    
    nc.close()


if __name__ == '__main__':
    from os.path import expanduser
    home = expanduser("~")

    make_frc(rootdir=os.path.join(home, 'Projects/shelfplume/project/test'))
