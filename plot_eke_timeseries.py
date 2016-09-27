
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import netCDF4

import octant
import octant.roms

from case_dictionary import case_dictionary

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory with top-level cases listed')
args = parser.parse_args()

cases = case_dictionary(args.directory)

##### USE THE UNION.

files_delta = [cases.find(Ri=Ri, delta=delta, N2=1e-4)
                for delta in [0.1, 0.2, 0.3, 0.5]
                for Ri in [1, 2, 3, 5, 10]]

files_S = [cases.find(Ri=Ri, S=S, N2=1e-4)
            for S in [0.1, 0.2, 0.3, 0.5]
            for Ri in [1, 2, 3, 5, 10]]

files = list(np.unique(np.asarray(files_S + files_delta)))

def polyfit_omega(n=6):
    'fit an order-n polynomial to the maximum growth rate as a function of delta, the slope parameter.'
    delta, mu = np.mgrid[-1.2:2.2:1001j, 0:4.2:1001j]
    
    tmu = np.tanh(mu)
    omega = np.sqrt( (1.0+delta)*(mu - tmu)/tmu 
                     -0.25*(delta/tmu + mu)**2 ).real
    omega2 = (1.0+delta)*(mu - tmu)/tmu - 0.25*(delta/tmu + mu)**2
    
    omega = np.ma.masked_where(np.isnan(omega), omega)
    
    omega_max = omega.max(axis=1)
    idx = np.where(~omega_max.mask)
    
    omega_max = omega_max[idx]
    delta = delta[:, 0][idx]
    
    p = np.polyfit(delta, omega_max, n)
    
    return p

omega_poly = polyfit_omega()

ref_timescale = 10.0 # days
ref_delta = 0.1
ref_Ri = 1.0
ref_f = 1e-4

omega = np.polyval(omega_poly, ref_delta) # non-dim
omega_dim = 86400.0 * omega * ref_f / np.sqrt(ref_Ri)  # rad/days
timescale_factor = ref_timescale * omega_dim


Ris = []
deltas = []
tkes = []
ekes = []
mkes = []
epes = []
mpes = []
tpes = []
times = []
timescales = []
omegas = []

for file in files:
    
    hisfilename = os.path.join(args.directory, file, 'shelfstrat_his.nc')
    params = cases[file]
    print hisfilename
    
    omega = np.polyval(omega_poly, params['delta'])
    omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri'])
    
    timescale = timescale_factor / omega_dim
    
    nc = netCDF4.Dataset(hisfilename)
    time = nc.variables['ocean_time'][:] / 86400.0
    
    pm = nc.variables['pm'][:]
    pn = nc.variables['pn'][:]
    dA = 1.0/(pm*pn)
    
    zw = octant.roms.nc_depths(nc, 'w')
    zr = octant.roms.nc_depths(nc, 'rho')
    
    tke = np.zeros(len(time))
    eke = np.zeros(len(time))
    mke = np.zeros(len(time))
    epe = np.zeros(len(time))
    mpe = np.zeros(len(time))
    tpe = np.zeros(len(time))
    
    R0 = nc.variables['R0'][0]
    Scoef = nc.variables['Scoef'][0]
    
    for n in range(len(time)):
        u = nc.variables['u'][n, :, 1:-1, :]
        v = nc.variables['v'][n, :, :, 1:-1]
        u, v = octant.tools.shrink(u, v)

        umean = u.mean(axis=-1)[..., None]
        up = u - umean
        vp = v

        dV = (dA*np.diff(zw[n], axis=0))[:, 1:-1, 1:-1]
        sum_dV = np.sum(dV)
        
        tke[n] = R0*np.sum(0.5*(u**2 + v**2)*dV) / sum_dV
        eke[n] = R0*np.sum(0.5*(up**2 + vp**2)*dV) / sum_dV
        mke[n] = R0*np.sum(0.5*(umean**2)*dV) / sum_dV
        
        salt = nc.variables['salt'][n, :, 1:-1, 1:-1]
        rho_s = R0*(Scoef*(salt-35.0))
        rho_s_mean = rho_s.mean(axis=-1)[..., None]
        rho_s_prime = rho_s - rho_s_mean
        
        z = zr[n][:, 1:-1, 1:-1]
        epe[n] = np.sum(rho_s_prime*9.8*z*dV) / sum_dV
        mpe[n] = np.sum(rho_s_mean*9.8*z*dV) / sum_dV
        tpe[n] = np.sum(rho_s*9.8*z*dV) / sum_dV
        
        
    Ris.append(params['Ri'])
    deltas.append(params['delta'])
    tkes.append(tke)
    ekes.append(eke)
    mkes.append(mke)

    epes.append(epe)
    mpes.append(mpe)
    tpes.append(tpe)
    
    times.append(time)
    timescales.append(timescale)
    omegas.append(omega_dim)


np.save('eke_Ris', Ris)
np.save('eke_deltas', deltas)
np.save('eke_tkes', tkes)
np.save('eke_ekes', ekes)
np.save('eke_mkes', mkes)
np.save('eke_epes', epes)
np.save('eke_mpes', mpes)
np.save('eke_tpes', tpes)
np.save('eke_times', times)
np.save('eke_timescales', timescales)
np.save('eke_omegas', omegas)
    
    

