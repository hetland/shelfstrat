
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import netCDF4

import octant

from case_dictionary import case_dictionary

# supress warnings that occur in division in the polyfit_omega function
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory with top-level cases listed')
args = parser.parse_args()

cases = case_dictionary(args.directory)


# files = [cases.find(Ri=Ri, delta=delta, N2=1e-4)
#             for delta in [0.1, 0.2, 0.3, 0.5]
#             for Ri in [1, 2, 3, 5, 10]]
# fig, axs = plt.subplots(4, 5, figsize=(13.5, 6))
# normalize_time = False


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


times = []
timescales = []
stresses = []
ubar = []
S = []
Ri = []

for file in files:
    print file
    
    diafilename = os.path.join(args.directory, file, 'shelfstrat_dia.nc')
    params = cases[file]
    
    
    nd = netCDF4.Dataset(diafilename)
    time = nd.variables['ocean_time'][:]/86400.0
    
    ubar_bstr = nd.variables['ubar_bstr'][:, 1:51, :]
    vbar_bstr = nd.variables['vbar_bstr'][:, 1:50, :] 
    
    # ubar_bstr, vbar_bstr = octant.tools.shrink(ubar_bstr, vbar_bstr)
    # stress_mag = np.sqrt(ubar_bstr**2 + vbar_bstr**2)

    stress_mag = np.sqrt(ubar_bstr**2)
    
    stress = np.sqrt((stress_mag**2).mean(axis=-1).mean(axis=-1))
    
    omega = np.polyval(omega_poly, params['delta'])
    omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri']) # rad/days
    
    timescale = np.sqrt(params['S']) / omega_dim
    
    ubar.append( params['Ri']**(-1./2.) * np.sqrt(1e-4) * 50.0 / 2.0 )
    times.append(time)
    timescales.append(timescale)
    stresses.append(stress)
    S.append(params['S'])
    Ri.append(params['Ri'])



# frictional timescales
ubar = np.array(ubar)
tau = np.array([_[1] for _ in stresses])
tau = np.array(tau)
Tf = ubar/tau/86400.0

# Intstability timescales
Ti = np.array(timescales)

np.save('friction_timescales', vstack((Ti, Tf, S, Ri)))
