
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

for file in files:
    
    hisfilename = os.path.join(args.directory, file, 'shelfstrat_his.nc')
    diafilename = os.path.join(args.directory, file, 'shelfstrat_dia.nc')

    params = cases[file]
    # print hisfilename
    
    omega = np.polyval(omega_poly, params['delta'])
    omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri'])
    
    timescale = timescale_factor / omega_dim
    
    nd = netCDF4.Dataset(diafilename)
    # nc = netCDF4.Dataset(hisfilename)
    
    time = nd.variables['ocean_time'][:]/86400.0
    
    bstr = np.sqrt((nd.variables['ubar_bstr'][:, 1:51, :]**2).mean(axis=-1).mean(axis=-1))
    # ubar = params['M2'] * 50 / params['f'] / 2.0
    # print '%4.2f, %6.3f, %6.2f' % (params['S'], params['Ri'], ubar / bstr / 86400.0)
    
    plt.plot(time/timescale/np.sqrt(params['S']), bstr, '-k')
    
    
# foo =  [[0.10,  1.000,   5.05],
#         [0.07,  1.988,   4.93],
#         [0.06,  2.993,   4.21],
#         [0.04,  5.018,   3.71],
#         [0.03,  9.986,   3.23],
#         [0.30, 10.058,  14.38],
#         [0.50,  3.025,  28.19],
#         [0.50,  2.012,  26.90],
#         [0.30,  4.995,  32.54],
#         [0.20, 10.014,   9.50],
#         [0.30,  3.008,  30.15],
#         [0.50,  1.000,  25.59],
#         [0.35,  2.002,  27.77],
#         [0.29,  2.993,  31.39],
#         [0.22,  4.995,  10.96],
#         [0.16,  9.986,   8.53],
#         [0.20,  4.982,  37.37],
#         [0.30,  1.991,  28.77],
#         [0.20,  2.993,  34.88],
#         [0.10, 10.014,   6.05],
#         [0.10,  9.942,   5.38],
#         [0.30,  1.000,  21.28],
#         [0.21,  2.001,  14.02],
#         [0.17,  3.002,   9.51],
#         [0.13,  5.005,   7.09],
#         [0.20,  1.995,  15.56],
#         [0.10,  5.005,   5.88],
#         [0.09,  5.018,   5.29],
#         [0.06,  9.986,   4.44],
#         [0.20,  1.000,  11.21],
#         [0.14,  1.999,   7.46],
#         [0.12,  3.000,   6.19],
#         [0.10,  3.004,   5.77],
#         [0.50, 10.014,  28.76],
#         [0.10,  2.001,   6.34],
#         [0.50,  5.005,  30.40]]
    
    

