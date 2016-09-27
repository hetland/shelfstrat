
import os, sys
import argparse
import itertools

import numpy as np
import netCDF4

import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory of case to check time')
args = parser.parse_args()

def param_dict(directory):
    rootdir = directory.strip('/').split('/')[-1]
    case = rootdir.split('_')[1:]  # remove leading description 'shelfstrat'
    keys = case[::2]
    vals = [float(val) for val in case[1::2]]
    params = dict(zip(keys, vals))
    params['Ri'] = params['N2'] * params['f']**2 / params['M2']**2
    params['S'] = params['M2'] * 1e-3 * np.sqrt(params['Ri']) / params['f']**2
    params['delta'] = np.sqrt(params['Ri']) * params['S']
    return params

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

hisfilename = os.path.join(args.directory, 'shelfstrat_his.nc')
nc = netCDF4.Dataset(hisfilename)
integration_time = nc.variables['ocean_time'][-1]/86400.0

params = param_dict(args.directory)
omega = np.polyval(omega_poly, params['delta'])
omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri']) # rad/days

timescale = 50 * np.sqrt(params['S']) / omega_dim

if (params['S'] > 0.97) or (params['Ri'] > 11):
    sys.exit()

# print args.directory
# print '%6.2f %5.3f %7.2f %7.2f ' % (params['Ri'], params['S'], timescale, integration_time),
# print ''
if timescale <= integration_time:
     # print ''
    pass
else:
    print args.directory
    print '%6.2f %5.3f %7.2f %7.2f ' % (params['Ri'], params['S'], timescale, integration_time),
    print ''




