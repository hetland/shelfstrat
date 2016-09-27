
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import netCDF4

import octant

from case_dictionary import case_dictionary

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory with top-level cases listed')
args = parser.parse_args()

cases = case_dictionary(args.directory)

##### TRY BOTH FILES SETS TO EXAMINE DIFFERENCE. MAYBE USE THE UNION.

files_delta = [cases.find(Ri=Ri, delta=delta, N2=1e-4)
                for delta in [0.1, 0.2, 0.3, 0.5]
                for Ri in [1, 2, 3, 5, 10]]

files_S = [cases.find(Ri=Ri, S=S, N2=1e-4)
            for S in [0.1, 0.2, 0.3, 0.5]
            for Ri in [1, 2, 3, 5, 10]]

files = list(np.unique(files_S + files_delta))
files = [[file] for file in files]

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
ekes = []
mkes = []
tkes = []
norms = []
len_flags = []
timescales = []
times = []

for file in files:
    
    hisfilename = os.path.join(args.directory, file[0], 'shelfstrat_his.nc')
    params = cases[file[0]]
    print hisfilename
    
    omega = np.polyval(omega_poly, params['delta'])
    omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri'])
    
    timescale = timescale_factor / omega_dim
    timescales.append(timescale)
    
    nc = netCDF4.Dataset(hisfilename)
    time = nc.variables['ocean_time'][:] / 86400.0

    tidx = np.where( time >= timescale )[0]
    if len(tidx) == 0:
        tidx = len(time) - 1
        len_flags.append(0)
    else:
        tidx = tidx.min()
        len_flags.append(1)

    print('   time index: {0:d}/{1:d} -- {2:f}'.format(tidx, len(time), time[tidx]))
    
    u = nc.variables['u'][:tidx, -1, :, :]
    v = nc.variables['v'][:tidx, -1, :, :]
    time = time[:tidx]
    times.append(time[-1])
    
    u, v = octant.tools.shrink(u, v)
    
    umean = u.mean(axis=-1)[:, :, None]
    up = u - umean
    vp = v
    
    tke = 0.5*(u**2 + v**2)
    eke = 0.5*(up**2 + v**2)
    mke = 0.5*(umean**2)
    
    norm = mke.mean(axis=-1).mean(axis=-1)[0]
        
    Ris.append(params['Ri'])
    deltas.append(params['delta'])

    ekes.append( (eke.mean(axis=-1).mean(axis=-1)/norm).max() )
    mkes.append( (mke.mean(axis=-1).mean(axis=-1)).max() )
    tkes.append( (tke.mean(axis=-1).mean(axis=-1)).max() )
    
    # #  energy at days == ref_timescale
    # ekes.append( (eke.mean(axis=-1).mean(axis=-1)/norm)[tidx] )
    # mkes.append( (mke.mean(axis=-1).mean(axis=-1))[tidx] )
    # tkes.append( (tke.mean(axis=-1).mean(axis=-1))[tidx] )
    
    norms.append( norm )
    

Ri = np.array(Ris)
delta = np.array(deltas)
eke = np.array(ekes)
S = delta/np.sqrt(Ri)
len_flag = np.array(len_flags)
timescale = np.array(timescales)
time = np.array(times)

np.save('Ri', Ri)
np.save('S', S)
np.save('eke', eke)
np.save('len_flag', len_flag)
np.save('timescale', timescale)
np.save('time', time)


idx_time = (time/timescale) > 0.95

idx = np.isclose(S, 0.1, rtol=1e-2)
idx = idx | np.isclose(S, 0.2, rtol=1e-2)
idx = idx | np.isclose(S, 0.3, rtol=1e-2)
idx = idx | np.isclose(S, 0.5, rtol=1e-2)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.semilogx(S, eke, 'k.')
ax.semilogx(S[idx_time], eke[idx_time], 'ro')

p = np.polyfit(np.log(S[idx_time]), eke[idx_time], 1)

logS = np.linspace(-5, 0, 50)
ax.plot(np.exp(logS), p[0]*logS + p[1], 'k--')
ax.set_xlim(1e-2, 1)
ax.set_ylim(0, 2)
r = np.corrcoef(np.log(S), eke)[0, 1]
ax.set_xlabel('S')
ax.set_ylabel('Normalized eddy kinetic energy')

plt.show()
