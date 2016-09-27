
import os, sys
import argparse
import itertools

import numpy as np
import netCDF4

# import octant

import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory of case to extend')
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

# ref_timescale = 6.0 # days
# ref_delta = 0.1
# ref_Ri = 1.0
# ref_f = 1e-4
# omega = np.polyval(omega_poly, ref_delta) # non-dim
# omega_dim = 86400.0 * omega * ref_f / np.sqrt(ref_Ri)  # rad/days
# timescale_factor = ref_timescale * omega_dim

hisfilename = os.path.join(args.directory, 'shelfstrat_his.nc')
nc = netCDF4.Dataset(hisfilename)
integration_time = nc.variables['ocean_time'][-1]/86400.0

params = param_dict(args.directory)
omega = np.polyval(omega_poly, params['delta'])
omega_dim = 86400.0 * omega * params['f'] / np.sqrt(params['Ri']) # rad/days

# Old way, using just omega to estimate integration time.
# timescale = timescale_factor / omega_dim

# Here, use the sqrt(S) factor, and integrate to t omega / sqrt(S) == 50
timescale = 50 * np.sqrt(params['S']) / omega_dim

print timescale, integration_time

if integration_time >= timescale:
    print ' ### %s EXITING -- already there' % args.directory
    sys.exit()

preamble = '''#!/bin/bash

### Set the job name
#PBS -N shelfstrat

### Declare myprogram non-rerunable
#PBS -r n

### Set the queue to "default"
#PBS -q default

### Run on
#PBS -l nodes=1:ppn=8

### You can override the default 1 hour real-world time limit.
#PBS -l walltime=172:00:00

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo HOSTNAME: $HOSTNAME
/usr/mpi/gcc/openmpi-1.4.3/bin/mpiexec -np 8 ./project/coawstM %s/ocean_extend.in > %s/ocean_extend.out
''' % (args.directory.strip('/'), args.directory.strip('/'))


class ROMS_in(object):
    """docstring for ROMS_IN"""
    def __init__(self, infile):
        self.infile = infile
        f = open(self.infile)

        self.variables = {}
        self._varlist = []       # to keep the order of variables

        for line in f.readlines():
            # strip comments
            s = line.split('!')[0].strip()
            if len(s) == 0: continue     # comment only
            vals = s.split('=')
            self.variables[vals[0].strip()] = vals[-1].strip()
            self._varlist.append(vals[0].strip())

    def write(self, filename):
        """docstring for write"""
        f = open(filename, 'w')
        for key in self._varlist:
            f.write('%s == %s\n' % (key, str(self.variables[key])))
        f.close()

    def __setitem__(self, key, val):
        """docstring for __setitem__"""
        if key not in self._varlist:
            self._varlist.append(key)
            warnings.warn('%s not previously in variable list.' % key)
        self.variables[key] = str(val)


rin = ROMS_in('%s/ocean_shelfstrat.in' % args.directory.strip('/'))

rin.variables['ININAME'] = './%s/shelfstrat_his.nc' % args.directory.strip('/')

ntimes = np.ceil(timescale) * 86400.0 / float(rin.variables['DT'])
print ' ### from NTIMES=',rin.variables['NTIMES'], ' to NTIMES=', str(int(ntimes))

rin.variables['NTIMES'] = str(int(ntimes))

rin.variables['LDEFOUT'] = 'F'
rin.variables['NRREC'] = '-1'

rin.write('%s/ocean_extend.in' % args.directory.strip('/'))

casestr = args.directory.split('shelfstrat')[-1].strip('/')


qfile = open('shelfstrat%s.sh' % casestr, 'w')
qfile.writelines(preamble)
qfile.close()

os.system('qsub shelfstrat%s.sh' % casestr)



