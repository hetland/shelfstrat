#!/bin/env python2.7

import os
import argparse
import itertools

import numpy as np

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
#PBS -l walltime=72:00:00

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR 

echo HOSTNAME: $HOSTNAME
'''


def get_suffix(**kwargs):
    keys = kwargs.keys()
    keys.sort()
    suffix=''
    for key in keys:
        suffix += '_%s_%4.2e' % (str(key), kwargs[key])
    return suffix
    
def qsub_case(sustr=0.0, svstr=0.0, dryrun=False):
    suffix = get_suffix(svstr=svstr)
    rootdir = './sims_cross/shelfstrat%s' % suffix
    
    print suffix
    M2 = 1e-6
    N2 = 1e-4
    f = 1e-4
    
    Ri = N2 * f**2 / M2**2
    phi = M2 * 1e-3 * np.sqrt(Ri) / f**2
    delta = np.sqrt(Ri) * phi
    
    print 'M2 = %4.2e;   N2 = %4.2e;    f = %4.2e;    Ri=%4.2f;   phi=%4.2f,   sustr=%0.6f,   sustr=%0.6f' % (M2, N2, f, Ri, phi, sustr, svstr)
    
    if os.path.exists(rootdir):
        print ' ### DIRECTORY ', rootdir, ' EXISTS'
        return
    
    if not dryrun:
        qfile = open('shelfstrat%s.sh' % suffix, 'w')
        qfile.writelines(preamble)
        qfile.write('/anaconda/bin/python2.7 ./run_case.py --rootdir %s --M2 %4.2e --N2 %4.2e --f %f --sustr %0.6f --svstr %0.6f\n' % (rootdir, M2, N2, f, sustr, svstr))
        qfile.write('/usr/mpi/gcc/openmpi-1.4.3/bin/mpiexec -np 8 ./project/coawstM %s/ocean_shelfstrat.in > %s/ocean_shelfstrat.out\n' % (rootdir, rootdir))
        qfile.close()
        os.system('qsub shelfstrat%s.sh' % suffix)


for wind_speed in [-40, -30, -20, -10, 10, 20, 30, 40]:
    svstr = 1e-3 * wind_speed * abs(wind_speed)
    qsub_case(svstr=svstr)
