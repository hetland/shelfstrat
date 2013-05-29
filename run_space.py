#!/bin/env python2.7

import os
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
    
    
for M2 in np.logspace(-6, -8, 3):
    for N2 in np.logspace(-4, -6, 3):
        suffix = get_suffix(M2=M2, N2=N2)
        rootdir = './simulations/shelfstrat%s' % suffix
        f = open('shelfstrat%s.sh' % suffix, 'w')
        f.writelines(preamble)
        f.write('python2.7 ./run_case.py --rootdir %s --M2 %4.2e --N2 %4.2e\n' % (rootdir, M2, N2))
        f.write('/usr/mpi/gcc/openmpi-1.4.3/bin/mpiexec -np 8 ./project/coawstM %s/ocean_shelfstrat.in > %s/ocean_shelfstrat.out\n' % (rootdir, rootdir))
        f.close()
        os.system('qsub shelfstrat%s.sh' % suffix)



