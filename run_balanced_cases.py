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
    
def qsub_case(f, M2, N2):
    suffix = get_suffix(f=f, M2=M2, N2=N2)
    rootdir = './simulations/shelfstrat%s' % suffix
    
    Ri = N2 * f**2 / M2**2
    phi = M2 * 1e-3 * np.sqrt(Ri) / f**2
    
    print 'M2 = %4.2e;   N2 = %4.2e;    f = %4.2e;    Ri=%4.2f;   phi=%4.2f' % (M2, N2, f, Ri, phi)
    qfile = open('shelfstrat%s.sh' % suffix, 'w')
    qfile.writelines(preamble)
    qfile.write('/anaconda/bin/python2.7 ./run_case.py --rootdir %s --M2 %4.2e --N2 %4.2e --f %f \n' % (rootdir, M2, N2, f))
    qfile.write('/usr/mpi/gcc/openmpi-1.4.3/bin/mpiexec -np 8 ./project/coawstM %s/ocean_shelfstrat.in > %s/ocean_shelfstrat.out\n' % (rootdir, rootdir))
    qfile.close()
    os.system('qsub shelfstrat%s.sh' % suffix)




qsub_case(M2=1.00e-06, N2=1.00e-04, f=1.00e-04)
qsub_case(M2=5.00e-07, N2=1.00e-04, f=5.00e-05)
qsub_case(M2=2.00e-07, N2=1.00e-04, f=2.00e-05)
qsub_case(M2=1.00e-07, N2=1.00e-04, f=1.00e-05)
qsub_case(M2=7.07e-07, N2=1.00e-04, f=1.00e-04)
qsub_case(M2=3.54e-07, N2=1.00e-04, f=5.00e-05)
qsub_case(M2=1.41e-07, N2=1.00e-04, f=2.00e-05)
qsub_case(M2=7.07e-08, N2=1.00e-04, f=1.00e-05)
qsub_case(M2=4.47e-07, N2=1.00e-04, f=1.00e-04)
qsub_case(M2=2.24e-07, N2=1.00e-04, f=5.00e-05)
qsub_case(M2=8.94e-08, N2=1.00e-04, f=2.00e-05)
qsub_case(M2=4.47e-08, N2=1.00e-04, f=1.00e-05)
qsub_case(M2=3.16e-07, N2=1.00e-04, f=1.00e-04)
qsub_case(M2=1.58e-07, N2=1.00e-04, f=5.00e-05)
qsub_case(M2=6.32e-08, N2=1.00e-04, f=2.00e-05)
qsub_case(M2=3.16e-08, N2=1.00e-04, f=1.00e-05)

