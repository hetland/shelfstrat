
import os
import argparse
import itertools

import numpy as np

import octant


parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='directory with top-level cases listed')
parser.add_argument('--factor', type=float, default=0.0,
                        help='Factor by which to extend runs')
parser.add_argument('--days', type=float, default=0.0,
                        help='Factor by which to extend runs')
args = parser.parse_args()




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

rin.variables['ININAME'] = '%s/shelfstrat_his.nc' % args.directory.strip('/')

if args.factor > 0:
    ntimes = args.factor * int(rin.variables['NTIMES'])
else:
    ntimes = args.days * 86400.0 / float(rin.variables['DT'])
    print ' ### from NTIMES=',rin.variables['NTIMES'], ' to NTIMES=', ntimes

rin.variables['NTIMES'] = str(int(ntimes))
    
rin.variables['LDEFOUT'] = 'F'
rin.variables['NRREC'] = '-1'

rin.write('%s/ocean_extend.in' % args.directory.strip('/'))

casestr = args.directory.split('shelfstrat')[-1].strip('/')


qfile = open('shelfstrat%s.sh' % casestr, 'w')
qfile.writelines(preamble)
qfile.close()
os.system('qsub shelfstrat%s.sh' % casestr)



