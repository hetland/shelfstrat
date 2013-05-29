#!/bin/env python2.7
# encoding: utf-8
"""
run_warner_case.py

Created by Rob Hetland on 2007-10-15.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
Release under MIT license.
"""

import os

from grd import make_grd
from frc import make_frc
from ini import make_ini


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


def run_case(rootdir='./project',
             M20=1e-6, N20=1e-4,
             shp = (30, 128, 256)):
    
    # Make grid
    grdshp = (shp[1]+3, shp[2]+3)
    make_grd(rootdir, Hmin=5.0, alpha=0.001, f=1e-4,
             dx=1e3, dy=1e3, shp=grdshp)
    
    make_frc(rootdir, s_rho=shp[0], sustr0=0.0, svstr0=0.0, Tramp=3.0)
    
    make_ini(rootdir, s_rho=30, theta_s = 3.0, theta_b = 0.4, hc = 5.0,
             Vtransform=1, Vstretching=1,
             R0=1027.0, T0=25.0, S0=35.0, TCOEF=1.7e-4, SCOEF=7.6e-4,
             M20=M20, M2_yo=50e3, M2_r=5e3,
             N20=N20, N2_zo=50.0, N2_r=50.0)
    
    # run 3D case
    rin_3d = ROMS_in('./project/ocean_shelfstrat.in')
    rin_3d['GRDNAME'] = os.path.join(rootdir, 'shelfstrat_grd.nc')
    rin_3d['FRCNAME'] = os.path.join(rootdir, 'shelfstrat_frc.nc')
    rin_3d['HISNAME'] = os.path.join(rootdir, 'shelfstrat_his.nc')
    rin_3d['AVGNAME'] = os.path.join(rootdir, 'shelfstrat_avg.nc')
    rin_3d['DIANAME'] = os.path.join(rootdir, 'shelfstrat_dia.nc')
    rin_3d['ININAME'] = os.path.join(rootdir, 'shelfstrat_ini.nc')
    rin_3d['RSTNAME'] = os.path.join(rootdir, 'shelfstrat_rst.nc')
    rin_3d['VARNAME'] = './project/varinfo.dat'

    rin_3d['Lm'] = shp[2]
    rin_3d['Mm'] = shp[1]
    rin_3d['N'] = shp[0]
    
    infile = os.path.join(rootdir, 'ocean_shelfstrat.in')
    outfile = os.path.join(rootdir, 'ocean_shelfstrat.out')
    rin_3d.write(infile)
    # os.system('/usr/mpi/gcc/openmpi-1.4.3/bin/mpiexec -np 8 ./project/coawstM %s' % infile)
    # os.system('./project/coawstS < %s' % infile)
    
    print ' ### Running 3D ROMS...'
    return 0
    


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--M2', type=float, default=1e-6, 
                        help='Horizontal stratification parameter (default=1e-6)')
    parser.add_argument('--N2', type=float, default=1e-4, 
                        help='Vertical stratification parameter (default=1e-6)')
    parser.add_argument('--rootdir', type=str, default='./project/test', 
                        help='Simulation root directory.')
    args = parser.parse_args()
    
    run_case(rootdir=args.rootdir, M20=args.M2, N20=args.N2)
