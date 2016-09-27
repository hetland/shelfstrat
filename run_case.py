#!/bin/env python2.7
# encoding: utf-8
"""
run_warner_case.py

Created by Rob Hetland on 2007-10-15.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
Release under MIT license.
"""

import matplotlib
matplotlib.use('Agg')

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
             M20=1e-6, N20=1e-4, f=1e-4, z0=0.003,
             sustr0=0.0, svstr0=0.0,
             shp = (30, 128, 256),
             dt=30.0, ndays=60):
    
    # Make grid
    grdshp = (shp[1]+3, shp[2]+3)
    make_grd(rootdir, Hmin=5.0, alpha=0.001, f=f,
             dx=1e3, dy=1e3, shp=grdshp)
    
    make_frc(rootdir, s_rho=shp[0], sustr0=sustr0, svstr0=svstr0, Tramp=1.0)
    
    print ' RUN CASE sustr = ', sustr0
    
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
    
    rin_3d['NTIMES'] = int(86400 * ndays / dt)
    rin_3d['DT'] = dt
    rin_3d['NDTFAST'] = int(30)
    rin_3d['NHIS'] = int((86400.0 / dt) / 8.0)
    rin_3d['NAVG'] = int((86400.0 / dt) * 3.0)
    rin_3d['NDIA'] = int((86400.0 / dt) * 1.0)
    
    rin_3d['Zob'] = z0
    
    infile = os.path.join(rootdir, 'ocean_shelfstrat.in')
    outfile = os.path.join(rootdir, 'ocean_shelfstrat.out')
    rin_3d.write(infile)
    print infile
    
    print ' ### Running 3D ROMS...'
    os.system('/usr/mpi/gcc/openmpi-1.4.3/bin/mpiexec -np 8 ./project/coawstM %s > %s &' % (infile, outfile))
    # os.system('./project/coawstS < %s' % infile)
    
    return 0
    


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--z0', type=float, default=0.003, 
                        help='Bottom roughness parameter (default=0.003)')
    parser.add_argument('--M2', type=float, default=1e-6, 
                        help='Horizontal stratification parameter (default=1e-6)')
    parser.add_argument('--N2', type=float, default=1e-4, 
                        help='Vertical stratification parameter (default=1e-6)')
    parser.add_argument('--f', type=float, default=1e-4, 
                        help='Coreolis parameter (default=1e-4)')
    parser.add_argument('--sustr', type=float, default=0.0, 
                        help='Along-shore wind stress (default=0.0)')
    parser.add_argument('--svstr', type=float, default=0.0, 
                        help='Along-shore wind stress (default=0.0)')
    parser.add_argument('--Lm', type=int, default=128, 
                        help='Number of x-grid points')
    parser.add_argument('--rootdir', type=str, default='./project/test', 
                        help='Simulation root directory.')
    args = parser.parse_args()
    
    print args
    
    run_case(rootdir=args.rootdir, M20=args.M2, N20=args.N2, f=args.f, sustr0=args.sustr, svstr0=args.svstr, z0=args.z0)
