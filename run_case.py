#!/usr/bin/env python2.7
# encoding: utf-8
"""
run_warner_case.py

Created by Rob Hetland on 2007-10-15.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
Release under MIT license.
"""

from grd import make_grd
# from analysis import make_xsec_movie

import os

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
                    Uriver=0.1, 
                    Utide=1.0,
                    Href=10.0,
                    Zob=0.001,
                    shp=(30, 2, 256),
                    Hmin=3.0, Hmax=15.0, 
                    Lh=100e3, Le=300e3):
    
    # Make grid
    grdshp = (shp[1]+3, shp[2]+3)
    make_grd(rootdir=rootdir, shp=grdshp, Hmin=Hmin, Hmax=Hmax, Le=Le, Lh=Lh)
    
    # run 3D case
    rin_3d = ROMS_in('./project/External/coawst_estuary.in')
    rin_3d['GRDNAME'] = os.path.join(rootdir, 'estuary_grd.nc')
    rin_3d['FRCNAME'] = os.path.join(rootdir, 'estuary_frc.nc')
    rin_3d['HISNAME'] = os.path.join(rootdir, 'estuary_his.nc')
    rin_3d['AVGNAME'] = os.path.join(rootdir, 'estuary_avg.nc')
    rin_3d['DIANAME'] = os.path.join(rootdir, 'estuary_dia.nc')
    rin_3d['ININAME'] = './project/External/estuary_ini.nc'
    rin_3d['RSTNAME'] = os.path.join(rootdir, 'estuary_rst.nc')
    rin_3d['VARNAME'] = './project/External/varinfo.dat'
    rin_3d['NUSER'] = '3'
    rin_3d['USER'] = '%f %f %f' % (Uriver, Utide, Href)
    rin_3d['Lm'] = shp[2]
    rin_3d['Mm'] = shp[1]
    rin_3d['N'] = shp[0]
    
    dt = 30.0
    rin_3d['NTIMES'] = 181400
    rin_3d['DT'] = dt
    rin_3d['NDTFAST'] = 10
    
    rin_3d['NHIS'] = int(86400 / dt / 24 / 4)
    rin_3d['NDIA'] = int(86400 / dt / 24 / 4)
    rin_3d['NDEFHIS'] = int(86400 / dt / 2)
    rin_3d['NDEFDIA'] = int(86400 / dt / 2)
    rin_3d['NAVG'] = int(86400 / dt / 2)
    
    rin_3d['Zob'] = 0.001
    
    infile = os.path.join(rootdir, 'ocean_estuary.in')
    outfile = os.path.join(rootdir, 'ocean_estuary.out')
    rin_3d.write(infile)
    print ' ### Running 3D ROMS...'
    os.system('mpirun -np 8 ./project/coawstS < %s > %s' % (infile, outfile))


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--foo', type=float, default=42.0, help='FOO!')
    args = parser.parse_args()
    print args
