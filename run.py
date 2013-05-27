#!/usr/local/bin/env python
# encoding: utf-8
"""
run.py

Created by Rob Hetland on 2007-10-16.
Copyright (c) 2007 Texas A&M Univsersity. All rights reserved.
"""

import sys
import getopt

import os
from run_warner_case import run_warner_case 

help_message = '''
The help message goes here.
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hv", 
                         ["help", "Uriver=", "Utide=", "Hmin=", "Hmax=", "Le=", "Lh=", "Zob="])
        except getopt.error, msg:
            raise Usage(msg)
    
        optiondict = {"Hmin": 3.0, "Hmax": 15.0, "Zob": 0.001, 
                      "Le": 300e3, "Lh": 100e3}
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("--Uriver"):
                optiondict['Uriver'] = float(value)
            if option in ("--Utide"):
                optiondict['Utide'] = float(value)
            if option in ("--Hmin"):
                optiondict['Hmin'] = float(value)
            if option in ("--Hmax"):
                optiondict['Hmax'] = float(value)
            if option in ("--Le"):
                optiondict['Le'] = float(value)
            if option in ("--Lh"):
                optiondict['Lh'] = float(value)
            if option in ("--Zob"):
                optiondict['Zob'] = float(value)
    
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    
    rootdir = './archive/warner_Uriver%4.2f_Utide%4.2f_Zob%4.3f'  \
              % (optiondict['Uriver'],
                 optiondict['Utide'],
                 optiondict['Zob'])
    optiondict['rootdir'] = rootdir
    
    print 'Runing case:'
    for key, val in optiondict.iteritems():
        print '    ', key, val
    
    run_warner_case(**optiondict) 


if __name__ == "__main__":
    sys.exit(main())
