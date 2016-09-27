# README #

### Shelfstrat

commands to run the shelf eddy scenarios. Key files are

    grd/make_grd.py
    ini/make_ini.py
    frc/make_frc.py

These files are run automatically by the script

    run_case.py

This file takes arguments and runs a case:

    usage: run_case.py [-h] [--z0 Z0] [--M2 M2] [--N2 N2] [--f F] [--sustr SUSTR]
                   [--svstr SVSTR] [--Lm LM] [--rootdir ROOTDIR]

    optional arguments:
      -h, --help         show this help message and exit
      --z0 Z0            Bottom roughness parameter (default=0.003)
      --M2 M2            Horizontal stratification parameter (default=1e-6)
      --N2 N2            Vertical stratification parameter (default=1e-6)
      --f F              Coreolis parameter (default=1e-4)
      --sustr SUSTR      Along-shore wind stress (default=0.0)
      --svstr SVSTR      Along-shore wind stress (default=0.0)
      --Lm LM            Number of x-grid points
      --rootdir ROOTDIR  Simulation root directory.

E.g., the 'base' case can be run with the following command:

    mkdir simulations
    ./run_case.py --rootdir simulations

The case can be modified by something like

    ./run_case.py --M2 1e-5 --N2 1e-4 --f 3.33e-3 --rootdir simulations

This script can then be run over a wide parameter space running the script

    run_space.py

Edit this script to create the parameter space you want. This script uses qsub to 
run the various cases, so the qsub script should be edited as well.