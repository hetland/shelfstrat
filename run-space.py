#!env python

import os

preamble = '''#!/bin/sh

### Set the job name
#PBS -N constr 

### Declare myprogram non-rerunable
#PBS -r n

### Set the queue to "default"
#PBS -q default

### Run on 
#PBS -l nodes=1:ppn=1

### You can override the default 1 hour real-world time limit.
#PBS -l walltime=72:00:00

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR 
'''

for Wc in [200, 300, 400]:
    for Tide in [0.00, 0.25, 0.5, 1.0, 2.0]:
        Qbar = int(500.0 * Wc / 300.0) 
        f = open('constr-Wc%d-Tide%4.2f-Qbar%d.sh' % (Wc, Tide, Qbar), 'w')
        f.writelines(preamble)
        f.write('./run.py --Wc=%d --Tide=%4.2f --Qbar=%d' % (Wc, Tide, Qbar))
        f.close()
        os.system('qsub constr-Wc%d-Tide%4.2f-Qbar%d.sh' % (Wc, Tide, Qbar))
