#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Created on Tuesday Oct 27 11:29am 2020

Converts the .bestprof file into a par file

"""
from __future__ import division, print_function
import numpy as np
import pathlib
import Lv0_dirs

Lv0_dirs.global_par() #obtaining the global parameters

def bestprof2par(bestprof,parfile):
    """
    Converts the .bestprof text file into a .par file

    bestprof - path to the bestprof file
    parfile - path to the parfile
    """
    contents = open(bestprof,'r').read().split('\n')

    for line in contents:
        if line.startswith('# Epoch_bary'):
            pepoch = float(line.split()[4])
        if line.startswith('# P_bary'):
            p0 = float(line.split()[4])/1000.0
        if line.startswith("# P'_bary"):
            p1 = float(line.split()[4])
        if line.startswith("# P''_bary"):
            p2 = float(line.split()[4])

    f0 = 1/p0
    f1 = -p1/p0**2
    f2 = (2*p1**2-p0*p2)/p0**3

    par_write = open(parfile,'w')
    par_write.write("PEPOCH " + str(pepoch) + '\n')
    par_write.write("F0 " + repr(f0) + " 1" + '\n')
    par_write.write("F1 " + repr(f1) + " 1" + '\n')
    par_write.write("F2 " + repr(f2))

    par_write.close()

    return

if __name__ == "__main__":

    bestprof = Lv0_dirs.NICER_DATADIR + 'igrj17494-3030/accelsearch_GTIs/2010270036_bary_GTI000002_E0200-1200_ACCEL_Cand_1.pfd.bestprof'
    parfile = Lv0_dirs.NICER_DATADIR + 'igrj17494-3030/igrj17494.par'
    bestprof2par(bestprof,parfile)
