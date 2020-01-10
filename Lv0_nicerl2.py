#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thurs May 30 9:31am 2019

Program for doing nicerl2 (mix of nicercal, niprefilter2, nimaketime, nicermergeclean)

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs
from tqdm import tqdm
import os
import subprocess
import glob

Lv0_dirs.global_par()

def nicerl2(obsdir,nicerl2_flags):
    """
    Running nicerl2 to do initial filtering of the ufa file!

    obsdir - NICER data directory containing all the data files (e.g., path_to_NICER_dir/1034070101)
    nicerl2_flags - a LIST of input flags for nicerl2
    """
    if type(nicerl2_flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    logfile = obsdir + '/nicerl2.log'
    with open(logfile,'w') as logtextfile:
        output = subprocess.run(['nicerl2',obsdir]+nicerl2_flags,capture_output=True,text=True)
        logtextfile.write(output.stdout)
        logtextfile.write('*------------------------------* \n')
        logtextfile.write(output.stderr)
        logtextfile.close()

    return

if __name__ == "__main__":
    obsid = '1034070101'
    obsdir = Lv0_dirs.NICER_DATADIR + obsid
    nicerl2(obsdir,['clobber=YES'])
