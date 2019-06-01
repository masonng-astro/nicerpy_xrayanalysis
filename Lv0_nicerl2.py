#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs May 30 9:31am 2019

Program for doing nicerl2 (mix of nicercal, niprefilter2, nimaketime, nicermergeclean)

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs
from tqdm import tqdm
import os
import subprocess
import glob

Lv0_dirs.global_par()

def nicerl2(obsid,nicerl2_flags):
    """
    Running nicerl2 to do initial filtering of the ufa file!

    obsid - Observation ID of the object of interest (10-digit str)
    flags - a LIST of input flags for nicerl2
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(nicerl2_flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    logfile = obsid + '_nicerl2.log'
    nicerdata_dir = Lv0_dirs.NICER_DATADIR + obsid

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['nicerl2',nicerdata_dir]+nicerl2_flags))
        logtextfile.close()

    subprocess.check_call(['mv',obsid+'_nicerl2.log',nicerdata_dir])

    return
