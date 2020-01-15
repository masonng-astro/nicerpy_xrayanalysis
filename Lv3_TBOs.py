#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 5:44pm 2019

Script that analyzes Leahy-normalized power spectra of overlapping, sliding
windows of a time series to find burst oscillations

1/15/20: Moved methods to Lv2_TBOs_method, and revamped Lv3_TBOs to be a more
executive script than a methods script. I'll be moving search_window under the for
loop soon, when I interface it with a script that will grab search intervals...

"""
from __future__ import division, print_function
import numpy as np
import subprocess
from astropy.io import fits
from scipy import stats
from tqdm import tqdm

import Lv0_dirs,Lv2_TBOs_method
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
import mplcursors

Lv0_dirs.global_par()

eventfiles = ['/Volumes/Samsung_T5/NICERsoft_outputs/2584010501_pipe/ni2584010501_nicersoft_bary.evt']

do_search = True # Searching for burst candidates!
do_plots = True # Creating dynamic power spectra and contours

search_window = [5560,5660] #start/end times of burst oscillation search
T = np.array([10]) #window sizes
dt = T/10  #i.e., do overlapping windows of T seconds in steps of dt seconds
tbin_size = 0.001 #size of time bins in seconds
df = 10 #tolerance of 10 Hz, say
f_central = 401 #central frequency
mode = "show"

for i in range(len(eventfiles)):
    if do_search == True:
        Lv2_TBOs_method.burst_cands(eventfile)

    if do_plots == True:
        Lv2_TBOs_method.dynamic_ps(eventfile,search_window,T,dt,tbin_size,df,f_central,mode)
