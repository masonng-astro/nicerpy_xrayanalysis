#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues April 14 10:57am 2020

Very quick analysis of time-resolved spectroscopy results. I need to generalize this
and the previous scripts (used for NGC300 ULX-1) during the summer...!

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv1_data_gtis,Lv2_presto_subroutines,Lv2_mkdir
import os
from scipy import stats
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess
import pathlib
import glob

Lv0_dirs.global_par()

tbin_size = 16
eventfile = Lv0_dirs.NICER_DATADIR + 'xtej1812/2003250926_bary.evt'
times = fits.open(eventfile)[1].data['TIME']
trunc = times - times[0]
t_bins = np.linspace(0,np.ceil(trunc[-1]),np.ceil(trunc[-1])*1/tbin_size+1)
summed_data, bin_edges, binnumber = stats.binned_statistic(trunc,np.ones(len(times)),statistic='sum',bins=t_bins)

spectra_GTIs = np.array([1,2,3,4,7,8,9,11,12,13])
gti_start,gti_stop = np.loadtxt('/Volumes/Samsung_T5/NICER-data/xtej1812/bunchedgtis.txt',usecols=(0,1),unpack=True)

gti_used_start = np.array([gti_start[i-1] for i in spectra_GTIs]) - times[0]
gti_used_stop = np.array([gti_stop[i-1] for i in spectra_GTIs]) - times[0]
#for i in range(len(gti_used_start)):
#    print(gti_used_stop[i]-gti_used_start[i])
gti_used_centroid = gti_used_start + (gti_used_stop-gti_used_start)/2

contents = open('/Volumes/Samsung_T5/NICER-data/xtej1812/accelsearch_GTIs/tbabs-bbodyrad.txt','r').read().split('\n')
phoindex = []
phoindex_err = []
norm = []
norm_err = []

for i in range(len(contents)): #for each line
    if "kT" in contents[i]:
        line = [x for x in contents[i].split(' ') if x]
        phoindex.append(float(line[5]))
        phoindex_err.append(float(line[7]))
    elif "norm" in contents[i]:
        line = [x for x in contents[i].split(' ') if x]
        norm.append(float(line[4]))
        norm_err.append(float(line[6]))

fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

time_bins = t_bins[:-1]
count_bins = summed_data/tbin_size
ax1.plot(time_bins[count_bins>0],count_bins[count_bins>0],'rx') #plot the light curve
for i in range(len(gti_used_start)):
    ax1.axvline(x=gti_used_start[i],alpha=0.5,lw=0.5)
    ax1.axvline(x=gti_used_stop[i],alpha=0.5,lw=0.5)

ax1.set_ylabel('Counts/s',fontsize=12)

ax2.errorbar(x=gti_used_centroid,y=phoindex,yerr=phoindex_err,fmt='x-')
ax2.set_ylabel('kT (keV)',fontsize=12)

ax3.errorbar(x=gti_used_centroid,y=norm,yerr=norm_err,fmt='x-')
ax3.set_ylabel('norm',fontsize=12)
ax3.set_xlabel('Time from first event (s)',fontsize=12)

plt.show()
