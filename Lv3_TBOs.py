#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 5:44pm 2019

Script that analyzes Leahy-normalized power spectra of overlapping, sliding
windows of a time series to find burst oscillations

"""
from __future__ import division, print_function
import numpy as np
import subprocess
from astropy.io import fits
from scipy import stats
from tqdm import tqdm

import Lv0_dirs,Lv2_ps_method,Lv3_detection_level
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
import mplcursors

Lv0_dirs.global_par()

### write with the old functions first (e.g., with ObsIDs and old barycorr) then
### alter next week or so... now just get Lv3_TBOs.py working!

### preprocessing the data and inspecting it

### specify parameters here
event_file = Lv0_dirs.NICERSOFT_DATADIR + '2584010501_pipe/ni2584010501_nicersoft_bary.evt'
MJDREFI = fits.open(event_file)[0].header['MJDREFI']
MJDREFF = fits.open(event_file)[0].header['MJDREFF']
source_name = 'J1808.4-3658'
obsid = '2584010501'
search_start = 5560
search_end = 5660

T = np.array([10]) #window sizes
dt = T/10 #i.e., do overlapping windows of T seconds in steps of dt seconds
tbin_size = 0.001 #size of time bins in seconds
significance = 3 #in no. of standard deviations

df = 10 #tolerance of say, 10 Hz
f_central = 401 #central frequency

### get the time series and zero-ise it
#define an array of start times? So do in steps of dt from search_start to search_end
times = fits.open(event_file)[1].data['TIME']
T_zeroized = times-times[0]
counts = np.ones(len(T_zeroized))

T_bins = np.linspace(0,np.ceil(T_zeroized[-1]),np.ceil(T_zeroized[-1])*1/tbin_size+1)
binned_counts, bin_edges, binnumber = stats.binned_statistic(T_zeroized,counts,statistic='sum',bins=T_bins) #binning the photons

for i in tqdm(range(len(T))): #for every window size:
    output_file = open(Lv0_dirs.NICERSOFT_DATADIR + '2584010501_pipe/TBO_search_' + str(T[i]) + 's.txt','w')
    output_file.write('Source name: ' + source_name + ' ; ObsID: ' + obsid + '\n')
    output_file.write('Window size: T = ' + str(T[i]) + 's, stepping size = ' + str(dt[i]) + ' ; dt = ' + str(tbin_size) + '\n')

    T_start = np.arange(search_start,search_end,dt[i]) #start time of each sliding window
    T_end = T_start + T[i] #end time of each sliding window
    N = T[i]/tbin_size #number of trials for each window
    sig3 = Lv3_detection_level.power_for_sigma(3,N,1,1)
    sig4 = Lv3_detection_level.power_for_sigma(4,N,1,1)
    sig5 = Lv3_detection_level.power_for_sigma(5,N,1,1)
    power_required = Lv3_detection_level.power_for_sigma(significance,N,1,1)

    output_file.write('Power needed for: 3 sigma - ' + str(sig3) + ' ; 4 sigma - ' + str(sig4) + ' ; 5 sigma - ' + str(sig5) + '\n')
    output_file.write('Starting/Ending MJD of TBO search scheme: ' + str(MJDREFI+MJDREFF+(times[0]+search_start)/86400) + '/' + str(MJDREFI+MJDREFF+(times[0]+search_end)/86400) + '\n')

    fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True) #dynamic power spectrum #define a 2x1 subplot or something
    fig.subplots_adjust(hspace=0)

    f_max = [] #corresponding frequencies to the maximum power
    ps_max = [] #to store the maximum power from each power spectrum of each sliding time series
    for j in tqdm(range(len(T_start))): #for every sliding window
        T_search = T_bins[:-1][(T_bins[:-1]>=T_start[j])&(T_bins[:-1]<=T_end[j])] #time series to search for burst oscillations
        binned_search = binned_counts[(T_bins[:-1]>=T_start[j])&(T_bins[:-1]<=T_end[j])]

        f,ps = Lv2_ps_method.manual(T_search,binned_search,[False,400,500],[False,400],False,[False,5]) #calculating Leahy-normalized power spectra
        f_window = f[(f>=f_central-df)&(f<=f_central+df)] #get frequency values 'df' Hz about f_central
        ps_window = ps[(f>=f_central-df)&(f<=f_central+df)] #get powers 'df' Hz about f_central

        scatt = ax1.scatter(x=np.ones(len(f_window))*T_start[j],y=f_window,s=12,c=ps_window,marker='o',cmap=cm.gist_heat,vmin=1,vmax=50)

        f_max.append(f_window[ps_window==np.max(ps_window)][0])
        ps_max.append(ps_window[ps_window==np.max(ps_window)][0])

        output_file.write('Start time for this window: zeroized - ' + str(T_start[j]) + ' ; MJD - ' + str(MJDREFI+MJDREFF+(times[0]+T_start[j])/86400) + '\n')
        for k in range(len(f_window)):
            output_file.write(str(f_window[k]) + ' ' + str(ps_window[k]) + ' ' + str(Lv3_detection_level.signal_significance(1,1,ps_window[k])) + '\n')

    output_file.close()

    mplcursors.cursor(hover=True)
    ax1.set_title('Window size: ' + str(T[i]) + 's, dt='+str(dt[i])+'s \n' + 'Central freq. = '+str(f_central) + 'Hz, df = ' + str(df) + 'Hz \n Power required for ' + str(significance) + ' sigma: ' + str(power_required),fontsize=12)
    ax1.set_ylabel('Frequency (Hz)',fontsize=12)
    ax1.set_ylim([f_central-df,f_central+df])


    ax2.set_ylabel('Frequency (Hz)',fontsize=12)

    scat = ax2.scatter(x=T_start,y=f_max,s=12,c=ps_max,marker='o',cmap=cm.gist_heat,vmin=np.min(ps_max),vmax=np.max(ps_max),edgecolors='k')
    mplcursors.cursor(hover=True)

    #fig.colorbar(scat,ax=ax1)
    #fig.colorbar(scat,ax=ax2)
    ax2.set_ylim([f_central-df,f_central+df])

    ps_contour = ax3.tricontour(T_start,f_max,ps_max,levels=30,linewidths=0.5,colors='k')
    ax3.clabel(ps_contour,fontsize=8)
    ax3.set_xlabel('Time (s)')
    ax3.set_ylim([f_central-df,f_central+df])
    mplcursors.cursor(hover=True)

plt.show()

### calculate the significance - print a couple of things:
### source name
### obsid

### produce a text file which has something like:
### text file name: $ObsID.txt
### window size - xx seconds, dt - xx seconds,
### 3 sigma - [Leahy power] ; 4 sigma - [Leahy power] ; 5 sigma - [Leahy power]
### starting MJD, end MJD
### ##### [start time 1 - in zeroized units]
### frequency - Leahy power - significance
### frequency - Leahy power - significance
### ...
### [start time 2 - in zeroized units]
