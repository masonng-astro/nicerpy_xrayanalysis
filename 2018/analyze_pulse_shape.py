#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 19:08:43 2018

Looking at pulse shapes!

GX349+2 - Obs ID 1034090111

Columns of data: TIME, RAWX, RAWY, PHA, PHA_FAST, DET_ID, DEADTIME, EVENT_FLAGS,
TICK, MPU_A_TEMP, MPU_UNDER_COUNT, PI_FAST, PI, PI_RATIO

1013010125 - Crab pulsar #frequency of ~30Hz
sumcounts_1013010125_1s.txt - 1s bins
sumcounts_1013010125_0p01s.txt - 0.01s bins
sumcounts_1013010125_0p001s.txt - 0.001s bins

1034070104 - Cen X-3 #frequency of ~0.207Hz 
sumcounts_1034070104_0p001s.txt - 0.001s bins
sumcounts_1034070104_0p01s.txt - 0.01s bins

1034070105 - Cen X-3 as well
sumcounts_1034070105_0p1s.txt - 0.1s bins
sumcounts_1034070105_1s.txt - 1s bins

1034090111 - GX 349+2
"""

from __future__ import division
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from scipy import stats
from scipy import signal

timestart = time.time()

work_dir = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/'
obsid = '0034070101' 
t_interval = '0p01s'
place = 100 #for t_interval = 0p1s
bary = True
if bary == True:
    bary = '_bary'
else:
    bary = ''

times, counts_data = np.loadtxt('sumcounts_'+obsid+'_'+t_interval+bary+'.txt', usecols=(0,1), unpack=True) 
#times_unb, counts_data_unb = np.loadtxt('sumcounts_'+obsid+'_'+t_interval+'.txt', usecols=(0,1), unpack=True)

#plt.figure(figsize=(10,8))
#plt.plot(times,counts_data,'rx')
#plt.xlabel('Time (s)')
#plt.ylabel('Counts')
#plt.xlim([1750000,2000000])
#
## Analyze pulse shape:
#Frequency - 29.6297Hz ; T = 0.03374991984s BEFORE barycentering
#Frequency - 29.63271Hz # AFTER barycentering
#
#plt.plot(times[0:1750000],counts_data[0:1750000],'rx')

times_1 = times
counts_1 = counts_data
#
#T = times_1[-1]-times_1[0]
#dt = T/(len(times_1))
##
#plt.figure(figsize=(10,8))
#mean_corrected = counts_data-np.mean(counts_data)
#power_spec = 2.0/len(times)*np.abs(np.fft.fft(mean_corrected))**2
#freqs = np.fft.fftfreq(mean_corrected.size, dt)
#N = len(freqs)
###freqs = np.linspace(0,1/(2.0*dt),int(len(test_times)/2)) #correct method! Same as that of fftfreq
#plt.semilogy(freqs[1:int(N/2)],power_spec[1:int(N/2)]/np.mean(power_spec[1:int(N/2)]),'r-')
#plt.xlabel('Hz')
#plt.ylabel('Some normalized power spectrum')
#plt.xlim([29.632,29.634])
#plt.axvline(x=29.63271,color='k',alpha=0.5,lw=0.5)

###############################################################################
############################ FOLDING LIGHT CURVE ##############################
###############################################################################

#https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/folding.html
from PyAstronomy.pyasl import foldAt
#f_pulsation = 29.63271
f_pulsation = 0.208
period = 1.0/f_pulsation
phases = foldAt(times_1,period,T0=0.45*period)
#
index_sort = np.argsort(phases)
phases = list(phases[index_sort]) + list(phases[index_sort]+1)
counts = list(counts_1[index_sort])*2

phase_bins = np.linspace(0,2,1001)
summed_phases, bin_edges, binnumber = stats.binned_statistic(phases,counts,statistic='sum',bins=phase_bins)
#
plt.figure(figsize=(10,8))
#
plt.plot(phase_bins[:-1],summed_phases,'rx')


timeend = time.time()


print (timeend-timestart)