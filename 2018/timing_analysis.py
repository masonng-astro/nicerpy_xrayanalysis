#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 22:13:37 2018

Simple test routines to find pulsation signals? First do it with 
GX349+2 - Obs ID 1034090111

Columns of data: TIME, RAWX, RAWY, PHA, PHA_FAST, DET_ID, DEADTIME, EVENT_FLAGS,
TICK, MPU_A_TEMP, MPU_UNDER_COUNT, PI_FAST, PI, PI_RATIO

1013010125 - Crab pulsar #frequency of ~30Hz
sumcounts_1013010125_1s.txt - 1s bins
sumcounts_1013010125_0p01s.txt - 0.01s bins
sumcounts_10130101250p001s.txt - 0.001s bins

0034070103 - Cen X-3 [Commissioning phase]

1034070104 - Cen X-3 #frequency of ~0.207Hz 
sumcounts_1034070104_0p001s.txt - 0.001s bins
sumcounts_1034070104_0p01s.txt - 0.01s bins
sumcounts_1034070104_0p1s.txt - 0.1s bins

1034070105 - Cen X-3 as well
sumcounts_1034070105_0p1s.txt - 0.1s bins
sumcounts_1034070105_1s.txt - 1s bins

1034090111 - GX 349+2
"""
from __future__ import division
import datetime
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from scipy import stats
from scipy import signal

timestart = time.time()

print(datetime.datetime.now())

work_dir = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/'
obsid = '0034070104' 
t_interval = '0p01s'
extra = '_angdist0p035'
place = 100 #for t_interval = 0p1s

event = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl.evt'
event = fits.open(event) #see event.info() for each card
#event[1].header for events; event[2].header for GTIs

#first filter out the nans in the data
#pi_fast = event[1].data['PI_FAST']
#print 'Done with PI_FAST'
#truncated_pi_fast = pi_fast[pi_fast>=0]
#print 'Done with truncated PI_FAST'
times_d = event[1].data['TIME'][pi_fast>=0]
print 'Done with TIME'
#pi = event[1].data['PI'][pi_fast>=0]
#print 'Done with PI'
#pi_ratio = event[1].data['PI_RATIO'][pi_fast>=0]
#print 'Done with PI_RATIO'
#flags = event[1].data['EVENT_FLAGS'][pi_fast>=0]
#print 'Done with FLAGS'
##plt.figure(figsize=(10,8))
##plt.plot(times[0:10000000], truncated_pi_fast[0:10000000],'x')
##
#shifted_t = times_d - times_d[0]
#t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])*place+1) #gives bins of 0.01s! 
#####
##plt.figure(figsize=(10,8))
#summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,truncated_pi_fast,statistic='sum',bins=t_bins)
##plt.plot(t_bins[:-1],summed_data,'x')
###plt.xlim([0,1000])
##timestep = 1/(t_bins[-1]) #because the start time is 0! 
####
#countsfile = open("sumcounts_"+obsid+"_"+t_interval+extra+'.txt', "w+")
#for i in tqdm(range(len(summed_data))):
#    countsfile.write(str(i/place) + ' ' + str(summed_data[i]) + '\n')
####
#countsfile.close()
##
#plt.figure(figsize=(10,8))
times, counts = np.loadtxt('sumcounts_'+obsid+'_'+t_interval+extra+'.txt', usecols=(0,1), unpack=True) 
#plt.plot(times,counts,'rx')
#plt.xlabel('Time (s)')
#plt.ylabel('Counts')
#plt.xlim([0,5000])
#plt.ylim([1750000,2000000])

#counts_1 = np.zeros(len(counts))
##set counts = 1
#for i in tqdm(range(len(counts))):
#    if counts[i] > 0:
#        counts_1[i] = 1
#
#plt.figure(figsize=(10,8))
#plt.plot(times,counts_1,'rx')


################################ OVERSAMPLING #################################
#pad_zeros = np.zeros(len(times)*4) #i.e., add 4T worth of 0s onto the data
#counts = np.array(list(counts) + list(pad_zeros))
#times = np.arange(times[0],times[-1]*5+0.005,0.01)

times_1 = times
counts_1 = counts
####
####
T = times_1[-1]-times_1[0]
dt = T/(len(times_1))
####
###################################### PERIODOGRAM #################################
####
plt.figure(figsize=(10,8))
freq_sample = 1.0/dt
f,pxx = signal.periodogram(counts_1,fs=freq_sample)
plt.semilogy(f[1:],pxx[1:]/np.mean(pxx[1:]))
plt.xlabel('Hz')
plt.ylabel('Some normalized power spectrum')
plt.xlim([0,1])
plt.axvline(x=0.2082,color='k',alpha=0.5,lw=0.5)

### find the frequencies with some sort of peaks
##cutoff = f[pxx/np.mean(pxx)>10]
##harmonic_f = cutoff[cutoff>5]
################################# MANUAL FFT ###################################
####
plt.figure(figsize=(10,8))
mean_corrected = counts_1-np.mean(counts_1)
power_spec = 2.0/len(times_d)**2*np.abs(np.fft.fft(mean_corrected))**2
freqs = np.fft.fftfreq(mean_corrected.size, dt)
N = len(freqs)
####freqs = np.linspace(0,1/(2.0*dt),int(len(test_times)/2)) #correct method! Same as that of fftfreq
plt.semilogy(freqs[1:int(N/2)],power_spec[1:int(N/2)]/np.mean(power_spec[1:int(N/2)]),'r-')
plt.xlabel('Hz')
plt.ylabel('Some normalized power spectrum')
plt.xlim([0,1])
plt.axvline(x=0.2082,color='k',alpha=0.5,lw=0.5)

#
timeend = time.time()


print (timeend-timestart)

################################ TRIAL #######################################
#
#import numpy as np
#import matplotlib.pyplot as plt
#import scipy.fftpack

# Number of samplepoints
#N = 800
## sample spacing
#T = 1.0 / 800.0
#x = np.linspace(0.0, N*T, N)
#y = 10 + np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
#yf = np.fft.fft(y)
##yf = scipy.fftpack.fft(y)
#xf = np.linspace(0.0, 1.0/(2.0*T), int(N/2))
#
#plt.figure(figsize=(10,8))
##plt.subplot(2, 1, 1)
##plt.plot(xf, 2.0/N * np.abs(yf[0:int(N/2)]))
#plt.subplot(2, 1, 2)
#plt.plot(xf[1:], 2.0/N * np.abs(yf[0:int(N/2)])[1:])

## sample binning algorithm
#
#trial_t = np.array([0.11,0.13,0.15,0.19,0.21,0.22,0.23,0.25,0.28,0.31,0.33,0.335,0.35,0.39,0.395,0.412,0.44,0.456,0.462,0.49])
#trial_values = np.array([1,2,3,4, 5,6,7,8,9, 10,11,12,13,14,15, 16,17,18,19,20])
#bin_edge = np.linspace(0.1,0.5,5)
#
#binned_data = stats.binned_statistic(trial_t,trial_values,statistic='sum',bins=bin_edge)