#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 18:16:53 2018

Let us do count rate vs angular offset between pointing and source

1013010125 - Crab pulsar #frequency of ~30Hz
sumcounts_1013010125_1s.txt - 1s bins
sumcounts_1013010125_0p01s.txt - 0.01s bins
sumcounts_1013010125_0p001s.txt - 0.001s bins

1034070101 - Cen X-3 #frequency of ~0.207Hz

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
place = 10 

print(datetime.datetime.now())

work_dir = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/'
#obsid = '1034070101' 
obsid = '1034070103'

event = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl.evt'
event = fits.open(event) 

###############################################################################
###############################################################################

## file containing the attitude information. If the observation includes the slew,
## both slew and pointing attitude are present
#aux_att = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/auxil/ni' + obsid + '.att'
#aux_att = fits.open(aux_att)
### Has 2 cards, "ATTITUDE" - data from ST and "INST_ATTITUDE" - data from INST; I assume instrument attitude?
#
################################################################################
################################################################################
#
### file containing the list of files output of the processing that are included
### within the sequences in the archive
#aux_cat = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/auxil/ni' + obsid + '.cat'
#aux_cat = fits.open(aux_cat)
#
################################################################################
################################################################################
#
### Make filter file containing a subset of HK parameters from the individual
### instrument, the subsystems, as well as from the spacecraft. This information
### is used to screen the data
aux_mkf = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/auxil/ni' + obsid + '.mkf'
aux_mkf = fits.open(aux_mkf)
#
################################################################################
################################################################################
#
### File containing the orbit information
#aux_orb = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/auxil/ni' + obsid + '.orb'
#aux_orb = fits.open(aux_orb)
#
################################################################################
################################################################################
#
#hk_mpu = '0mpu4' #from 0mpu0 to 0mpu6
#housekeeping = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/xti/hk/ni' + obsid + '_' + hk_mpu + '.hk'
#housekeeping = fits.open(housekeeping)
#
################################################################################
################################################################################
#
event = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl.evt'
event = fits.open(event) 
##see event.info() for each card
##event[1].header for events; event[2].header for GTIs
#
################################################################################
################################################################################
#
### Get data - counts, angular distance, and time from aux_mkf
tot_xray_count = aux_mkf[1].data['TOT_XRAY_COUNT'] #sum of undercounts, xray counts and overcounts
ang_dist = aux_mkf[1].data['ANG_DIST']
times_aux = aux_mkf[1].data['TIME']
#
##first filter out the nans in the data
pi_fast = event[1].data['PI_FAST']
print 'Done with PI_FAST'
truncated_pi_fast = pi_fast[pi_fast>=0]
print 'Done with truncated PI_FAST'
times = event[1].data['TIME'][pi_fast>=0]
print 'Done with TIME'
pi = event[1].data['PI'][pi_fast>=0]
print 'Done with PI'
pi_ratio = event[1].data['PI_RATIO'][pi_fast>=0]
print 'Done with PI_RATIO'
flags = event[1].data['EVENT_FLAGS'][pi_fast>=0]
print 'Done with FLAGS'
#
#
################################################################################
### Plot TOT_XRAY_COUNT vs TIME in blue; then TOT_XRAY_COUNT vs TIME in red for GTIs
#
#plt.figure(figsize=(10,8))
#plt.plot(times_aux,tot_xray_count,'bx')
#plt.plot(times_aux[(times_aux>=times[0])&(times_aux<=times[-1])], tot_xray_count[(times_aux>=times[0])&(times_aux<=times[-1])], 'rx')
#plt.xlabel('MJD Time (s)',fontsize=15)
#plt.ylabel('XRAY counts',fontsize=15)
#plt.legend(('All counts', 'Counts corresponding to GTIs'),loc='best')
#plt.xlim([1.12403642e+08,times_aux[-1]])
#
################################################################################
### PLOT ANG_DIST vs TIME in blue; then ANG_DIST vs TIME in red for GTIs
#
#plt.figure(figsize=(10,8))
#plt.plot(times_aux,ang_dist,'bx')
#plt.plot(times_aux[(times_aux>=times[0])&(times_aux<=times[-1])], ang_dist[(times_aux>=times[0])&(times_aux<=times[-1])], 'rx')
#plt.xlabel('MJD Time (s)',fontsize=15)
#plt.ylabel('Angular Distance (Deg)',fontsize=15)
#plt.legend(('All counts', 'Counts corresponding to GTIs'),loc='best')
#plt.xlim([1.12403642e+08,times_aux[-1]])
#plt.ylim([0,0.15])
#
################################################################################
### PLOT TOT_XRAY_COUNT vs ANG_DIST in blue; then TOT_XRAY_COUNT vs ANG_DIST in red for GTIs
#
#plt.figure(figsize=(10,8))
#plt.plot(ang_dist,tot_xray_count,'bx')
#plt.plot(ang_dist[(times_aux>=times[0])&(times_aux<=times[-1])], tot_xray_count[(times_aux>=times[0])&(times_aux<=times[-1])], 'rx')
#plt.xlabel('Angular Distance (Deg)',fontsize=15)
#plt.ylabel('XRAY counts',fontsize=15)
#plt.legend(('All counts', 'Counts corresponding to GTIs'),loc='best')
#plt.xlim([0,0.2])
#
################################################################################
#bg_counts = tot_xray_count[(tot_xray_count>200)&(tot_xray_count<750)&(ang_dist>0.05)&(ang_dist<0.11)] ## THIS IS PER SECOND!!!
#bg_time = times_aux[(tot_xray_count>200)&(tot_xray_count<750)&(ang_dist>0.05)&(ang_dist<0.11)]
#bg_ang_dist = ang_dist[(tot_xray_count>200)&(tot_xray_count<750)&(ang_dist>0.05)&(ang_dist<0.11)]
#
################################################################################
#### TRY A QUICK FOURIER TRANSFORM
#
#truncated_t = times-times[0] #to shift the times
#t_bins = np.linspace(0,int(truncated_t[-1]),int(truncated_t[-1])*place+1) #gives bins of 0.1s! 
#summed_data, bin_edges, binnumber = stats.binned_statistic(truncated_t,truncated_pi_fast,statistic='sum',bins=t_bins)
#subtracted_counts = summed_data - np.mean(bg_counts)/place 
#
#T = times[-1]-times[0]
#dt = T/(len(times))
#plt.figure(figsize=(10,8))
##mean_corrected = subtracted_counts
#mean_corrected = subtracted_counts-np.mean(subtracted_counts)
#power_spec = 2.0/len(times)*np.abs(np.fft.fft(mean_corrected))**2
#freqs = np.fft.fftfreq(mean_corrected.size, dt)
#N = len(freqs)
#####freqs = np.linspace(0,1/(2.0*dt),int(len(test_times)/2)) #correct method! Same as that of fftfreq
#plt.semilogy(freqs[1:int(N/2)],power_spec[1:int(N/2)]/np.mean(power_spec[1:int(N/2)]),'r-')
#plt.xlabel('Hz')
#plt.ylabel('Some normalized power spectrum')
#plt.xlim([0,5])
#plt.axvline(x=0.2082,color='k',alpha=0.5,lw=0.5)

timeend = time.time()


print (timeend-timestart)


