#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 11:51:21 2018

@author: masonng

----------------------- !!! SPECTRAL ANALYSIS !!! -----------------------------

GX349+2 - Obs ID 1034090111

Columns of data: TIME, RAWX, RAWY, PHA, PHA_FAST, DET_ID, DEADTIME, EVENT_FLAGS,
TICK, MPU_A_TEMP, MPU_UNDER_COUNT, PI_FAST, PI, PI_RATIO

Recall that PI=energy/10 eV, pi=110 is 1.10keV

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

1103010117 - GRS 1915+105
Just to check my code for creating the hardness ratios:
https://arxiv.org/pdf/1806.02342.pdf?fbclid=IwAR0H78FWyi_sIyIWk0Z0oO09nyzs329xS9sJHiVBK1OHfNX2MFcA-POxkOg
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
obsid = '1034070104'  #Cen X-3
#obsid = '1103010117' ## GRS 1915+105, from that paper by Neilsen et al. 2018
t_interval = '0p01s'
extra = ''
place = 100 #for t_interval = 0p1s
cutoff = 4.0 #keV

event = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl.evt'
event = fits.open(event) #see event.info() for each card
#event[1].header for events; event[2].header for GTIs

#recall: energy stuff is given in PI, so don't do PI_FAST truncation
pi_fast = event[1].data['PI_FAST']
print 'Done with PI_FAST'
times = event[1].data['TIME']
print 'Done with TIME'
pi_data = event[1].data['PI']
print 'Done with PI'
pi_ratio = event[1].data['PI_RATIO']
print 'Done with PI_RATIO'
flags = event[1].data['EVENT_FLAGS']
print 'Done with FLAGS'
#plt.figure(figsize=(10,8))
#plt.plot(times[0:10000000], truncated_pi_fast[0:10000000],'x')

######################### GETTING THE LIGHT CURVE #############################

pi_fast_noneg = pi_fast
np.place(pi_fast_noneg,pi_fast_noneg<0,0) #replace negative values by 0

sort_ind = np.argsort(pi_data)
sorted_E = pi_data[sort_ind]*10/1000
sorted_counts = pi_fast_noneg
plt.figure(figsize=(10,8))
plt.plot(sorted_E,sorted_counts,'rx')

def get_cutoff(sorted_E,sorted_counts):
    #given a SORTED array of counts vs energy (in keV), find the corresponding
    #E value such that you get 50% of counts on either side
    
    half_counts = int(0.5*sum(sorted_counts)) 
    tally = 0
    for i in tqdm(range(len(sorted_counts))):
        if tally < half_counts:
            tally += sorted_counts[i]
            energy = sorted_E[i]
        else:
            break
    
    return energy, i

cutoff,index = get_cutoff(sorted_E,sorted_counts)
print cutoff,index
    

shifted_t = times - times[0]
t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1]+2))
#bin into counts per second
summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,pi_fast,statistic='sum',bins=t_bins)

#plt.plot(t_bins[:-1],summed_data,'rx')

######################### GETTING THE ENERGIES ################################
def soft_counts(cutoff,shifted_t,pi_data,pi_fast_noneg_soft):
    cutoff_PI = cutoff*1000/10 #get the cut-off in terms of PI 

    #we do assume that there will be NO bins in which you get ZERO counts
    np.place(pi_fast_noneg_soft,pi_data>=cutoff_PI,0) #get the counts for soft photons

    return pi_fast_noneg_soft

#print soft_counts(cutoff,shifted_t,pi_data,pi_fast_noneg)[0:50]

#event = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl.evt'
#event = fits.open(event)
#pi_fast = event[1].data['PI_FAST']
#print 'Done with PI_FAST'
#print pi_fast[0:50]

#pi_fast_noneg = pi_fast
#np.place(pi_fast_noneg,pi_fast_noneg<0,0) #replace negative values by 0

def hard_counts(cutoff,shifted_t,pi_data,pi_fast_noneg_hard):
    cutoff_PI = cutoff*1000/10 #get the cut-off in terms of PI 

    #we do assume that there will be NO bins in which you get ZERO counts
    np.place(pi_fast_noneg_hard,pi_data<cutoff_PI,0) #get the counts for soft photons

    return pi_fast_noneg_hard

#print hard_counts(cutoff,shifted_t,pi_data,pi_fast_noneg)[0:50]


def color_plot(cutoff,shifted_t,pi_data,pi_fast_noneg):
    # cutoff is the energy cutoff for the hardness ratio/color
    # shifted_t is where the time series starts at t=0
    # pi_data is the corresponding energy tag for each event
    # pi_fast_noneg is the corresponding counts 

    #for example, for the soft photons below some cutoff energy, replace each
    #element in the array corresponding to COUNTS by 0, IF the corresponding
    #energy value is ABOVE the cutoff (i.e., filter out higher E photons)

    pi_fast_soft = soft_counts(cutoff,shifted_t,pi_data,pi_fast_noneg)

    event = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/' + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl.evt'
    event = fits.open(event)
    pi_fast = event[1].data['PI_FAST']
    pi_fast_noneg = pi_fast
    np.place(pi_fast_noneg,pi_fast_noneg<0,0) #replace negative values by 0

    pi_fast_hard = hard_counts(cutoff,shifted_t,pi_data,pi_fast_noneg)

    t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1]+2)) #get 1-second bins
    sum_soft, bin_edges_soft, binnumber_soft = stats.binned_statistic(shifted_t,pi_fast_soft,statistic='sum',bins=t_bins)
    sum_hard, bin_edges_hard, binnumber_hard = stats.binned_statistic(shifted_t,pi_fast_hard,statistic='sum',bins=t_bins)
    np.place(sum_soft,sum_soft==0,1) #so that you get 0/1 instead

#    color = (sum_hard-sum_soft)/(sum_soft+sum_hard)
    color = sum_hard/sum_soft

    plt.figure(figsize=(10,8))
    plt.plot(t_bins[:-1],color,'r-')
    plt.ylim([0,20])
    plt.xlim([11100,11900])
    plt.ylabel('Color')
    plt.xlabel('Time(s)')
    
color_plot(cutoff,shifted_t,pi_data,pi_fast_noneg)

#energies = pi*10/1000 #to convert into keV 
#E_max = max(energies)
#E_min = min(energies)
#
#pi_bins = np.linspace(pi_min,pi_max,(pi_max-pi_min)+1)
#
#altered_pi_fast = pi_fast
#altered_pi_fast[np.where(pi_fast<0)]=0

#summed_data, bin_edges, binnumber = stats.binned_statistic(pi,pi_fast,statistic='sum',bins=pi_bins)
#
#E_bins = pi_bins*10/1000
#
#plt.figure(figsize=(10,8))
#plt.plot(E_bins[:-1],summed_data)
#plt.ylabel('Counts')
#plt.xlabel('keV')


timeend = time.time()


print (timeend-timestart)