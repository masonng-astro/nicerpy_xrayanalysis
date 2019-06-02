#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 9:59pm 2019

Calculating the duty cycle in the data.

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
from scipy import stats
import glob

Lv0_dirs.global_par()

def duty_cycle(obsid,tbin,segment_length,duty_cycle_bin,threshold):
    """
    To determine the percentage of "data used"/total amount of data. Have two
    types of values:
    1) % of bins (of size duty_cycle_bin) with data over the ENTIRE observation
    2) % of bins with data over POTENTIALLY USABLE segments (recall that if the event
       file is completely empty, then nicerfits2presto would not process it!)
    3) % of bins with data ACTUALLY USED (so with a threshold defined)

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    duty_cycle_bin - binning used to calculate the duty cycle
    threshold - if amount of data in the segment is more than threshold IN PERCENTAGE, use the data
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsdir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    dat_files = sorted(glob.glob(obsdir+'*GTI*'+str(segment_length)+'s.dat')) #grab the .dat files with a specified segment length

    event = obsdir + 'ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data
    obs_duration = gtis[-1][1] - gtis[0][0] #to get the duration of the ENTIRE observation

    available_data = 0 #non-empty bins (binned by duty_cycle_bin seconds long)
    useful_data = 0
    segment_duration = 0
    for i in range(len(dat_files)):
        binned_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
        dat_times = np.arange(0,tbin*len(binned_data),tbin)

        duty_cycle_times = np.arange(0,tbin*len(binned_data)+tbin,duty_cycle_bin)
        summed_data, binedges, binnumber = stats.binned_statistic(dat_times,binned_data,statistic='sum',bins=duty_cycle_times)
        #plt.plot(duty_cycle_times[:-1],summed_data)
        #plt.title(dat_files[i])
        #plt.show()

#        print(stats.mode(summed_data))
        usable_data = len(summed_data[(summed_data!=stats.mode(summed_data)[0][0])&(summed_data>0)])
        available_data += usable_data
        if usable_data/len(summed_data)*100 >= threshold:
            useful_data += usable_data
        segment_duration += len(summed_data)

    print(available_data,useful_data,segment_duration,obs_duration)
    print("Percentage of bins with data over the ENTIRE observation: " + str(available_data/obs_duration*100))
    print("Percentage of bins over POTENTIALLY USABLE segments: " + str(available_data/segment_duration*100))
    print("Percentage of bins over ACTUALLY USED segments: " + str(useful_data/segment_duration*100))

    return

#print('0034070101')
#duty_cycle('0034070101',0.0005,100,1,10)
#print('1034090111')
#duty_cycle('1034090111',0.00025,1000,1,10)

## so for each ObsID,
## for each .dat file (so nicerfits2presto would already have taken out empty data)
## bin the data (total time is segment_length!) by time, and then
## do the threshold thing
## and calculate the duty cycle from there
## maybe also calculate the total amount of data/total observation time!
