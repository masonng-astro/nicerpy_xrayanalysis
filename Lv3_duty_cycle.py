#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 9:59pm 2019

Calculating the duty cycle in the data.

Jun 6 - Added duty_cycle_tE and duty_cycle_tE_dist. Need to add duty_cycle_E
and duty_cycle_E_dist in the future!

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
from scipy import stats
import subprocess
import glob

Lv0_dirs.global_par()

def duty_cycle(obsid,tbin,segment_length,duty_cycle_bin,threshold):
    """
    To determine the percentage of "data used"/total amount of data. Have two
    types of values:
    1) % of bins (of size duty_cycle_bin) with data over the ENTIRE observation
    2) % of bins with data over the GTIs

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    duty_cycle_bin - binning used to calculate the duty cycle
    threshold - if amount of data in the segment is more than threshold IN PERCENTAGE, use the data
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsdir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(obsdir+'*GTI*'+str(segment_length)+'s*.dat')) #grab the .dat files with a specified segment length

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data
    obs_duration = gtis[-1][1] - gtis[0][0] #to get the duration of the ENTIRE observation

    total_gti = sum([ (gtis[i][1]-gtis[i][0]) for i in range(len(gtis))])

    useful_data = 0
    segment_duration = 0
    for i in range(len(dat_files)):
        tbin = np.float(tbin)
        binned_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
        dat_times = np.arange(0,tbin*len(binned_data),tbin)

        duty_cycle_times = np.arange(0,tbin*len(binned_data)+tbin,duty_cycle_bin)
        summed_data, binedges, binnumber = stats.binned_statistic(dat_times,binned_data,statistic='sum',bins=duty_cycle_times)
        #plt.plot(duty_cycle_times[:-1],summed_data)
        #plt.title(dat_files[i])
        #plt.show()

#        print(stats.mode(summed_data))
        usable_data = len(summed_data[(summed_data!=stats.mode(summed_data)[0][0])&(summed_data>0)])
        if usable_data/len(summed_data)*100 >= threshold:
            useful_data += usable_data
        segment_duration += len(summed_data)

    print('ObsID: ' + str(obsid))
    print('Segment Length: ' + str(segment_length) + 's')
    print('Useful Data: ' + str(useful_data) + 's ; threshold used: ' + str(threshold) + '%')
    print('Sum of GTIs: ' + str(total_gti) + 's')
    print('Total observation duration: ' + str(obs_duration) + 's')
    print("Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: " + str(useful_data/obs_duration*100))
    print("Percentage of bins with data (over threshold in each segment) over the GTIs: " + str(useful_data/total_gti*100))

    return

#print('0034070101')
#duty_cycle('0034070101',0.0005,100,1,10)
#print('1034090111')
#duty_cycle('1034090111',0.00025,200,1,10)

def duty_cycle_dist(obsid,tbin,segment_length,duty_cycle_bin,threshold):
    """
    To get the distribution of duty cycles over all segments, given an ObsID and
    a desired segment length!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    duty_cycle_bin - binning used to calculate the duty cycle
    threshold - if amount of data in the segment is more than threshold IN PERCENTAGE, use the data
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsdir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(obsdir+'*GTI*'+str(segment_length)+'s*.dat')) #grab the .dat files with a specified segment length

    duty_cycle_array = []
    print('Calculating the duty cycle!')
    tbin = np.float(tbin)
    for i in tqdm(range(len(dat_files))):
        binned_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
        dat_times = np.arange(0,tbin*len(binned_data),tbin)

        duty_cycle_times = np.arange(0,tbin*len(binned_data)+tbin,duty_cycle_bin)
        summed_data, binedges, binnumber = stats.binned_statistic(dat_times,binned_data,statistic='sum',bins=duty_cycle_times)
        # the duty_cycle_frac is NOT a perfect metric!!! There will be data equal to the average counts (because I've binned).
        duty_cycle_frac = len(summed_data[(summed_data!=stats.mode(summed_data)[0][0])&(summed_data>0)])/len(summed_data)
        duty_cycle_array.append(duty_cycle_frac)

    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)

    ax1.hist(duty_cycle_array,bins=100,log=True)
    ax1.set_xlabel('Duty cycle fraction',fontsize=12)
    ax1.set_ylabel('log10(Number of segments)',fontsize=12)

    ax2 = fig.add_subplot(111,sharex=ax1,frameon=False)
    ax2.hist(duty_cycle_array,bins=100)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('Number of segments',fontsize=12)

    plt.title('Segment length = ' + str(segment_length) + 's \n Number of segments here: ' + str(len(dat_files)),fontsize=12)
    plt.savefig(obsdir+str(segment_length)+'s_dutycycle.pdf',dpi=900,format='pdf')
    plt.close()

    return np.array(duty_cycle_array)

def duty_cycle_tE(obsid,tbin,segment_length,PI1,PI2,duty_cycle_bin,threshold):
    """
    To determine the percentage of "data used"/total amount of data. Have two
    types of values:
    1) % of bins (of size duty_cycle_bin) with data over the ENTIRE observation
    2) % of bins with data over the GTIs

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    PI1 - lower energy boundary (in units of PI)
    PI2 - upper energy boundary (in units of PI)
    duty_cycle_bin - binning used to calculate the duty cycle
    threshold - if amount of data in the segment is more than threshold IN PERCENTAGE, use the data
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsdir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(obsdir+'*GTI*'+str(segment_length)+'s*'+str(PI1)+'-'+str(PI2)+'*.dat')) #grab the .dat files with a specified segment length

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data
    obs_duration = gtis[-1][1] - gtis[0][0] #to get the duration of the ENTIRE observation

    total_gti = sum([ (gtis[i][1]-gtis[i][0]) for i in range(len(gtis))])

    useful_data = 0
    segment_duration = 0
    for i in range(len(dat_files)):
        tbin = np.float(tbin)
        binned_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
        dat_times = np.arange(0,tbin*len(binned_data),tbin)

        duty_cycle_times = np.arange(0,tbin*len(binned_data)+tbin,duty_cycle_bin)
        summed_data, binedges, binnumber = stats.binned_statistic(dat_times,binned_data,statistic='sum',bins=duty_cycle_times)
        #plt.plot(duty_cycle_times[:-1],summed_data)
        #plt.title(dat_files[i])
        #plt.show()

#        print(stats.mode(summed_data))
        usable_data = len(summed_data[(summed_data!=stats.mode(summed_data)[0][0])&(summed_data>0)])
        if usable_data/len(summed_data)*100 >= threshold:
            useful_data += usable_data
        segment_duration += len(summed_data)

    print('ObsID: ' + str(obsid) + ' ; PI = ' + str(PI1) + ' - ' + str(PI2))
    print('Segment Length: ' + str(segment_length) + 's')
    print('Useful Data: ' + str(useful_data) + 's ; threshold used: ' + str(threshold) + '%')
    print('Sum of GTIs: ' + str(total_gti) + 's')
    print('Total observation duration: ' + str(obs_duration) + 's')
    print("Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: " + str(useful_data/obs_duration*100))
    print("Percentage of bins with data (over threshold in each segment) over the GTIs: " + str(useful_data/total_gti*100))

    return

#duty_cycle_tE('1200250108',0.00025,200,30,250,1,10)

def duty_cycle_tE_dist(obsid,tbin,segment_length,PI1,PI2,duty_cycle_bin,threshold):
    """
    To get the distribution of duty cycles over all segments, given an ObsID and
    a desired segment length!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    PI1 - lower energy boundary (in units of PI)
    PI2 - upper energy boundary (in units of PI)
    segment_length - length of the individual segments
    duty_cycle_bin - binning used to calculate the duty cycle
    threshold - if amount of data in the segment is more than threshold IN PERCENTAGE, use the data
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsdir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(obsdir+'*GTI*'+str(segment_length)+'s*'+str(PI1)+'-'+str(PI2)+'*.dat')) #grab the .dat files with a specified segment length

    duty_cycle_array = []
    print('Calculating the duty cycle!')
    tbin = np.float(tbin)
    for i in tqdm(range(len(dat_files))):
        binned_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
        dat_times = np.arange(0,tbin*len(binned_data),tbin)

        duty_cycle_times = np.arange(0,tbin*len(binned_data)+tbin,duty_cycle_bin)
        summed_data, binedges, binnumber = stats.binned_statistic(dat_times,binned_data,statistic='sum',bins=duty_cycle_times)
        # the duty_cycle_frac is NOT a perfect metric!!! There will be data equal to the average counts (because I've binned).
        duty_cycle_frac = len(summed_data[(summed_data!=stats.mode(summed_data)[0][0])&(summed_data>0)])/len(summed_data)
        duty_cycle_array.append(duty_cycle_frac)

    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)

    ax1.hist(duty_cycle_array,bins=100,log=True)
    ax1.set_xlabel('Duty cycle fraction',fontsize=12)
    ax1.set_ylabel('log10(Number of segments)',fontsize=12)

    ax2 = fig.add_subplot(111,sharex=ax1,frameon=False)
    ax2.hist(duty_cycle_array,bins=100)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('Number of segments',fontsize=12)

    plt.title('Segment length = ' + str(segment_length) + 's \n Number of segments here: ' + str(len(dat_files)),fontsize=12)
    plt.savefig(obsdir+str(segment_length)+'s_dutycycle.pdf',dpi=900,format='pdf')
    plt.close()

    return np.array(duty_cycle_array)

def compare_segment_lengths(obsid,tbin,segment_lengths,duty_cycle_bin):
    """
    To get the distribution of duty cycles over all segments, given an ObsID and
    a desired segment length! Compare through different thresholds!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_lengths - array of length of the individual segments
    duty_cycle_bin - binning used to calculate the duty cycle
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    thresholds = np.array([0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    store_threshold_y = []
    for i in range(len(seg_lengths)):
        duty_cycle_array = duty_cycle_dist(obsid,tbin,segment_lengths[i],duty_cycle_bin,thresholds[i])
        above5 = len(duty_cycle_array[duty_cycle_array>0.05])/len(duty_cycle_array)
        above10 = len(duty_cycle_array[duty_cycle_array>0.1])/len(duty_cycle_array)
        above20 = len(duty_cycle_array[duty_cycle_array>0.2])/len(duty_cycle_array)
        above30 = len(duty_cycle_array[duty_cycle_array>0.3])/len(duty_cycle_array)
        above40 = len(duty_cycle_array[duty_cycle_array>0.4])/len(duty_cycle_array)
        above50 = len(duty_cycle_array[duty_cycle_array>0.5])/len(duty_cycle_array)
        above60 = len(duty_cycle_array[duty_cycle_array>0.6])/len(duty_cycle_array)
        above70 = len(duty_cycle_array[duty_cycle_array>0.7])/len(duty_cycle_array)
        above80 = len(duty_cycle_array[duty_cycle_array>0.8])/len(duty_cycle_array)
        above90 = len(duty_cycle_array[duty_cycle_array>0.9])/len(duty_cycle_array)
        threshold_y = np.array([above5,above10,above20,above30,above40,above50,above60,above70,above80,above90])
        store_threshold_y.append(threshold_y)

    for i in range(len(store_threshold_y)):
        plt.plot(thresholds,store_threshold_y[i],'x-',lw=0.5,alpha=0.5)
    plt.xlabel('Threshold',fontsize=12)
    plt.ylabel('Fraction of segments (for a given segment length)',fontsize=12)
    plt.legend((tuple(str(seg_lengths[i]) for i in range(len(seg_lengths)))),loc='best')
    plt.savefig('threshold_segmentlength.pdf',dpi=900,format='pdf')
    plt.close()

    subprocess.check_call(['mv','threshold_segmentlength.pdf',Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'])

#seg_lengths = [200,300,500,800,1000,1500,2000]
#for i in range(len(seg_lengths)):
#    duty_cycle_dist('1034090111',0.00025,seg_lengths[i],1,10)

#compare_segment_lengths('1034090111',0.00025,seg_lengths,1)
"""
For 1034090111:

Segment Length: 2000s
Available Data: 29169s
Useful Data: 28445s
Total observation duration: 84669.33687669039s
Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: 33.595397164178046
Percentage of bins with data (over threshold in each segment) over the GTIs: 97.5179128527

Segment Length: 1500s
Available Data: 29726s
Useful Data: 29549s
Total observation duration: 84669.33687669039s
Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: 34.899293049896194
Percentage of bins with data (over threshold in each segment) over the GTIs: 99.4045616632

Segment Length: 1000s
Available Data: 28590s
Useful Data: 28428s
Total observation duration: 84669.33687669039s
Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: 33.57531905724217
Percentage of bins with data (over threshold in each segment) over the GTIs: 99.4333683106

Segment Length: 800s
Available Data: 28709s
Useful Data: 28502s
Total observation duration: 84669.33687669039s
Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: 33.66271787566893
Percentage of bins with data (over threshold in each segment) over the GTIs: 99.278971751

Segment Length: 500s
Available Data: 28875s
Useful Data: 28741s
Total observation duration: 84669.33687669039s
Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: 33.944992437885084
Percentage of bins with data (over threshold in each segment) over the GTIs: 99.5359307359

Segment Length: 300s
Available Data: 28811s
Useful Data: 28780s
Total observation duration: 84669.33687669039s
Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: 33.99105397732621
Percentage of bins with data (over threshold in each segment) over the GTIs: 99.8924022075

Segment Length: 200s
Available Data: 28712s
Useful Data: 28664s
Total observation duration: 84669.33687669039s
Percentage of bins with data (over threshold in each segment) over the ENTIRE observation: 33.8540504241167
Percentage of bins with data (over threshold in each segment) over the GTIs: 99.8328225132

"""
