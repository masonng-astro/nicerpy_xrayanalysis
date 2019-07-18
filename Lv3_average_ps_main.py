#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs Jul 18 11:00am 2019

Time to execute the programs for averaged power spectra!

"""
from __future__ import division, print_function
import numpy as np
import time

import Lv0_dirs
import Lv2_average_ps_methods,Lv2_average_merge_ps_methods

import matplotlib.pyplot as plt

demod = False
merged = True
time_segments = False
time_energy_segments = True

##### For merged = False:
if merged == False:
    obsid = '0034070101' #observation ID
    segment_length = 100 #segment length
    par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1231-1411.par' #parameter file for demodulation
    tbin = 0.00025 #bin size in s
    threshold = 10 #threshold for counts in each segment
    W = 50 #number of consecutive Fourier bins to average over
    starting_freq = 10 #for noise_hist

##### For merged = True:
if merged == True:
    obsids = ['2060060363','2060060364','2060060365'] #list of ObsIDs to combine!
    merged_id = '000001' #need to be very careful that I know what the next one is!
    segment_length = 5000 #segment length
    par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1231-1411.par' #parameter file for demodulation
    PI1 = 30 #lower bound for PI
    PI2 = 200 #upper bound for PI
    tbin = 0.00025 #bin size in s
    threshold = 3 #threshold for counts in each segment
    W = 5000 #number of consecutive Fourier bins to average over
    starting_freq = 10 #for noise_hist

################################################################################

if merged == False:
    if demod == True:
        Lv2_average_ps_methods.do_demodulate(obsid,segment_length,par_file)
    Lv2_average_ps_methods.do_nicerfits2presto(obsid,tbin,segment_length)
    Lv2_average_ps_methods.edit_inf(obsid,tbin,segment_length)
    Lv2_average_ps_methods.edit_binary(obsid,tbin,segment_length)
    Lv2_average_ps_methods.realfft_segment(obsid,segment_length)
    Lv2_average_ps_methods.segment_threshold(obsid,segment_length,demod,tbin,threshold)
    f,ps = Lv2_average_ps_methods.average_ps(obsid,segment_length,demod,tbin,threshold,W)
    ps_bins,N_greaterthanP = Lv2_average_ps_methods.noise_hist(obsid,segment_length,demod,tbin,threshold,starting_freq,W)

    plt.figure(1)
    plt.plot(f,ps,'r-')
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)

    plt.figure(2)
    plt.semilogy(ps_bins,N_greaterthanP,'rx')
    plt.xlabel('Leahy-normalized power',fontsize=12)
    plt.ylabel('log[N(>P)',fontsize=12)

    plt.show()

if merged == True:
    """
    Lv2_average_merge_ps_methods.merging(obsids)
    Lv2_average_merge_ps_methods.merging_GTIs(obsids,merged_id)
    Lv2_average_merge_ps_methods.get_gti_file(merged_id,segment_length)
    if time_segments == True:
        Lv2_average_merge_ps_methods.niextract_gti_time(merged_id,segment_length)
    if time_energy_segments == True:
        Lv2_average_merge_ps_methods.niextract_gti_time_energy(merged_id,segment_length,PI1,PI2)

    if demod == True:
        Lv2_average_merge_ps_methods.do_demodulate(merged_id,segment_length,par_file)

    Lv2_average_merge_ps_methods.do_nicerfits2presto(merged_id,tbin,segment_length)
    Lv2_average_merge_ps_methods.edit_inf(merged_id,tbin,segment_length)
    Lv2_average_merge_ps_methods.edit_binary(merged_id,tbin,segment_length)
    Lv2_average_merge_ps_methods.realfft(merged_id,segment_length)
    """
    f,ps,ps_bins,N_greaterthanP = Lv2_average_merge_ps_methods.average_ps(merged_id,segment_length,demod,PI1,PI2,tbin,threshold,starting_freq,W)
    #ps_bins,N_greaterthanP = Lv2_average_merge_ps_methods.noise_hist(merged_id,segment_length,demod,PI1,PI2,tbin,threshold,starting_freq,W)

    plt.figure(1)
    plt.plot(f,ps,'r-')
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)
    plt.title('Energy range: ' + str(PI1) + ' - ' + str(PI2) + ', W = '+str(W),fontsize=12)

    plt.figure(2)
    plt.semilogy(ps_bins,N_greaterthanP,'rx')
    plt.xlabel('Leahy-normalized power',fontsize=12)
    plt.ylabel('log[N(>P)',fontsize=12)
    plt.title('Energy range: ' + str(PI1) + ' - ' + str(PI2) + ', W = ' + str(W),fontsize=12)

    plt.show()
