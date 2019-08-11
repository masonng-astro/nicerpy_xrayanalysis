#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs Jul 18 11:00am 2019

Time to execute the programs for averaged power spectra!

"""
from __future__ import division, print_function
import numpy as np
import time
import glob
import subprocess
from tqdm import tqdm

import Lv0_dirs
import Lv2_average_ps_methods,Lv2_average_merge_ps_methods
import Lv3_detection_level

import matplotlib.pyplot as plt

demod = True
merged = True
preprocessing = True
time_segments = False
time_energy_segments = True

##### For merged = False:
if merged == False:
    obsid = '0034070101' #observation ID
    segment_length = 100 #segment length
    par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1231-1411.par' #parameter file for demodulation
    tbin = 0.00025 #bin size in s
    N = Lv3_detection_level.N_trials(tbin,segment_length)
    threshold = 20 #threshold for counts in each segment
    W = 2 #number of consecutive Fourier bins to average over
    starting_freq = 10 #for noise_hist

##### For merged = True:
if merged == True:
    obsids = ['20600603'+str(i) for i in range(58,72)]
    #obsids = ['10301801'+str(i) for i in range(49,58)]

    merged_id = '000006' #need to be very careful that I know what the next one is!
    segment_length = 500 #segment length
    par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1231-1411.par' #parameter file for demodulation
    PI1 = 30 #lower bound for PI
    PI2 = 200 #upper bound for PI
    tbin = 0.001 #bin size in s
    N = Lv3_detection_level.N_trials(tbin,1)*100
    threshold = 20 #threshold for counts in each segment
    W = 1 #number of consecutive Fourier bins to average over
    starting_freq = 10 #for noise_hist

################################################################################

if merged == False:
    if preprocessing == True:
        if demod == True:
            Lv2_average_ps_methods.do_demodulate(obsid,segment_length,par_file)
        Lv2_average_ps_methods.do_nicerfits2presto(obsid,tbin,segment_length)
        Lv2_average_ps_methods.edit_inf(obsid,tbin,segment_length)
        Lv2_average_ps_methods.edit_binary(obsid,tbin,segment_length)
        Lv2_average_ps_methods.realfft_segment(obsid,segment_length)

    f,ps,ps_bins,N_greaterthanP,M = Lv2_average_ps_methods.average_ps(obsid,segment_length,demod,tbin,threshold,starting_freq,W)

    power_required_3 = Lv3_detection_level.power_for_sigma(3,N,M,W) #power required for significance
    power_required_4 = Lv3_detection_level.power_for_sigma(4,N,M,W) #power required for significance

    plt.figure(1)
    plt.plot(f,ps,'rx-')
    plt.axhline(y=power_required_3,lw=0.8,alpha=0.5,color='b')
    plt.axhline(y=power_required_4,lw=0.8,alpha=0.5,color='k')
    plt.xlim([0.3,0.5])
    plt.ylim([0,800])
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)
    plt.title('W = '+ str(W) + ', Threshold = '+str(threshold) + '%' + '\n' + 'Segment Length: ' + str(segment_length) + 's, No. Segments = ' + str(M) + '\n' + 'Demodulated: ' + str(demod),fontsize=12)
    plt.legend(('Power Spectrum','3 sigma','4 sigma'),loc='best')

    plt.figure(2)
    plt.semilogy(ps_bins,N_greaterthanP,'rx')
    plt.xlabel('Leahy-normalized power',fontsize=12)
    plt.ylabel('log[N(>P)]',fontsize=12)
    plt.title('W = ' + str(W),fontsize=12)

    plt.show()

if merged == True:
    if preprocessing == True:
        Lv2_average_merge_ps_methods.merging(obsids)
        Lv2_average_merge_ps_methods.merging_GTIs(obsids,merged_id)
        Lv2_average_merge_ps_methods.get_gti_file(merged_id,segment_length)
        if time_segments == True:
            Lv2_average_merge_ps_methods.niextract_gti_time(merged_id,segment_length)
        if time_energy_segments == True:
            Lv2_average_merge_ps_methods.niextract_gti_time_energy(merged_id,segment_length,PI1,PI2)

        if demod == True:
            Lv2_average_merge_ps_methods.do_demodulate(merged_id,segment_length,par_file,PI2)

        Lv2_average_merge_ps_methods.do_nicerfits2presto(merged_id,tbin,segment_length)
        Lv2_average_merge_ps_methods.edit_inf(merged_id,tbin,segment_length)
        Lv2_average_merge_ps_methods.edit_binary(merged_id,tbin,segment_length)
        Lv2_average_merge_ps_methods.realfft(merged_id,segment_length)

    f,ps,ps_bins,N_greaterthanP,M = Lv2_average_merge_ps_methods.average_ps(merged_id,segment_length,demod,PI1,PI2,tbin,threshold,starting_freq,W)

    power_required_3 = Lv3_detection_level.power_for_sigma(3,N,M,W) #power required for significance
    power_required_4 = Lv3_detection_level.power_for_sigma(4,N,M,W) #power required for significance

    plt.figure(1)
    plt.plot(f,ps,'rx-')
    plt.axhline(y=power_required_3,lw=0.8,alpha=0.5,color='b')
    plt.axhline(y=power_required_4,lw=0.8,alpha=0.5,color='k')
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)
    #plt.xlim([621.5,622.5])
    #plt.ylim([1,4])
    plt.axvline(x=271.453,lw=0.5,alpha=0.5)
    plt.title('W = ' + str(W) + ', Threshold = ' + str(threshold) + '%' + '\n' + 'Segment Length: ' + str(segment_length) + 's, No. Segments = ' + str(M) + '\n' + 'Demodulated: ' + str(demod) + ' ; St.D = ' + str(np.std(ps)), fontsize=12)
    plt.legend(('Power Spectrum','3 sigma','4 sigma'),loc='best')
    #pngname = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000005/W_dir/W_300.png'
    #plt.savefig(pngname,dpi=900)
    #plt.close()

    plt.figure(2)
    plt.semilogy(ps_bins,N_greaterthanP,'rx')
    plt.xlabel('Leahy-normalized power',fontsize=12)
    plt.ylabel('log[N(>P)]',fontsize=12)
    plt.title('Energy range: ' + str(PI1) + ' - ' + str(PI2) + ', W = ' + str(W),fontsize=12)

    plt.show()
