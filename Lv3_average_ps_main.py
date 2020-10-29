#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thurs Jul 18 11:00am 2019

Time to execute the programs for averaged power spectra!

1/9/2020 - I was so confused before as to why I didn't have niextract for merged = False,
then I realized it was because I probably did the segments in presto_subroutines
(or presto_segments in the previous version). Will put it in this time because of
the integration with Lv3_quicklook.py!

"""
from __future__ import division, print_function
import numpy as np
import time
import glob
import subprocess
from tqdm import tqdm

import Lv0_dirs
import Lv2_average_ps_methods,Lv2_presto_subroutines,Lv2_merging_events
import Lv3_detection_level

import matplotlib.pyplot as plt

Lv0_dirs.global_par()

demod = True
merged = True
preprocessing = False
time_segments = False
time_energy_segments = False

##### For merged = False:
if merged == False:
    #eventfile = Lv0_dirs.NICERSOFT_DATADIR + '1034070101_pipe/ni1034070101_nicersoft_bary.evt'
    eventfile = Lv0_dirs.NICER_DATADIR + '/rxj0209/rxj0209kgfilt_bary.evt'
    segment_length = 10000 #segment length
    PI1 = 30 #lower bound for PI
    PI2 = 200 #upper bound for PI
    par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1231-1411.par' #parameter file for demodulation
    tbin = 0.025 #bin size in s
    N = Lv3_detection_level.N_trials(tbin,segment_length)
    threshold = 10 #threshold for counts in each segment
    W = 1 #number of consecutive Fourier bins to average over
    starting_freq = 1 #for noise_hist
    mode = 't'
    xlims = np.array([])
    plot_mode = 'show'

##### For merged = True:
if merged == True:
    obsids = ['20600603'+str(i) for i in range(61,66)]

    merged_id = '000013' #need to be very careful that I know what the next one is!
    eventfile = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/merged' + merged_id + '/merged' + merged_id + '_nicersoft_bary.evt'
    segment_length = 10000 #segment length
    mode = 't'
    par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1231-1411.par' #parameter file for demodulation
    PI1 = 30 #lower bound for PI
    PI2 = 200 #upper bound for PI
    tbin = 0.001 #bin size in s
    N = Lv3_detection_level.N_trials(tbin,segment_length)
    threshold = 2 #threshold for counts in each segment
    W = 1 #number of consecutive Fourier bins to average over
    starting_freq = 10 #for noise_hist
    xlims = np.array([])
    plot_mode = 'show'

################################################################################

if merged == False:
    if preprocessing == True:
        if time_segments == True or time_energy_segments == True:
            Lv2_presto_subroutines.get_gti_file(eventfile,segment_length)
        if time_segments == True:
            Lv2_presto_subroutines.niextract_gti_time(eventfile,segment_length)
        if time_energy_segments == True:
            Lv2_presto_subroutines.niextract_gti_time_energy(eventfile,segment_length,PI1,PI2)

        if demod == True:
            Lv2_average_ps_methods.do_demodulate(eventfile,segment_length,mode,par_file)

        Lv2_average_ps_methods.do_nicerfits2presto(eventfile,tbin,segment_length)
        Lv2_average_ps_methods.edit_inf(eventfile,tbin,segment_length)
        Lv2_average_ps_methods.edit_binary(eventfile,tbin,segment_length)
        Lv2_average_ps_methods.realfft(eventfile,segment_length)

    Lv2_average_ps_methods.plotting(eventfile,segment_length,demod,tbin,threshold,PI1,PI2,starting_freq,W,N,xlims,plot_mode)

if merged == True:
    if preprocessing == True:
        Lv2_merging_events.merging(obsids,merged_id)
        if time_segments == True or time_energy_segments == True:
            Lv2_presto_subroutines.get_gti_file(eventfile,segment_length)
        if time_segments == True:
            Lv2_presto_subroutines.niextract_gti_time(eventfile,segment_length)
        if time_energy_segments == True:
            Lv2_presto_subroutines.niextract_gti_time_energy(eventfile,segment_length,PI1,PI2)

        if demod == True:
            Lv2_average_ps_methods.do_demodulate(eventfile,segment_length,mode,par_file)

        Lv2_average_ps_methods.do_nicerfits2presto(eventfile,tbin,segment_length)
        Lv2_average_ps_methods.edit_inf(eventfile,tbin,segment_length)
        Lv2_average_ps_methods.edit_binary(eventfile,tbin,segment_length)
        Lv2_average_ps_methods.realfft(eventfile,segment_length)

    Lv2_average_ps_methods.plotting(eventfile,segment_length,demod,tbin,threshold,PI1,PI2,starting_freq,W,N,xlims,plot_mode)
