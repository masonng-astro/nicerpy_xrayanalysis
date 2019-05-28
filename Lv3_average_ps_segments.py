#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 7:06pm 2019

Getting averaged power spectra from M segments to the whole data, where the data
was pre-processed using NICERsoft!

"""
from __future__ import division, print_function
import numpy as np
from scipy import stats, signal
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
import Lv0_dirs,Lv2_ps_method,Lv2_phase

Lv0_dirs.global_par()

def get_event(obsid,par_list):
    """
    Get a data dictionary from the NICERsoft pipe folders for a desired ObsID -
    data was pre-processed by NICERsoft!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    data_dict = {}
    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    for i in range(len(par_list)):
        data_dict[par_list[i]] = event[1].data[par_list[i]]

    return data_dict

def binned_data(obsid,par_list,tbin_size):
    """
    Get binned (by tbin_size in s) data for a given ObsID - data was pre-processed
    by NICERsoft!
    """
    data_dict = get_event(obsid,par_list)
    times = data_dict['TIME']
    truncated_times = times-times[0]
    counts = np.ones(len(times))

    startt = 0
    endt = int(truncated_times[-1])

    print('Binning started.')
    t_bins = np.linspace(startt,endt,(endt-startt)*1/tbin_size+1) #getting an array of time values for the bins
    summed_data, bin_edges, binnumber = stats.binned_statistic(truncated_times,counts,statistic='sum',bins=t_bins) #binning the counts in the data
    print('Binning finished.')

    return t_bins[:-1], summed_data
"""
event = Lv0_dirs.NICERSOFT_DATADIR + '1034090111_pipe/GTI_1000s/ni1034090111_nicersoft_bary_GTI1_1000s.evt'
event = fits.open(event)
times = event[1].data['TIME']
counts = np.ones(len(times))

truncated_times = times-times[0]
startt = 0
endt = int(truncated_times[-1])

t_bins = np.linspace(startt,endt,(endt-startt)*1/0.00025+1)
summed_data,bin_edges,binnumber = stats.binned_statistic(truncated_times,counts,statistic='sum',bins=t_bins)

f,ps = Lv2_ps_method.manual(t_bins,summed_data,[False,0,1],[False,1],False,[False,5])

plt.figure(1)
plt.plot(f,ps)
"""

def segments_FFT(obsid,par_list,tbin_size,desired_length):
    """
    Obtain the dictionary of time segments + corresponding counts, and do the manual Fourier transform

    desire_length - desired length of the segments
    """
    segment_dict = {}

    t_bins, summed_data = binned_data(obsid,par_list,tbin_size)
    time_segments = np.arange(0,t_bins[-1],desired_length)
    for i in tqdm(range(len(time_segments)-1)):
        truncated_t = t_bins[(t_bins>=time_segments[i])&(t_bins<=time_segments[i+1])]
        truncated_counts = summed_data[(t_bins>=time_segments[i])&(t_bins<=time_segments[i+1])]

        f,ps = Lv2_ps_method.manual(truncated_t,truncated_counts,[False,0,1],[False,1],False,[False,5])

        key = 'segment' + str(i+1)
        segment_dict[key] = (f,ps)

    return segment_dict

def average_ps_segments(obsid,par_list,tbin_size,desired_length):
    """
    Make sure I call barycentered data!
    """
    segment_dict = segments_FFT(obsid,par_list,tbin_size,desired_length)
    power_spectra = np.zeros(len(segment_dict['segment1'][0])) #initialize power spectra
    for i in range(len(segment_dict)):
        key = 'segment' + str(i+1)
        power_spectra = power_spectra + segment_dict[key][1]

    freqs = segment_dict['segment1'][0]
    averaged_ps = power_spectra/len(segment_dict)

    return freqs,power_spectra,averaged_ps

################################################################################
### PULSE PROFILE

def get_pulse_profile(obsid,par_list,tbin_size,f_pulse,shift,no_phase_bins):
    """
    Get pulse profile given binned data and frequency to fold

    obsid - Observation ID of the object of interest (10-digit str)
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    f_pulse - the frequency of the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how the pulse profile is 'displaced'.
    no_phase_bins - number of phase bins desired
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    t_bins, summed_data = binned_data(obsid,par_list,tbin_size)
    print("Making plots now.")
    if type(f_pulse) == np.ndarray or type(f_pulse) == list:
        for i in range(len(f_pulse)):
            print(f_pulse[i])
            phase, phase_bins, summed_profile = Lv2_phase.pulse_profile(f_pulse[i],t_bins,summed_data,shift,no_phase_bins)
            print(len(phase_bins),len(summed_profile))
            plt.figure(i)
            plt.title('Frequency: ' + str(f_pulse[i]) + ' Hz',fontsize=12)
            plt.xlabel('Phase',fontsize=12)
            plt.ylabel('Counts/s (binned by 0.0005s)',fontsize=12)
            plt.step(phase_bins[:-1],summed_profile)
            plt.show()
    else:
        phase, phase_bins, summed_profile = Lv2_phase.pulse_profile(f_pulse,t_bins,summed_data,shift,no_phase_bins)
        plt.step(phase_bins[:-1],summed_profile)
        plt.title('Frequency: ' + str(f_pulse) + ' Hz',fontsize=12)
        plt.xlabel('Phase',fontsize=12)
        plt.ylabel('Counts/s (binned by 0.0005s)',fontsize=12)
        plt.show()

freqs_to_try = np.array([63.8638, 367.219, 306.683, 427.271, 977.339, 918.264, 982.823, #THIS LINE: ALL FROM AVERAGED POWER SPECTRA
                         210.651, 149.450, 115.677, 86.642, 315.014, 925.222, 203.422, 139.321, 262.021, 79.801, #THIS LINE: ALL FROM BREAKING UP ENERGY RANGE
                         156.290, 370.349]) #THIS LINE: ALL FROM BREAKING UP GTIs. They're from non-flare regions!

obsid = '1034090111'
par_list = ['TIME','PI']
tbin_size = 0.00025
desired_length = 1000 # seconds
shift = 0.3
no_phase_bins = 50

freqs,ps,averaged_ps = average_ps_segments(obsid,par_list,tbin_size,desired_length)

plt.figure(2)
plt.plot(freqs,ps,'r-')
plt.figure(3)
plt.plot(freqs,averaged_ps,'b-')

plt.show()
#get_pulse_profile(obsid,par_list,tbin_size,freqs_to_try,shift,no_phase_bins)
