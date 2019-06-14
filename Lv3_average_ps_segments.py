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
import Lv0_dirs,Lv0_call_nicersoft_eventcl,Lv2_ps_method,Lv2_phase

Lv0_dirs.global_par()

def binned_data(obsid,par_list,tbin_size):
    """
    Get binned (by tbin_size in s) data for a given ObsID - data was pre-processed
    by NICERsoft!

    obsid - Observation ID of the object of interest (10-digit str)
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - size of the time bin
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    data_dict = Lv0_call_nicersoft_eventcl.get_eventcl(obsid,[True,'','','','',''],par_list)
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

def segments_FFT(obsid,par_list,tbin_size,desired_length,threshold):
    """
    Obtain the dictionary of time segments + corresponding counts, and do the manual Fourier transform

    obsid - Observation ID of the object of interest (10-digit str)
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - size of the time bin
    desired_length - desired length of the segments
    threshold - if data is under threshold (in percentage), then throw OUT the segment!
    """
    segment_dict = {}

    t_bins, summed_data = binned_data(obsid,par_list,tbin_size)
    time_segments = np.arange(0,t_bins[-1],desired_length)
    for i in tqdm(range(len(time_segments)-1)):
        truncated_t = t_bins[(t_bins>=time_segments[i])&(t_bins<=time_segments[i+1])]
        truncated_counts = summed_data[(t_bins>=time_segments[i])&(t_bins<=time_segments[i+1])]

        ### applying a threshold
        t_bins_truncated = np.arange(truncated_t[0],truncated_t[-1],1) #should be 1 second!
        summed_truncated, bin_edges_trunc, binnumber_trunc = stats.binned_statistic(truncated_t,truncated_counts,statistic='sum',bins=t_bins_truncated)

        if len(summed_truncated[summed_truncated>0])/len(summed_truncated)*100 >= threshold:
            f,ps = Lv2_ps_method.manual(truncated_t,truncated_counts,[False,0,1],[False,1],False,[False,5])

            key = 'segment' + str(i+1).zfill(2)
            segment_dict[key] = (f,ps)
        else:
            pass

    print('Will therefore use ' + str(len(segment_dict)) + ' out of ' + str(len(time_segments)) + ' segments.')

    return segment_dict

def average_ps_segments(obsid,par_list,tbin_size,desired_length,threshold):
    """
    Make sure I call barycentered data! Does the averaging of the power spectra

    obsid - Observation ID of the object of interest (10-digit str)
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - size of the time bin
    desired_length - desired length of the segments
    threshold - if data is under threshold (in percentage), then throw OUT the segment!
    """
    segment_dict = segments_FFT(obsid,par_list,tbin_size,desired_length,threshold)
    segment_dict_keys = sorted(segment_dict.keys())
    power_spectra = np.zeros(len(segment_dict[segment_dict_keys[0]][0])) #initialize power spectra

    plt.figure(10)
    freqs1 = segment_dict[segment_dict_keys[2]][0]
    ps1 = segment_dict[segment_dict_keys[2]][1]
    plt.title(segment_dict_keys[2],fontsize=12)
    plt.hist(np.log10(ps1[freqs1>1]),bins=100,log=True)
    plt.xlabel('Leahy-normalized power')
    print(np.mean(ps1[freqs1>10]))
    plt.figure(11)
    plt.plot(freqs1,ps1)
    plt.xlim([1,1000])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')

    plt.figure(12)
    freqs2 = segment_dict[segment_dict_keys[5]][0]
    ps2 = segment_dict[segment_dict_keys[5]][1]
    plt.title(segment_dict_keys[5],fontsize=12)
    plt.hist(np.log10(ps2[freqs2>1]),bins=100,log=True)
    plt.xlabel('Leahy-normalized power')
    print(np.mean(ps2[freqs2>10]))
    plt.figure(13)
    plt.plot(freqs2,ps2)
    plt.xlim([1,1000])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')

    for i in range(len(segment_dict)):
        key = segment_dict_keys[i]
        power_spectra = power_spectra + segment_dict[key][1]

    freqs = segment_dict[segment_dict_keys[0]][0] #the frequencies are all the same, so use first segment's
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
#obsid = '1200250108' #use threshold of 0.1%...
par_list = ['TIME','PI']
tbin_size = 0.0005
desired_length = 1000 # seconds
shift = 0.3
threshold = 10
no_phase_bins = 50

freqs,ps,averaged_ps = average_ps_segments(obsid,par_list,tbin_size,desired_length,threshold)
ps_to_use = averaged_ps[freqs>1]

ps_bins = np.linspace(min(ps_to_use),max(ps_to_use),10000)
N_greaterthanP = []
for i in range(len(ps_bins)):
    array_greaterthan = ps_to_use[ps_to_use>ps_bins[i]]
    N_greaterthanP.append(len(array_greaterthan))

plt.figure(15)
#plt.bar(ps_bins,np.log10(N_greaterthanP),width=ps_bins[1]-ps_bins[0])
plt.semilogx(ps_bins,N_greaterthanP,'rx')
plt.figure(18)
plt.semilogy(ps_bins,N_greaterthanP,'rx')
plt.figure(16)
plt.plot(freqs[freqs>1],ps_to_use,'bx')

#print(np.mean(averaged_ps[freqs>10]))

plt.figure(2)
plt.plot(freqs,ps,'r-')
plt.xlim([1,1000])
plt.figure(3)
plt.plot(freqs,averaged_ps,'b-')
plt.xlim([1,1000])
plt.figure(4)
plt.hist(np.log10(averaged_ps[freqs>1]),bins=100,log=True) #log number of bars!

plt.show()
#get_pulse_profile(obsid,par_list,tbin_size,freqs_to_try,shift,no_phase_bins)
