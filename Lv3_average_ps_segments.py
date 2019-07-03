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
import binary_psr
import glob
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
    gtis = Lv0_call_nicersoft_eventcl.open_fits(obsid,[True,'','','','',''])[2].data
    times = data_dict['TIME']

    #demodulating!
    #times = binary_psr.binary_psr("/Volumes/Samsung_T5/NICERsoft_outputs/J1231-1411.par").demodulate_TOAs(times)

    truncated_times = times-gtis[0][0]
    counts = np.ones(len(times))

    startt = 0
    endt = int(truncated_times[-1])

    print('Binning started.')
    t_bins = np.linspace(startt,endt,(endt-startt)*1/tbin_size+1) #getting an array of time values for the bins
    summed_data, bin_edges, binnumber = stats.binned_statistic(truncated_times,counts,statistic='sum',bins=t_bins) #binning the counts in the data
    print('Binning finished.')

    return t_bins[:-1], summed_data

def presto_dat(obsid,segment_length):
    """
    Obtain the dat files that were generated from PRESTO

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    """
    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(segment_dir + '*.dat'))

    return dat_files

def presto_FFT(obsid,segment_length):
    """
    Obtain the FFT files that were generated from PRESTO

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    """
    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    fft_files = sorted(glob.glob(segment_dir + '*.fft'))

    return fft_files

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
#        plt.plot(t_bins_truncated[:-1],summed_truncated,'bx-')
#        print(len(summed_truncated[summed_truncated>0])/len(summed_truncated)*100)


        if len(summed_truncated[summed_truncated>0])/len(summed_truncated)*100 >= threshold:
            f,ps = Lv2_ps_method.manual(truncated_t,truncated_counts,[False,0,1],[False,1],False,[False,5])

            key = 'segment' + str(i+1).zfill(2)
            segment_dict[key] = (f,ps)
        else:
            pass

    print('Will therefore use ' + str(len(segment_dict)) + ' out of ' + str(len(time_segments)) + ' segments.')

    return segment_dict

test_data = '/Volumes/Samsung_T5/NICERsoft_outputs/1060020113_pipe/ni1060020113_nicersoft_bary.evt'
event = fits.open(test_data)
times = event[1].data['TIME']

times_zero = times - event[2].data[0][0]
counts = np.ones(len(times_zero))

#time_bins = np.arange(0,int(times_zero[-1]),0.00025)
time_bins = np.arange(0,int(times_zero[-1]),1)
summed_threshold,bin_edges,binnumber = stats.binned_statistic(times_zero,counts,statistic='sum',bins=time_bins)
plt.plot(time_bins[:-1],summed_threshold,'rx-')

#segments_FFT('1060020113',['TIME','PI','PI_FAST'],0.00025,1000,20)

def average_ps_presto_segments(obsid,segment_length,threshold):
    """
    Do averaged power spectra from the FFT files that were generated from PRESTO!

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    threshold - if data is under threshold (in percentage), then throw OUT the segment!
    """
    fft_files = presto_FFT(obsid,segment_length) #get the sorted list of .fft files
    dat_files = presto_dat(obsid,segment_length) #get the sorted list of .dat files

    counts_test = np.fromfile(dat_files[0],dtype='<f',count=-1)
    t_test = np.linspace(0,segment_length,len(counts_test))

    t_bins_threshold = np.arange(0,segment_length,1) #1 second bins
    summed_threshold, bin_edges_thres, binnumber_thres = stats.binned_statistic(t_test,counts_test,statistic='sum',bins=t_bins_threshold)

    initial_binned_data = np.fromfile(fft_files[0],dtype='complex64',count=-1) #just to initialize averaged power spectrum
    averaged_ps = np.zeros(len(initial_binned_data))
    t_bins_threshold = np.arange(0,segment_length,1) #1 second bins

    segment_no = 0
    multiplier = [0,11,12,16,17,18,1,21,22,23,24,27,28,29,33,34,38,39,40,43,44,45,46,
                49,50,51,55,56,57,5,60,61,62,66,67,68,6,71,72,73,74,77,78,79,7,82,83,
                84]
    for i in tqdm(range(len(fft_files))):
        dat_file = np.fromfile(dat_files[i],dtype='<f',count=-1)
        t_data = np.linspace(0,segment_length,len(dat_file))+(multiplier[i]*1000)
        #plt.plot(t_data,dat_file,'bx-')
        summed_threshold, bin_edges_thres, binnumber_thres = stats.binned_statistic(t_data,dat_file,statistic='sum',bins=t_bins_threshold)
        binsize = segment_length/len(dat_file)

        tbins_here = t_bins_threshold+(multiplier[i]*1000)
        plt.plot(tbins_here[:-1],summed_threshold,'bx-')
        #bin the dat_file to have the threshold!
#        print(len(summed_threshold),len(summed_threshold[summed_threshold>0])/len(summed_threshold)*100)
        if len(summed_threshold[summed_threshold>0])/len(summed_threshold)*100 >= threshold:
            binned_data = np.fromfile(fft_files[i],dtype='complex64',count=-1)
            averaged_ps += 2/sum(dat_file)*np.abs(binned_data)**2
            segment_no += 1

    freqs = np.fft.fftfreq(averaged_ps.size,binsize)
    N = len(freqs)

    print('We used ' + str(segment_no) + ' out of ' + str(len(fft_files)) + ' possible segments.')
    return freqs[1:int(N/2)], averaged_ps[1:int(N/2)]/segment_no

f,ps = average_ps_presto_segments('1060020113',1000,20)
print(np.mean(ps[f>1000]))
plt.figure(3)
plt.plot(f,ps,'r-')
plt.show()

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
    freqs1 = segment_dict[segment_dict_keys[0]][0]
    ps1 = segment_dict[segment_dict_keys[0]][1]
    plt.title(segment_dict_keys[0],fontsize=12)
    plt.hist(np.log10(ps1[freqs1>1]),bins=100,log=True)
    plt.xlabel('Leahy-normalized power',fontsize=12)
    plt.ylabel('log[N(=P)]',fontsize=12)
    print(np.mean(ps1[freqs1>10]))
    plt.figure(11)
    plt.plot(freqs1,ps1)
    plt.xlim([1,1000])
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)

    plt.figure(12)
    freqs2 = segment_dict[segment_dict_keys[1]][0]
    ps2 = segment_dict[segment_dict_keys[1]][1]
    plt.title(segment_dict_keys[1],fontsize=12)
    plt.hist(np.log10(ps2[freqs2>1]),bins=100,log=True)
    plt.xlabel('Leahy-normalized power',fontsize=12)
    plt.ylabel('log[N(=P)]',fontsize=12)
    print(np.mean(ps2[freqs2>10]))
    plt.figure(13)
    plt.plot(freqs2,ps2)
    plt.xlim([1,1000])
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)

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


if __name__ == "__main__":

    obsid = '1060020113'
    #obsid = '1060060170'
    #obsid = '1200250108' #use threshold of 0.1%...
    par_list = ['TIME','PI']
    tbin_size = 0.00025
    desired_length = 1000 # seconds
    shift = 0.3
    threshold = 20
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
    plt.xlabel('log[Leahy-normalized power]',fontsize=12)
    plt.ylabel('N(>P)')
    plt.figure(18)
    plt.semilogy(ps_bins,N_greaterthanP,'rx')
    plt.xlabel('Leahy-normalized power',fontsize=12)
    plt.ylabel('log[N(>P)]')
    plt.figure(16)
    plt.plot(freqs[freqs>1],ps_to_use,'bx')
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)

    print(np.mean(averaged_ps[freqs>10]))

    plt.figure(2)
    plt.plot(freqs,ps,'r-')
    plt.xlim([1,1000])
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)
    plt.figure(3)
    plt.plot(freqs,averaged_ps,'b-')
    plt.xlim([1,1000])
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Leahy-normalized power',fontsize=12)
    plt.figure(4)
    plt.hist(np.log10(averaged_ps[freqs>1]),bins=100,log=True) #log number of bars!
    plt.xlabel('Leahy-normalized power',fontsize=12)
    plt.ylabel('log[N(=P)]',fontsize=12)

    plt.show()
    #get_pulse_profile(obsid,par_list,tbin_size,freqs_to_try,shift,no_phase_bins)
