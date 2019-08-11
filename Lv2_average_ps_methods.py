#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jul 16 1:48pm 2019

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
import subprocess
import os
import glob
import Lv0_dirs

Lv0_dirs.global_par()

def do_demodulate(obsid,segment_length,par_file):
    """
    Do orbital demodulation on the original events.

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    par_file - orbital parameter file for input into binary_psr
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    eventfiles = sorted(glob.glob(segment_dir + '*.evt')) #get absolute paths of all event FITS files
    for i in range(len(eventfiles)): #for every event file (e.g., for each segment)
        oldfile = eventfiles[i] #old event FITS file
        newfile = eventfiles[i][:-4]+'_demod.evt' #new event FITS file, to be demodulated
        subprocess.check_call(['cp',oldfile,newfile])
        with fits.open(newfile,mode='update') as fitsfile_demod:
            MJDREFI = fitsfile_demod[1].header['MJDREFI'] #integer for MJD reference
            MJDREFF = fitsfile_demod[1].header['MJDREFF'] #float decimal for MJD reference

            times = fitsfile_demod[1].data['TIME'] #original time series
            gtis_start = fitsfile_demod[2].data['START'] #original GTI start times
            gtis_stop = fitsfile_demod[2].data['STOP'] #original GTI end times

            times_MJD = MJDREFI + MJDREFF + times/86400 #converting METs to MJD
            gtis_start_MJD = MJDREFI + MJDREFF + gtis_start/86400 #converting GTIs in METs to MJD
            gtis_stop_MJD = MJDREFI + MJDREFF + gtis_stop/86400 #converting GTIs in METs to MJD

            times_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(times_MJD) #demodulated event times
            gtis_start_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(gtis_start_MJD) #demodulated GTI start times
            gtis_stop_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(gtis_stop_MJD) #demodulated GTI end times

            fitsfile_demod[1].data['TIME'] = (times_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs
            fitsfile_demod[2].data['START'] = (gtis_start_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs
            fitsfile_demod[2].data['STOP'] = (gtis_stop_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs

            fitsfile_demod.flush()

    return

def do_nicerfits2presto(obsid,tbin,segment_length):
    """
    Using nicerfits2presto.py to bin the data, and to convert into PRESTO-readable format.

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments for combining power spectra
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    eventfiles = sorted(glob.glob(segment_dir + '*.evt')) #get absolute paths of all demodulated event FITS files
    print('Now converting NICER event FITS files into the PRESTO-readable binary format!')
    for i in tqdm(range(len(eventfiles))):
        try:
            subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),eventfiles[i]])
        except (ValueError,subprocess.CalledProcessError):
            pass

    obsid_files = glob.glob('*'+obsid+'*')
    for i in range(len(obsid_files)):
        subprocess.check_call(['mv',obsid_files[i],segment_dir])

def edit_inf(obsid,tbin,segment_length):
    """
    Editing the .inf file, as it seems like accelsearch uses some information from the .inf file!
    Mainly need to edit the "Number of bins in the time series".
    This is only for when we make segments by time though!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    inf_files = sorted(glob.glob(segment_dir + '*.inf')) #not the .evt file; some .evt files will be empty

    no_desired_bins = float(segment_length)/float(tbin)

    for i in range(len(inf_files)):
        inf_file = open(inf_files[i],'r')
        contents = inf_file.read()
        contents = contents.split('\n')
        inf_file.close()

        nobins_equal = contents[9].index('=') #find the '=' sign for the "Number of bins..." line)
        newstring = contents[9][:nobins_equal+1] + '  ' + str(int(no_desired_bins)) #replace old line with new line containing updated number of bins!

        inf_file = open(inf_files[i],'w')
        for j in range(len(contents)):
            if j != 9:
                inf_file.write(contents[j]+'\n')
            else:
                inf_file.write(newstring+'\n')
        inf_file.close()

    return

def edit_binary(obsid,tbin,segment_length):
    """
    To pad the binary file so that it will be as long as the desired segment length.
    The value to pad with for each time bin, is the average count rate in THAT segment!
    Jul 10: Do zero-padding instead... so that number of counts is consistent!
    Again, this is only for when we make segments by time!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(segment_dir + '*.dat')) #not that order matters here I think, but just in case
    no_desired_bins = float(segment_length)/float(tbin) #TOTAL number of desired bins for the segment
    for i in tqdm(range(len(dat_files))):
        bins = np.fromfile(dat_files[i],dtype='<f',count=-1) #reads the binary file ; converts to little endian, count=-1 means grab everything

        no_padded = int(no_desired_bins - len(bins)) #number of bins needed to reach the TOTAL number of desired bins
        if no_padded >= 0:
            #padding = np.ones(no_padded,dtype=np.float32)*average_count_rate #generate the array of (averaged) counts needed to pad the original segment
            padding = np.zeros(no_padded,dtype=np.float32) #just in case this is ever needed...
            new_bins = np.array(list(bins) + list(padding))
            new_bins.tofile(dat_files[i]) #don't need to do mv since obsdir already has absolute path to the SSD
        else:
            new_bins = bins[:int(no_desired_bins)] #truncate the original series; say we had a 1000s segment, but
            #nicerfits2presto went up to 1008s, so take that last 8s away because there's no data in it anyways...
            new_bins.tofile(dat_files[i])

    return

def realfft_segment(obsid,segment_length):
    """
    Performing PRESTO's realfft on the binned data (.dat)

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(segment_dir+'*.dat')) #not that order matters here I think, but just in case
    # recall that un-truncated data is "*bary.dat", so "*bary_*.dat" is truncated data!
    logfile = segment_dir + 'realfft.log'

    with open(logfile,'w') as logtextfile:
        for i in range(len(dat_files)):
            logtextfile.write(subprocess.check_output(['realfft',dat_files[i]]))
        logtextfile.close()

    return

def presto_dat(obsid,segment_length,demod):
    """
    Obtain the dat files that were generated from PRESTO

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(segment_dir + '*.dat'))
    demod_files = sorted(glob.glob(segment_dir + '*demod.dat'))
    if demod == True:
        return np.array(demod_files)
    else:
        return np.array([datfile for datfile in dat_files if datfile not in set(demod_files)])

def presto_fft(obsid,segment_length,demod):
    """
    Obtain the FFT files that were generated from PRESTO

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    segment_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    fft_files = sorted(glob.glob(segment_dir + '*.fft'))
    demod_files = sorted(glob.glob(segment_dir + '*demod.fft'))
    if demod == True:
        return np.array(demod_files)
    else:
        return np.array([fftfile for fftfile in fft_files if fftfile not in set(demod_files)])

def segment_threshold(obsid,segment_length,demod,tbin_size,threshold):
    """
    Using the .dat files, rebin them into 1s bins, to weed out the segments below
    some desired threshold. Will return a *list* of *indices*! This is so that I
    can filter out the *sorted* array of .dat and .fft files that are below threshold!

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    tbin_size - size of the time bin
    threshold - if data is under threshold (in percentage), then don't use the segment!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    dat_files = presto_dat(obsid,segment_length,demod)
    rebin_t = np.arange(segment_length+1)*1 #1-second bins

    passed_threshold = []
    print('Now finding the number of segments that can be used...')
    for i in tqdm(range(len(dat_files))):
        dat_file_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
        data_t = np.arange(len(dat_file_data))*tbin_size
        rebin_sum,rebin_edges,rebin_trunc = stats.binned_statistic(data_t,dat_file_data,statistic='sum',bins=rebin_t)
        #print(len(rebin_sum[rebin_sum>0])/len(rebin_sum)*100)
        if len(rebin_sum[rebin_sum>0])/len(rebin_sum)*100 >= threshold:
            passed_threshold.append(i)

    print('Will use ' + str(len(passed_threshold)) + ' out of ' + str(len(dat_files)) + ' segments.')

    return passed_threshold, len(passed_threshold)

def average_ps(obsid,segment_length,demod,tbin_size,threshold,starting_freq,W):
    """
    Given the full list of .dat and .fft files, and the indices where the PRESTO-binned
    data is beyond some threshold, return the averaged power spectrum!

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    tbin_size - size of the time bin
    threshold - if data is under threshold (in percentage), then don't use the segment!
    W - number of consecutive frequency bins to AVERAGE over
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    dat_files = presto_dat(obsid,segment_length,demod) #sorted array of .dat files
    fft_files = presto_fft(obsid,segment_length,demod) #sorted array of .fft files
    passed_threshold,M = segment_threshold(obsid,segment_length,demod,tbin_size,threshold)
    #list of indices where the rebinned .dat files are beyond the threshold

    dat_threshold = dat_files[passed_threshold] #.dat files that passed the threshold
    fft_threshold = fft_files[passed_threshold] #corresponding .fft files that passed the threshold

    freqs = np.fft.fftfreq(int(segment_length/tbin_size),tbin_size)
    N = len(freqs)

    average_ps = np.zeros(int(segment_length/(2*tbin_size)))
    print('Calculating the averaged spectrum...')
    for i in tqdm(range(len(dat_threshold))):
        dat_threshold_data = np.fromfile(dat_threshold[i],dtype='<f',count=-1)
        no_photons = sum(dat_threshold_data)

        fft_threshold_data = np.fromfile(fft_threshold[i],dtype='complex64',count=-1)
        ps = 2/no_photons * np.abs(fft_threshold_data)**2
        average_ps += ps

    print('The mean Leahy power of the latter 90% of the power spectrum is ' + str(np.mean(average_ps[np.int(0.1*len(average_ps)):])/len(passed_threshold)))

    if W == 1:
        f = freqs[1:int(N/2)]
        ps = average_ps[1:]/len(passed_threshold)

        ps_to_use = ps[f>starting_freq]
        ps_bins = np.linspace(min(ps_to_use),max(ps_to_use),1000)
        N_greaterthanP = []
        print('Creating the noise histogram [N(>P)]...')
        for i in tqdm(range(len(ps_bins))):
            array_greaterthan = ps_to_use[ps_to_use>ps_bins[i]]
            N_greaterthanP.append(len(array_greaterthan))

        return f,ps,ps_bins,N_greaterthanP,M

    else:
        pre_f = freqs[1:int(N/2)] #frequency array corresponding to W = 1
        pre_ps = average_ps[1:]/len(passed_threshold) #power array corresponding to W = 1

        consec_f = pre_f[::W] #frequency bins AFTER averaging W consecutive frequency bins
        consec_ps,consec_edges,consec_binnumber = stats.binned_statistic(pre_f,pre_ps,statistic='mean',bins=consec_f)

        f = consec_f[:-1]
        ps = consec_ps

        ps_to_use = ps[f>starting_freq]
        ps_bins = np.linspace(min(ps_to_use),max(ps_to_use),1000)
        N_greaterthanP = []
        print('Creating the noise histogram [N(>P)]...')
        for i in tqdm(range(len(ps_bins))):
            array_greaterthan = ps_to_use[ps_to_use>ps_bins[i]]
            N_greaterthanP.append(len(array_greaterthan))

        return f,ps,ps_bins,N_greaterthanP,M

def noise_hist(obsid,segment_length,demod,tbin_size,threshold,starting_freq,W):
    """
    Given the average spectrum for an ObsID, return the histogram of powers, such
    that you have N(>P). This is for powers corresponding to frequencies larger
    than some starting frequency (perhaps to avoid red noise).

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    tbin_size - size of the time bin
    threshold - if data is under threshold (in percentage), then don't use the segment!
    starting_freq - frequency to start constructing the histogram of powers from
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    f,ps = average_ps(obsid,segment_length,demod,tbin_size,threshold,W)

    ps_to_use = ps[f>starting_freq]
    ps_bins = np.linspace(min(ps_to_use),max(ps_to_use),1000)
    N_greaterthanP = []
    print('Creating the noise histogram [N(>P)]...')
    for i in range(len(ps_bins)):
        array_greaterthan = ps_to_use[ps_to_use>ps_bins[i]]
        N_greaterthanP.append(len(array_greaterthan))

    return ps_bins, N_greaterthanP

if __name__ == "__main__":
    #print(presto_dat('0034070101',100,False))
    #print(presto_fft('0034070101',100,False))
    #print(segment_threshold('0034070101',100,False,0.00025,3))
    f,ps = average_ps('0034070101',100,False,0.00025,10,3)

    plt.figure(1)
    plt.plot(f,ps,'r-')

    plt.figure(2)
    ps_bins, N_greaterthanP = noise_hist('0034070101',100,False,0.00025,10,10,3)
    plt.semilogy(ps_bins,N_greaterthanP,'rx')

    plt.show()
