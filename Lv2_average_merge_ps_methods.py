#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 3:35pm 2019

Getting averaged power spectra from M segments to MERGED OBSERVATIONS, where the data
was pre-processed using NICERsoft!

By MERGED OBSERVATIONS, we mean event files that have been created from merging
event times from multiple ObsIDs!

"""
from __future__ import division, print_function
import numpy as np
from scipy import stats, signal
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
import binary_psr
import subprocess
import glob
import os
import Lv0_dirs,Lv0_nicersoft_evt_filename,Lv2_mkdir

Lv0_dirs.global_par()

def merging(obsids):
    """
    Given a list of ObsIDs, create a file which merges all the rows in the EVENTS
    extension of the FITS files!

    obsids - list (or array) of ObsIDs
    """
    if type(obsids) != list and type(obsids) != np.ndarray:
        raise TypeError("obsids should either be a list or an array!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/'
    no_existing_merged = len(glob.glob(all_merged_dir+'merged*')) #number of existing merged event files/directories
    merged_id = str(no_existing_merged+1).zfill(6) #e.g., if there were 8 merged files, the ID for the next file is 00009
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    Lv2_mkdir.makedir(merged_dir)

    inputfile_string = ''
    for i in range(len(obsids)):
        abs_filename = Lv0_dirs.NICERSOFT_DATADIR + str(obsids[i]) + '_pipe/ni' + str(obsids[i]) + '_nicersoft_bary.evt'
        if i != len(obsids)-1:
            inputfile_string += abs_filename +','
        else:
            inputfile_string += abs_filename

    ##### Writing the ObsIDs into the text file
    merged_text_filename = merged_dir + 'merged' + merged_id + '.txt'
    merged_text = open(merged_text_filename,'w')
    merged_text.write('merged'+merged_id+': '+ ' '.join(obsids))
    merged_text.close()

    subprocess.check_call(['ftmerge',inputfile_string,merged_dir+'merged'+merged_id + '_nicersoft_bary.evt','clobber=YES'])

    return

def merging_GTIs(obsids,merged_id):
    """
    Given a list of ObsIDs and the merged_id, create the final event file, which
    already has all the rows in the EVENTS extension of the FITS files, BUT also
    including the GTI extension of the FITS files!

    obsids - list (or array) of ObsIDs
    merged_id - 6-digit ID for the merged event file
    """
    if type(obsids) != list and type(obsids) != np.ndarray:
        raise TypeError("obsids should either be a list or an array!")
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/'
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'

    inputfile_string = merged_dir+'merged' + merged_id + '_nicersoft_bary.evt[gti],'
    for i in range(1,len(obsids)):
        abs_filename = Lv0_dirs.NICERSOFT_DATADIR + str(obsids[i]) + '_pipe/ni' + str(obsids[i]) + '_nicersoft_bary.evt[gti]'
        if i != len(obsids)-1:
            inputfile_string += abs_filename +','
        else:
            inputfile_string += abs_filename

    subprocess.check_call(['ftmerge',inputfile_string,merged_dir+'merged'+merged_id + '_nicersoft_bary.evt','clobber=YES'])

    return

def get_gti_file(merged_id,segment_length):
    """
    Creating the individual .gti files for my data segments!

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the individual segments for combining power spectra
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_filename = all_merged_dir + 'merged' + merged_id + '/merged' + merged_id + '_nicersoft_bary.evt'
    event = fits.open(merged_filename)
    gtis = event[2].data

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length)
    #array of time values, of size defined by the desired segment length

    gti_folder = all_merged_dir + 'merged' + merged_id + '/' + 'accelsearch_' + str(segment_length).zfill(5) + 's/gtis/'
    Lv2_mkdir.makedir(gti_folder)
    print('Generating the GTI files.')
    for i in tqdm(range(len(segment_times)-1)):
        output_file = gti_folder + str(segment_length).zfill(5) + 's_GTI' + str(i).zfill(4) + '.gti'
        subprocess.check_call(['mkgti.py','--gtiname',output_file,str(segment_times[i]),str(segment_times[i+1])])

    return

def niextract_gti_time(merged_id,segment_length):
    """
    Using niextract-events to get segmented data based on the [segment_length]-length
    GTIs that were created above!

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the individual segments for combining power spectra
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_filename = all_merged_dir + 'merged' + merged_id + '/merged' + merged_id + '_nicersoft_bary.evt'
    event = fits.open(merged_filename)
    gtis = event[2].data

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length)
    #array of time values, of size defined by the desired segment length

    segment_folder = all_merged_dir + 'merged' + merged_id + '/' + 'accelsearch_' + str(segment_length).zfill(5) + '/'
    gti_folder = segment_folder + 'gtis/'
    for i in tqdm(range(len(segment_times)-1)):
        input_file = merged_filename
        output_file = segment_folder + 'merged' + merged_id + '_nicersoft_bary_GTI' + str(i).zfill(4) + '_' + str(segment_length).zfill(5) + 's.evt'

        subprocess.check_call(['niextract-events',input_file,output_file,'timefile='+gti_folder+str(segment_length).zfill(5)+'s_GTI'+str(i).zfill(4)+'.gti'])

    return

def niextract_gti_energy(merged_id,PI1,PI2):
    """
    Using niextract-events to get segmented data based on the energy range

    merged_id - 6-digit ID for the merged event file
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    merged_filename = merged_dir + 'merged' + merged_id + '_nicersoft_bary.evt'
    output_file = energy_trunc_file = all_merged_dir + 'merged' + merged_id + '/merged' + merged_id + '_nicersoft_bary_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '.evt'

    subprocess.check_call(['niextract-events',merged_filename + '[PI='+str(PI1)+':'+str(PI2)+']',energy_trunc_file])

def niextract_gti_time_energy(merged_id,segment_length,PI1,PI2):
    """
    Using niextract-events to get segmented data based on [segment_length]-length
    GTIs that were created above, AND energy range!

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the individual segments for combining power spectra
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_filename = all_merged_dir + 'merged' + merged_id + '/merged' + merged_id + '_nicersoft_bary.evt'
    event = fits.open(merged_filename)
    gtis = event[2].data

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length)
    #array of time values, of size defined by the desired segment length

    segment_folder = all_merged_dir + 'merged' + merged_id + '/' + 'accelsearch_' + str(segment_length).zfill(5) + 's/'
    gti_folder = segment_folder + 'gtis/'
    for i in tqdm(range(len(segment_times)-1)):
        output_filename = segment_folder + 'merged' + merged_id + '_nicersoft_bary_GTI' + str(i).zfill(4) + '_' + str(segment_length).zfill(5) + 's_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '.evt'
        subprocess.check_call(['niextract-events',merged_filename+'[PI='+str(PI1)+':'+str(PI2)+']',output_filename,'timefile='+gti_folder+str(segment_length).zfill(5)+'s_GTI'+str(i).zfill(4)+'.gti'])

    return

def do_nicerfits2presto(merged_id,tbin,segment_length):
    """
    Using nicerfits2presto.py to bin the data, and to convert into PRESTO-readable format.
    I can always move files to different folders to prevent repeats (especially for large files)

    merged_id - 6-digit ID for the merged event file
    tbin - size of the bins in time
    segment_length - length of the individual segments for combining power spectra
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    merged_segment = merged_dir + 'accelsearch_' + str(segment_length).zfill(5) +'s/'

    merged_dir_E_files = glob.glob(merged_dir + '*E*.evt')
    for i in range(len(merged_dir_E_files)):
        if not os.path.isfile(merged_dir_E_files[:-3] + 'dat'):
            try:
                subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),merged_dir_E_files[i]])
            except (ValueError,subprocess.CalledProcessError):
                pass

    merged_segment_files = glob.glob(merged_segment + '*GTI*.evt') #should be sufficient; won't accidentally take in .gti files, and will include GTI AND E files!
    for i in range(len(merged_segment_files)):
        if not os.path.isfile(merged_segment_files[i][:-3] + 'dat'):
            try:
                subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),merged_segment_files[i]])
            except (ValueError,subprocess.CalledProcessError):
                pass

    merged_files = glob.glob('*merged'+merged_id+'*') #so all the inf, dat, and events files!
    for i in range(len(merged_files)):
        subprocess.check_call(['mv',merged_files[i],merged_segment])

    return

def edit_inf(merged_id,tbin,segment_length):
    """
    Editing the .inf file, as it seems like accelsearch uses some information from the .inf file!
    Mainly need to edit the "Number of bins in the time series".
    This is only for when we make segments by time though!

    merged_id - 6-digit ID for the merged event file
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    merged_segment = merged_dir + 'accelsearch_' + str(segment_length).zfill(5) +'s/'
    inf_files = sorted(glob.glob(merged_segment + '*GTI*.inf')) #not the .evt file; some .evt files will be empty

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

    conv_inf_files = glob.glob('*merged'+merged_id+'*.inf')
    for i in range(len(conv_inf_files)):
        subprocess.check_call(['mv',conv_inf_files[i],merged_segment])

    return

def edit_binary(merged_id,tbin,segment_length):
    """
    To pad the binary file so that it will be as long as the desired segment length.
    The value to pad with for each time bin, is the average count rate in THAT segment!
    Jul 10: Do zero-padding instead... so that number of counts is consistent!
    Again, this is only for when we make segments by time!

    merged_id - 6-digit ID for the merged event file
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    merged_segment = merged_dir + 'accelsearch_' + str(segment_length).zfill(5) +'s/'
    dat_files = sorted(glob.glob(merged_segment+'*GTI*.dat')) #not that order matters here I think, but just in case

    print('Editing the binary files now!')
    for i in tqdm(range(len(dat_files))):
        bins = np.fromfile(dat_files[i],dtype='<f',count=-1) #reads the binary file ; converts to little endian, count=-1 means grab everything
#        bins_with_data = len(bins[bins>0]) #number of bins with data (NOT the ones with averaged count rate!)
#        average_count_rate = sum(bins)/len(bins)

        no_desired_bins = float(segment_length)/float(tbin) #TOTAL number of desired bins for the segment
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

    conv_dat_files = glob.glob('*merged'+merged_id+'*.dat')
    for i in range(len(conv_dat_files)):
        subprocess.check_call(['mv',conv_dat_files[i],merged_segment])

    return

def realfft(merged_id,segment_length):
    """
    Performing PRESTO's realfft on the binned data (.dat)

    merged_id - 6-digit ID for the merged event file
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    merged_segment = merged_dir + 'accelsearch_' + str(segment_length).zfill(5) +'s/'

    merged_dir_E_files = glob.glob(merged_dir + '*E*.dat')
    logfile_merged = merged_dir + 'realfft_all.log'
    with open(logfile_merged,'w') as logtextfile:
        for i in range(len(merged_dir_E_files)):
            logtextfile.write(subprocess.check_output(['realfft',merged_dir_E_files[i]]))
        logtextfile.close()

    merged_segment_files = glob.glob(merged_segment + '*GTI*.dat') #should be sufficient; won't accidentally take in .gti files, and will include GTI AND E files!
    logfile_segment = merged_segment + 'realfft_segment.log'
    with open(logfile_segment,'w') as logtextfile:
        for i in range(len(merged_segment_files)):
            subprocess.check_call(['realfft',merged_segment_files[i]])
        logtextfile.close()

    return

def do_demodulate(merged_id,segment_length,par_file):
    """
    Do orbital demodulation on the original events.

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the segments
    par_file - orbital parameter file for input into binary_psr

    Jul 18 - actually, since this is an "average_ps" thing, don't need to
    demodulate the energy truncation event files!
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    merged_segment = merged_dir + 'accelsearch_' + str(segment_length).zfill(5) +'s/'

    eventfiles = sorted(glob.glob(merged_segment + '*.evt')) #get absolute paths of all event FITS files
    for i in range(len(eventfiles)): #for every event file (e.g., for each segment)
        oldfile = eventfiles[i] #old event FITS file
        newfile = eventfiles[i][:-4]+'_demod.evt' #new event FITS file, to be demodulated
        subprocess.check_call(['cp',oldfile,newfile])
        #really to prevent myself from repeating the demodulation multiple times if I run the function again...
        with fits.open(newfile,mode='update') as fitsfile_demod:
            MJDREFI = fitsfile_demod[1].header['MJDREFI'] #integer for MJD reference
            MJDREFF = fitsfile_demod[1].header['MJDREFF'] #float decimal for MJD reference

            times = fitsfile_demod[1].data['TIME'] #original time series
            gtis_start = fitsfile_demod[2].data['START'] #original GTI start times
            gtis_stop = fitsfile_demod[2].data['STOP'] #original GTI end times

            times_MJD = MJDREFI + MJDREFF + times/86400 #converting METs to MJD
            gtis_start_MJD = MJDREFI + MJDREFF + gtis_start/86400 #converting GTIs in METs to MJD
            gtis_stop_MJD = MJDREFI + MJDREFF + gtis_stop/86400 #converting GTIs in METs to MJD

            try:
                times_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(times_MJD) #demodulated event times
                gtis_start_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(gtis_start_MJD) #demodulated GTI start times
                gtis_stop_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(gtis_stop_MJD) #demodulated GTI end times

                fitsfile_demod[1].data['TIME'] = (times_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs
                fitsfile_demod[2].data['START'] = (gtis_start_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs
                fitsfile_demod[2].data['STOP'] = (gtis_stop_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs

                fitsfile_demod.flush()

            except ValueError:
                pass
    return

def presto_dat(merged_id,segment_length,demod,PI1,PI2):
    """
    Obtain the dat files that were generated from PRESTO

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    merged_segment = merged_dir + 'accelsearch_' + str(segment_length).zfill(5) +'s/'

    dat_files = sorted(glob.glob(merged_segment + '*E'+str(PI1).zfill(4)+'-'+str(PI2).zfill(4)+'*.dat'))
    demod_files = sorted(glob.glob(merged_segment + '*E'+str(PI1).zfill(4)+'-'+str(PI2).zfill(4)+'*demod.dat'))
    if demod == True:
        return np.array(demod_files)
    else:
        return np.array([datfile for datfile in dat_files if datfile not in set(demod_files)])

def presto_fft(merged_id,segment_length,demod,PI1,PI2):
    """
    Obtain the FFT files that were generated from PRESTO

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    merged_segment = merged_dir + 'accelsearch_' + str(segment_length).zfill(5) +'s/'

    fft_files = sorted(glob.glob(merged_segment + '*E'+str(PI1).zfill(4)+'-'+str(PI2).zfill(4)+'*.fft'))
    demod_files = sorted(glob.glob(merged_segment + '*E'+str(PI1).zfill(4)+'-'+str(PI2).zfill(4)+'*demod.fft'))
    if demod == True:
        return np.array(demod_files)
    else:
        return np.array([fftfile for fftfile in fft_files if fftfile not in set(demod_files)])

def segment_threshold(merged_id,segment_length,demod,PI1,PI2,tbin_size,threshold):
    """
    Using the .dat files, rebin them into 1s bins, to weed out the segments below
    some desired threshold. Will return a *list* of *indices*! This is so that I
    can filter out the *sorted* array of .dat and .fft files that are below threshold!

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    tbin_size - size of the time bin
    threshold - if data is under threshold (in percentage), then don't use the segment!
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    dat_files = presto_dat(merged_id,segment_length,demod,PI1,PI2)
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

    return passed_threshold

def average_ps(merged_id,segment_length,demod,PI1,PI2,tbin_size,threshold,starting_freq,W):
    """
    Given the full list of .dat and .fft files, and the indices where the PRESTO-binned
    data is beyond some threshold, return the averaged power spectrum!

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    tbin_size - size of the time bin
    threshold - if data is under threshold (in percentage), then don't use the segment!
    W - number of consecutive frequency bins to AVERAGE over
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    dat_files = presto_dat(merged_id,segment_length,demod,PI1,PI2) #sorted array of .dat files
    fft_files = presto_fft(merged_id,segment_length,demod,PI1,PI2) #sorted array of .fft files
    passed_threshold = segment_threshold(merged_id,segment_length,demod,PI1,PI2,tbin_size,threshold)
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

        return f,ps,ps_bins,N_greaterthanP

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

        return f,ps,ps_bins,N_greaterthanP

def noise_hist(merged_id,segment_length,demod,PI1,PI2,tbin_size,threshold,starting_freq,W):
    """
    Given the average spectrum for an ObsID, return the histogram of powers, such
    that you have N(>P). This is for powers corresponding to frequencies larger
    than some starting frequency (perhaps to avoid red noise).

    merged_id - 6-digit ID for the merged event file
    segment_length - length of the segments
    demod - whether we're dealing with demodulated data or not!
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    tbin_size - size of the time bin
    threshold - if data is under threshold (in percentage), then don't use the segment!
    starting_freq - frequency to start constructing the histogram of powers from
    """
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")
    if demod != True and demod != False:
        raise ValueError("demod should either be True or False!")

    f,ps = average_ps(merged_id,segment_length,demod,PI1,PI2,tbin_size,threshold,W)

    ps_to_use = ps[f>starting_freq]
    ps_bins = np.linspace(min(ps_to_use),max(ps_to_use),1000)
    N_greaterthanP = []
    print('Creating the noise histogram [N(>P)]...')
    for i in tqdm(range(len(ps_bins))):
        array_greaterthan = ps_to_use[ps_to_use>ps_bins[i]]
        N_greaterthanP.append(len(array_greaterthan))

    return ps_bins, N_greaterthanP

if __name__ == "__main__":
    obsids = ['2060060363','2060060364','2060060365']
    merged_id = '000001'
    merging(obsids)
    merging_GTIs(obsids,merged_id)
