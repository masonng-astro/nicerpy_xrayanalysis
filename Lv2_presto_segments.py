#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 6:22pm 2019

Program for doing mkgti, niextract-events, nicerfits2presto, padding, calculate duty cycle,
realfft, accelsearch, prepfold, and ps2pdf! This is for when we want to split
the original time series up into segments (whether by energy or time)
Use Lv2_presto_all.py instead if you want to do the analysis for the whole time series.

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs
from tqdm import tqdm
import os
from os.path import relpath
import time
import subprocess
import glob

Lv0_dirs.global_par()

powers_of_two = [2**i for i in range(31)] # for deciding number of FFT bins

def get_gti_file(obsid,segment_length):
    """
    Creating the individual .gti files for my data segments!

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the individual segments for combining power spectra
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end,segment_length) #array of time values, starting

    #from the observation start time until the observation end time, in steps of 1000s.
    #This means that you'll get, for a 4392s observation, np.array([0,1000,2000,3000,4000])!
    #We'd lose 392s worth of counts, but it can't be used in combining the power spectra anyways.

    output_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    for i in tqdm(range(len(segment_times)-1)):
        output_file = output_dir + str(segment_length) + 's_GTI' + str(i) + '.gti'
        subprocess.check_call(['mkgti.py','--gtiname',output_file,str(segment_times[i]),str(segment_times[i+1])] )
        #no need to do 'startmet=', 'stopmet=', but make sure I do it in the right order!

    return

def niextract_gti_time(obsid,segment_length):
    """
    Using niextract-events to get segmented data based on the [segment_length]-length
    GTIs that were created above!

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the individual segments for combining power spectra
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end,segment_length) #array of time values, starting
    #from the observation start time until the observation end time, in steps of 1000s.
    #This means that you'll get, for a 4392s observation, np.array([0,1000,2000,3000,4000])!
    #We'd lose 392s worth of counts, but it can't be used in combining the power spectra anyways.

    working_dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/'
    for i in tqdm(range(len(segment_times)-1)):
        inputfile = working_dir+'ni' + obsid + '_nicersoft_bary.evt'
        outputfile = working_dir + 'ni' + obsid + '_nicersoft_bary_GTI'+str(i)+'_'+str(segment_length)+'s.evt'

        subprocess.check_call(['niextract-events',inputfile,outputfile,'timefile='+working_dir+str(segment_length)+'s_GTI'+str(i)+'.gti'])

    return

def niextract_gti_energy(obsid,PI1,PI2):
    """
    Using niextract-events to get segmented data based on the energy range

    obsid - Observation ID of the object of interest (10-digit str)
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'

    working_dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/'
    inputfile = working_dir + 'ni' + obsid + '_nicersoft_bary.evt'
    outputfile = working_dir + 'ni' + obsid + '_nicersoft_bary_E'+str(PI1)+'-'+str(PI2)+'.evt'

    subprocess.check_call(['niextract-events',inputfile+'[PI='+str(PI1)+':'+str(PI2)+']',outputfile])

    return

#niextract_gti_energy('1200250108',800,1200)

def niextract_gti_time_energy(obsid,segment_length,PI1,PI2):
    """
    Using niextract-events to get segmented data based on [segment_length]-length
    GTIs that were created above, AND energy range!

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the individual segments for combining power spectra
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end,segment_length) #array of time values, starting
    #from the observation start time until the observation end time, in steps of 1000s.
    #This means that you'll get, for a 4392s observation, np.array([0,1000,2000,3000,4000])!
    #We'd lose 392s worth of counts, but it can't be used in combining the power spectra anyways.

#    get_gti_file(obsid,segment_length) #to get GTI files...

    working_dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/'
    inputfile = working_dir + 'ni' + obsid + '_nicersoft_bary.evt'

    working_dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/'
    for i in tqdm(range(len(segment_times)-1)):
        outputfile = working_dir + 'ni' + obsid + '_nicersoft_bary_GTI'+str(i)+'_'+str(segment_length)+'s_' + 'E'+str(PI1)+'-'+str(PI2)+'.evt'
        subprocess.check_call(['niextract-events',inputfile+'[PI='+str(PI1)+':'+str(PI2)+']',outputfile,'timefile='+working_dir+str(segment_length)+'s_GTI'+str(i)+'.gti'])

    return

def do_nicerfits2presto(obsid,tbin,segment_length):
    """
    Using nicerfits2presto.py to bin the data, and to convert into PRESTO-readable format.
    I can always move files to different folders to prevent repeats (especially for large files)

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments for combining power spectra
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obs_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    E_files = glob.glob(obs_dir+'*E*.evt')

    time_files = glob.glob(obs_dir+'*GTI*'+str(segment_length)+'s.evt')

    if len(E_files) != 0: #if E_files is not empty
        for i in tqdm(range(len(E_files))):
#            if E_files[i] not in E_files_done:
            try:
                subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),E_files[i]])
            except (ValueError,subprocess.CalledProcessError):
                pass
    if len(time_files) != 0:
        for i in tqdm(range(len(time_files))):
            if time_files[i] not in E_files: #to prevent the files truncated by E AND time to be processed AGAIN.
                try:
                    subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),time_files[i]])
                except (ValueError,subprocess.CalledProcessError):
                    pass

    ##### will probably remove once I know to either NOT work in nicerpy_xrayanalysis,
    ##### and/or install my packages such that I can run the scripts from anywhere
    ##### in the terminal!
    obsid_files = glob.glob('*'+obsid+'*')
    for i in range(len(obsid_files)):
        subprocess.check_call(['mv',obsid_files[i],obs_dir])

    return

def edit_inf(obsid,tbin,segment_length):
    """
    Editing the .inf file, as it seems like accelsearch uses some information from the .inf file!
    Mainly need to edit the "Number of bins in the time series".
    This is only for when we make segments by time though!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    obs_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    inf_files = sorted(glob.glob(obs_dir + '*GTI*' + str(segment_length)+'s*.inf')) #not the .evt file; some .evt files will be empty

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
    Again, this is only for when we make segments by time!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsdir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    dat_files = sorted(glob.glob(obsdir+'*GTI*' + str(segment_length) + 's*.dat')) #not that order matters here I think, but just in case
    for i in tqdm(range(len(dat_files))):
        bins = np.fromfile(dat_files[i],dtype='<f',count=-1) #reads the binary file ; converts to little endian, count=-1 means grab everything
        bins_with_data = len(bins[bins>0]) #number of bins with data (NOT the ones with averaged count rate!)
        average_count_rate = sum(bins)/len(bins)

        no_desired_bins = float(segment_length)/float(tbin) #TOTAL number of desired bins for the segment
        no_padded = int(no_desired_bins - len(bins)) #number of bins needed to reach the TOTAL number of desired bins
        if no_padded >= 0:
            padding = np.ones(no_padded,dtype=np.float32)*average_count_rate #generate the array of (averaged) counts needed to pad the original segment
            #padding = np.zeros(no_padded,dtype=np.float32) #just in case this is ever needed...
            new_bins = np.array(list(bins) + list(padding))
            new_bins.tofile(dat_files[i]) #don't need to do mv since obsdir already has absolute path to the SSD
        else:
            new_bins = bins[:int(no_desired_bins)] #truncate the original series; say we had a 1000s segment, but
            #nicerfits2presto went up to 1008s, so take that last 8s away because there's no data in it anyways...
            new_bins.tofile(dat_files[i])

    return

def realfft(obsid):
    """
    Performing PRESTO's realfft on the binned data (.dat)

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    dat_files = sorted(glob.glob(nicersoft_output_folder+'*bary_*.dat')) #not that order matters here I think, but just in case
    # recall that un-truncated data is "*bary.dat", so "*bary_*.dat" is truncated data!
    logfile = nicersoft_output_folder + 'realfft.log'

    with open(logfile,'w') as logtextfile:
        for i in range(len(dat_files)):
            logtextfile.write(subprocess.check_output(['realfft',dat_files[i]]))
        logtextfile.close()

    return

def realfft_segment(obsid,segment_length):
    """
    Performing PRESTO's realfft on the binned data (.dat)

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    dat_files = sorted(glob.glob(nicersoft_output_folder+'*bary_*'+str(segment_length)+'s*.dat')) #not that order matters here I think, but just in case
    # recall that un-truncated data is "*bary.dat", so "*bary_*.dat" is truncated data!
    logfile = nicersoft_output_folder + 'realfft.log'

    with open(logfile,'w') as logtextfile:
        for i in range(len(dat_files)):
            logtextfile.write(subprocess.check_output(['realfft',dat_files[i]]))
        logtextfile.close()

    return

def accelsearch(obsid,flags):
    """
    Performing PRESTO's accelsearch on the FFT data (.fft)

    obsid - Observation ID of the object of interest (10-digit str)
    flags - a LIST of input flags for accelsearch
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    fft_files = sorted(glob.glob(nicersoft_output_folder+'*bary_*.fft')) #not that order matters here I think, but just in case
    logfile = nicersoft_output_folder + 'accelsearch_segments.log'

    with open(logfile,'w') as logtextfile:
        for i in range(len(fft_files)):
            logtextfile.write(subprocess.check_output(['accelsearch']+flags+[fft_files[i]]))
        logtextfile.close()

    return

def accelsearch_segment(obsid,flags,segment_length):
    """
    Performing PRESTO's accelsearch on the FFT data (.fft)

    obsid - Observation ID of the object of interest (10-digit str)
    flags - a LIST of input flags for accelsearch
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    fft_files = sorted(glob.glob(nicersoft_output_folder+'*bary_*'+str(segment_length)+'s*.fft')) #not that order matters here I think, but just in case
    logfile = nicersoft_output_folder + 'accelsearch_segments.log'

    with open(logfile,'w') as logtextfile:
        for i in range(len(fft_files)):
            logtextfile.write(subprocess.check_output(['accelsearch']+flags+[fft_files[i]]))
        logtextfile.close()

    return

def prepfold(obsid,zmax):
    """
    Performing PRESTO's prepfold on the pulsation candidates.

    obsid - Observation ID of the object of interest (10-digit str)
    zmax - maximum acceleration
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'

    ACCEL_files = sorted(glob.glob(nicersoft_output_folder+'ni'+obsid+'_nicersoft_bary_*ACCEL_'+str(zmax)))
    cand_files = [ACCEL_files[i] + '.cand' for i in range(len(ACCEL_files))]
    events_files = [cand_files[i][:-15]+'.events' for i in range(len(cand_files))]
    logfile = nicersoft_output_folder + 'prepfold_segments.log'
    log = open(logfile,'a')

# Bad test! There's no .cand file if there're no candidates...
#    if len(cand_files) != len(events_files):
#        raise ValueError("Note that the number of .cand files is not equal to the number of .events files! Likely because accelsearch was not run on all the .events files, or look at Lv2_presto_segments - might be only partially completed.")

    header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
    header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

    for i in range(len(ACCEL_files)): #for each successful ACCEL_zmax file:
        accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
        index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
        index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
        no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates
        cand_relpath = relpath(cand_files[i],nicersoft_output_folder) #relative path of .cand file ; PRESTO doesn't like using absolute paths
        events_relpath = relpath(events_files[i],nicersoft_output_folder) #relative path of .events file ; PRESTO doesn't like using absolute paths

        for j in range(no_cands):
            command = 'cd ' + nicersoft_output_folder + ' ; prepfold -double -events -noxwin -n 50 -accelcand ' + str(j+1) + ' -accelfile ' + cand_relpath + ' ' + events_relpath
            subprocess.Popen(command,shell=True)

    log.close()

    return

def prepfold_segment(obsid,zmax,segment_length):
    """
    Performing PRESTO's prepfold on the pulsation candidates.

    obsid - Observation ID of the object of interest (10-digit str)
    zmax - maximum acceleration
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    #nicersoft_output_folder = './'
    ############### EDIT THE FUNCTION
    ACCEL_files = sorted(glob.glob(nicersoft_output_folder+'ni'+obsid+'_nicersoft_bary_*ACCEL_'+str(zmax)))
    cand_files = [ACCEL_files[i] + '.cand' for i in range(len(ACCEL_files))]
    events_files = [cand_files[i][:-15]+'.events' for i in range(len(cand_files))]
    logfile = nicersoft_output_folder + 'prepfold.log'
    log = open(logfile,'a')

# Bad test! There's no .cand file if there're no candidates...
#    if len(cand_files) != len(events_files):
#        raise ValueError("Note that the number of .cand files is not equal to the number of .events files! Likely because accelsearch was not run on all the .events files, or look at Lv2_presto_segments - might be only partially completed.")

    header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
    header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

    for i in range(len(ACCEL_files)): #for each successful ACCEL_zmax file:
        accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
        index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
        index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
        no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates
        cand_relpath = relpath(cand_files[i],nicersoft_output_folder) #relative path of .cand file ; PRESTO doesn't like using absolute paths
        events_relpath = relpath(events_files[i],nicersoft_output_folder) #relative path of .events file ; PRESTO doesn't like using absolute paths

        for j in range(no_cands):
            command = 'cd ' + nicersoft_output_folder + ' ; prepfold -double -events -noxwin -n 50 -accelcand ' + str(j+1) + ' -accelfile ' + cand_relpath + ' ' + events_relpath
            subprocess.Popen(command,shell=True)

    log.close()

    return

def ps2pdf(obsid):
    """
    Converting from .ps to .pdf

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    ps_files = glob.glob(nicersoft_output_folder + '*.ps') #grabs a list of files with extension .ps
    for i in range(len(ps_files)):
        pdf_file = ps_files[i][:-2]+'pdf' #replacing .ps to .pdf
        ps2pdf_command = ['ps2pdf',ps_files[i],pdf_file]
        subprocess.check_call(ps2pdf_command) #using ps2pdf to convert from ps to pdf ; not just a simple change in extension

    return

def ps2pdf_segment(obsid,segment_length):
    """
    Converting from .ps to .pdf

    obsid - Observation ID of the object of interest (10-digit str)
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(segment_length) + 's/'
    ps_files = glob.glob(nicersoft_output_folder + '*.ps') #grabs a list of files with extension .ps
    for i in range(len(ps_files)):
        pdf_file = ps_files[i][:-2]+'pdf' #replacing .ps to .pdf
        ps2pdf_command = ['ps2pdf',ps_files[i],pdf_file]
        subprocess.check_call(ps2pdf_command) #using ps2pdf to convert from ps to pdf ; not just a simple change in extension

    return

if __name__ == "__main__":
    get_gti_file('1060060127',100)
    niextract_gti_time('1060060127',100)
    niextract_gti_time_energy('0034070101',100,300,800)
    do_nicerfits2presto('1200250108',0.00025,200)
    edit_inf('1060060127',0.0005,1000)
    edit_binary('1060060127',0.0005,1000)
    realfft('1060060127')
    accelsearch_flags = ['-numharm','8','-zmax','200','-photon','-flo','1','-fhi','1000']
    accelsearch('1060060127',accelsearch_flags)
    prepfold('0034070101',10,0)
    ps2pdf('1200250126')
