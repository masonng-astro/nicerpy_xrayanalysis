#!/usr/bin/env python
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
import Lv0_dirs,Lv2_mkdir
from tqdm import tqdm
import os
from os.path import relpath
import time
import pathlib
import subprocess
import glob

Lv0_dirs.global_par()

powers_of_two = [2**i for i in range(31)] # for deciding number of FFT bins

def get_gti_file(eventfile,segment_length):
    """
    Creating the individual .gti files for my data segments!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    segment_length - length of the individual segments for combining power spectra
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    gtis = fits.open(eventfile)[2].data

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length) #array of time values, starting
    #Jul 10 2019: also added Tobs_end+segment_length instead of Tobs_end, to get that final piece of data

    gti_folder = parent_folder + '/gtis/'
    Lv2_mkdir.makedir(gti_folder)
    for i in tqdm(range(len(segment_times)-1)):
        gti_file = gti_folder + str(segment_length) + 's_GTI' + str(i) + '.gti'
        #print(str(segment_times[i]),str(segment_times[i+1]))
        subprocess.run(['mkgti.py','--gtiname',gti_file,str(segment_times[i]),str(segment_times[i+1])] )
        #no need to do 'startmet=', 'stopmet=', but make sure I do it in the right order!

    return

def niextract_gti_time(eventfile,segment_length):
    """
    Using niextract-events to get segmented data based on the [segment_length]-length
    GTIs that were created above!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    segment_length - length of the individual segments for combining power spectra
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    gtis = fits.open(eventfile)[2].data

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time
#    Tobs_start = event[1].data['TIME'][0]
#    Tobs_end = event[1].data['TIME'][-1]

    segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length) #array of time values, starting
    #from the observation start time until the observation end time, in steps of segment length.
    #Also had Tobs_end+segment_length to get the last section of data. It's ok if it's short, will zero-pad in
    #edit_binary! Done Jul 10.

    niextract_folder = parent_folder + '/accelsearch_' + str(segment_length) + 's/'
    Lv2_mkdir.makedir(niextract_folder)
    for i in tqdm(range(len(segment_times)-1)):
        outputfile = niextract_folder + 'ni' + obsid + '_nicersoft_bary_GTI'+str(i)+'_'+str(segment_length)+'s.evt'
        subprocess.run(['niextract-events',eventfile,outputfile,'timefile='+parent_folder+'/gtis/'+str(segment_length)+'s_GTI'+str(i)+'.gti'])

    return

def niextract_gti_energy(eventfile,PI1,PI2):
    """
    Using niextract-events to get segmented data based on the energy range

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    gtis = fits.open(eventfile)[2].data

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    niextract_folder = parent_folder + '/accelsearch_E/'
    Lv2_mkdir.makedir(niextract_folder)
    output_file = niextract_folder + eventfile.split('/')[-1][:-4] + '_E'+str(PI1)+'-'+str(PI2)+'.evt'

    subprocess.run(['niextract-events',eventfile+'[PI='+str(PI1)+':'+str(PI2)+']',output_file])

    return

def niextract_gti_time_energy(eventfile,segment_length,PI1,PI2):
    """
    Using niextract-events to get segmented data based on [segment_length]-length
    GTIs that were created above, AND energy range!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    segment_length - length of the individual segments for combining power spectra
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    gtis = fits.open(eventfile)[2].data

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

#    Tobs_start = event[1].data['TIME'][0]
#    Tobs_end = event[1].data['TIME'][-1]

    segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length) #array of time values, starting
    #Jul 10: added Tobs_end + segment_length to ge tthat final piece of data

    niextract_folder = parent_folder + '/accelsearch_' + str(segment_length) + 's/'
    Lv2_mkdir.makedir(niextract_folder)

    for i in tqdm(range(len(segment_times)-1)):
        output_file = niextract_folder + 'ni' + obsid + '_nicersoft_bary_GTI'+str(i)+'_'+str(segment_length)+'s_' + 'E'+str(PI1)+'-'+str(PI2)+'.evt'
        subprocess.run(['niextract-events',eventfile+'[PI='+str(PI1)+':'+str(PI2)+']',output_file,'timefile='+parent_folder+'/gtis/'+str(segment_length)+'s_GTI'+str(i)+'.gti'])

    return

def do_nicerfits2presto(eventfile,tbin,segment_length):
    """
    Using nicerfits2presto.py to bin the data, and to convert into PRESTO-readable format.
    I can always move files to different folders to prevent repeats (especially for large files)

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    tbin - size of the bins in time
    segment_length - length of the individual segments for combining power spectra
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if segment_length == 0: #for non-segmented event file
        subprocess.run(['nicerfits2presto.py','--dt='+str(tbin),eventfile])

    else:
        E_files = glob.glob(parent_folder+'/accelsearch_E/*.evt')
        time_files = glob.glob(parent_folder+'/accelsearch_' + str(segment_length) + 's/*.evt')

        if len(E_files) != 0: #if E_files is not empty
            for i in tqdm(range(len(E_files))):
                try:
                    subprocess.run(['nicerfits2presto.py','--dt='+str(tbin),E_files[i]])
                except (ValueError,subprocess.CalledProcessError):
                    pass
        if len(time_files) != 0:
            for i in tqdm(range(len(time_files))):
                try:
                    subprocess.run(['nicerfits2presto.py','--dt='+str(tbin),time_files[i]])
                except (ValueError,subprocess.CalledProcessError):
                    pass

    ##### will probably remove once I know to either NOT work in nicerpy_xrayanalysis,
    ##### and/or install my packages such that I can run the scripts from anywhere
    ##### in the terminal!

    presto_files = glob.glob('*'+obsid+'*')
    for i in range(len(presto_files)):
        if segment_length == 0:
            subprocess.run(['mv',presto_files[i],parent_folder+'/'])
        elif 'GTI' in presto_files[i]:
            subprocess.run(['mv',presto_files[i],parent_folder+'/accelsearch_' + str(segment_length) + 's/'])
        else:
            subprocess.run(['mv',presto_files[i],parent_folder+'/accelsearch_E/'])

    return

def edit_inf(eventfile,tbin,segment_length):
    """
    Editing the .inf file, as it seems like accelsearch uses some information from the .inf file!
    Mainly need to edit the "Number of bins in the time series".
    This is only for when we make segments by time though!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    inf_files = sorted(glob.glob(parent_folder + '/accelsearch_' + str(segment_length) + 's/*GTI*' + str(segment_length)+'s*.inf')) #not the .evt file; some .evt files will be empty

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

def edit_binary(eventfile,tbin,segment_length):
    """
    To pad the binary file so that it will be as long as the desired segment length.
    The value to pad with for each time bin, is the average count rate in THAT segment!
    Jul 10: Do zero-padding instead... so that number of counts is consistent!
    Again, this is only for when we make segments by time!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    dat_files = sorted(glob.glob(parent_folder + '/accelsearch_' + str(segment_length) + 's/*GTI*' + str(segment_length) + 's*.dat')) #not that order matters here I think, but just in case
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

    return

def realfft(eventfile,segment_length,mode):
    """
    Performing PRESTO's realfft on the binned data (.dat)

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    segment_length - length of the individual segments
    mode - "all", "t" or "E" ; basically to tell the function where to access files to run realfft for
    """
    parent_folder = str(pathlib.Path(eventfile).parent)

    if mode == "all":
        dat_files = sorted(glob.glob(parent_folder+'/*.dat'))
        logfile = parent_folder + '/realfft_all.log'
    elif mode == "t":
        dat_files = sorted(glob.glob(parent_folder+'/accelsearch_'+str(segment_length)+'s/*.dat'))
        logfile = parent_folder + '/accelsearch_' + str(segment_length) + 's/realfft_t.log'
    elif mode == "E":
        dat_files = sorted(glob.glob(parent_folder+'/accelsearch_E/*.dat'))
        logfile = parent_folder + '/accelsearch_E/realfft_E.log'
    else:
        raise ValueError("mode should either of 'all', 't', or 'E'!")

    print("Running realfft now:")
    with open(logfile,'w') as logtextfile:
        for i in tqdm(range(len(dat_files))):
            output = subprocess.run(['realfft',dat_files[i]],capture_output=True,text=True)
            logtextfile.write(output.stdout)
            logtextfile.write('*------------------------------* \n')
            logtextfile.write(output.stderr)
        logtextfile.close()

    return

def accelsearch(eventfile,segment_length,mode,flags):
    """
    Performing PRESTO's accelsearch on the FFT data (.fft)

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    segment_length - length of the individual segments
    mode - "all", "t" or "E" ; basically to tell the function where to access files to run accelsearch for
    flags - a LIST of input flags for accelsearch
    """
    if type(flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    parent_folder = str(pathlib.Path(eventfile).parent)

    if mode == "all":
        fft_files = sorted(glob.glob(parent_folder+'/*.fft'))
        logfile = parent_folder + '/accelsearch_all.log'
    elif mode == "t":
        fft_files = sorted(glob.glob(parent_folder+'/accelsearch_'+str(segment_length)+'s/*.fft'))
        logfile = parent_folder + '/accelsearch_' + str(segment_length) + 's/accelsearch_t.log'
    elif mode == "E":
        fft_files = sorted(glob.glob(parent_folder+'/accelsearch_E/*.fft'))
        logfile = parent_folder + '/accelsearch_E/accelsearch_E.log'
    else:
        raise ValueError("mode should either of 'all', 't', or 'E'!")

    print("Running accelsearch now:")
    with open(logfile,'w') as logtextfile:
        for i in tqdm(range(len(fft_files))):
            command = ['accelsearch'] + flags + [fft_files[i]]
            output = subprocess.run(command,capture_output=True,text=True)
            logtextfile.write(output.stdout)
            logtextfile.write('*------------------------------* \n')
            logtextfile.write(output.stderr)
        logtextfile.close()

    return

def prepfold(eventfile,mode,zmax):
    """
    Performing PRESTO's prepfold on the pulsation candidates.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    mode - "all", "t" or "E" ; basically to tell the function where to access files to run prepfold for
    zmax - maximum acceleration
    """
    parent_folder = str(pathlib.Path(eventfile).parent)

    if mode == "all":
        ACCEL_files = sorted(glob.glob(parent_folder+'/*ACCEL_'+str(zmax)))
        logfile = parent_folder + '/prepfold_all.log'
    elif mode == "t":
        ACCEL_files = sorted(glob.glob(parent_folder+'/accelsearch_'+str(segment_length)+'s/*ACCEL_'+str(zmax)))
        logfile = parent_folder + '/accelsearch_' + str(segment_length) + 's/prepfold_t.log'
    elif mode == "E":
        ACCEL_files = sorted(glob.glob(parent_folder+'/accelsearch_E/*ACCEL_'+str(zmax)))
        logfile = parent_folder + '/accelsearch_E/prepfold.log'
    else:
        raise ValueError("mode should either of 'all', 't', or 'E'!")

    cand_files = [ACCEL_files[i] + '.cand' for i in range(len(ACCEL_files))]
    events_files = [cand_files[i][:-15]+'.events' for i in range(len(cand_files))]

    header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
    header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

    for i in range(len(ACCEL_files)): #for each successful ACCEL_zmax file:
        accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
        index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
        index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
        no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates
        cand_relpath = relpath(cand_files[i],parent_folder) #relative path of .cand file ; PRESTO doesn't like using absolute paths
        events_relpath = relpath(events_files[i],parent_folder) #relative path of .events file ; PRESTO doesn't like using absolute paths

        for j in range(no_cands):
            command = 'cd ' + nicersoft_output_folder + ' ; prepfold -double -events -noxwin -n 50 -accelcand ' + str(j+1) + ' -accelfile ' + cand_relpath + ' ' + events_relpath
            subprocess.Popen(command,shell=True)

    return

def ps2pdf(eventfile,mode):
    """
    Converting from .ps to .pdf

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    mode - "all", "t" or "E" ; basically to tell the function where to access files to run ps2pdf for
    """
    parent_folder = str(pathlib.Path(eventfile).parent)

    if mode == "all":
        ps_files = sorted(glob.glob(parent_folder+'/*ps'))
    elif mode == "t":
        ps_files = sorted(glob.glob(parent_folder+'/accelsearch_'+str(segment_length)+'s/*ps'))
    elif mode == "E":
        ps_files = sorted(glob.glob(parent_folder+'/accelsearch_E/*ps'))
    else:
        raise ValueError("mode should either of 'all', 't', or 'E'!")

    for i in range(len(ps_files)):
        pdf_file = ps_files[i][:-2]+'pdf' #replacing .ps to .pdf
        ps2pdf_command = ['ps2pdf',ps_files[i],pdf_file]
        subprocess.run(ps2pdf_command) #using ps2pdf to convert from ps to pdf ; not just a simple change in extension

    return

if __name__ == "__main__":
    eventfile = Lv0_dirs.NICERSOFT_DATADIR + '1034070101_pipe/ni1034070101_nicersoft_bary.evt'
    mode = 't'
    segment_length = 100
    PI1 = 100
    PI2 = 1000
    zmax=100

    tbin = 0.025

    accelsearch_flags = ['-numharm','4','-zmax','100','-photon','-flo','1','-fhi','1000']

    get_gti_file(eventfile,segment_length)
    niextract_gti_time(eventfile,segment_length)
    niextract_gti_energy(eventfile,PI1,PI2)
    niextract_gti_time_energy(eventfile,segment_length,PI1,PI2)
    do_nicerfits2presto(eventfile,tbin,segment_length)
    edit_inf(eventfile,tbin,segment_length)
    edit_binary(eventfile,tbin,segment_length)
    realfft(eventfile,segment_length,mode)
    accelsearch(eventfile,segment_length,mode,accelsearch_flags)
    prepfold(eventfile,mode,zmax)
    ps2pdf(eventfile,mode)
    #testfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034070101_pipe/niextract/ni1034070101_nicersoft_bary_GTI0_100s.dat'
    #bins = np.fromfile(testfile,dtype='<f',count=-1)
    #print(len(bins))
