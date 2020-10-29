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
import Lv0_dirs,Lv1_data_gtis,Lv2_mkdir
from scipy import stats
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
    times = fits.open(eventfile)[1].data['TIME']

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length) #array of time values, starting
    #Jul 10 2019: also added Tobs_end+segment_length instead of Tobs_end, to get that final piece of data

    binned_counts, bin_edges, binnumber = stats.binned_statistic(times,np.ones(len(times)),statistic='sum',bins=segment_times)
    #bin edges defined by left boundary

    gtilist = open(parent_folder + '/segment_GTI.list','w')

    gti_folder = parent_folder + '/gtis/'
    Lv2_mkdir.makedir(gti_folder)
    for i in tqdm(range(len(bin_edges)-1)):
        gtilist.write('GTI'+str(i).zfill(6)+ '      ' + str(bin_edges[i]) + '\n')
        gti_file = gti_folder + str(segment_length).zfill(5) + 's_GTI' + str(i).zfill(6) + '.gti'
        if binned_counts[i] != 0 and os.path.exists(gti_file)==False: #prevents from processing GTIs with no counts
            #print(str(segment_times[i]),str(segment_times[i+1]))
            subprocess.run(['mkgti.py','--gtiname',gti_file,str(segment_times[i]),str(segment_times[i+1])] )
            #no need to do 'startmet=', 'stopmet=', but make sure I do it in the right order!

    gtilist.close()

    return

def niextract_gti(eventfile,gap,gtifile):
    """
    Using niextract-events to get segmented data based on the individual GTIs created with
    GTI_bunching in Lv1_data_gtis.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    gap - maximum separation between end time of 1st GTI and start time of 2nd GTI allowed
    gtifile - name of GTI file
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    Lv1_data_gtis.GTI_bunching(eventfile,gap,gtifile)

    gtis = list(fits.open(parent_folder+'/'+gtifile)[1].data)
    niextract_folder = parent_folder + '/accelsearch_GTIs/'
    Lv2_mkdir.makedir(niextract_folder)
    for i in tqdm(range(len(gtis))):
        gti_start = gtis[i][0]
        gti_end = gtis[i][1]
        if os.path.exists(niextract_folder+'GTI'+str(i+1).zfill(6)+'.gti') == False:
            subprocess.run(['mkgti.py','--gtiname',niextract_folder+'GTI'+str(i+1).zfill(6)+'.gti',str(gti_start),str(gti_end)])

        niextract_output = niextract_folder+str(pathlib.Path(eventfile).name)[:-4]+'_GTI'+str(i+1).zfill(6)+'.evt'
        if os.path.exists(niextract_output)==False:
            subprocess.run(['niextract-events',eventfile,niextract_output,'timefile='+niextract_folder+'GTI'+str(i+1).zfill(6)+'.gti'])

    return

def niextract_gti_E(eventfile,gap,gtifile,PI1,PI2):
    """
    Using niextract-events to get segmented data based on the individual GTIs created with
    GTI_bunching in Lv1_data_gtis, AND with energy cuts

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    gap - maximum separation between end time of 1st GTI and start time of 2nd GTI allowed
    gtifile - name of GTI file
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    Lv1_data_gtis.GTI_bunching(eventfile,gap,gtifile)

    gtis = list(fits.open(parent_folder+'/'+gtifile)[1].data)
    niextract_folder = parent_folder + '/accelsearch_GTIs/'
    Lv2_mkdir.makedir(niextract_folder)
    for i in tqdm(range(len(gtis))):
        gti_start = gtis[i][0]
        gti_end = gtis[i][1]
        if os.path.exists(niextract_folder+'GTI'+str(i+1).zfill(6)+'.gti') == False:
            subprocess.run(['mkgti.py','--gtiname',niextract_folder+'GTI'+str(i+1).zfill(6)+'.gti',str(gti_start),str(gti_end)])

        niextract_output = niextract_folder+str(pathlib.Path(eventfile).name)[:-4]+'_GTI'+str(i+1).zfill(6)+'_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '.evt'
        if os.path.exists(niextract_output)==False:
            subprocess.run(['niextract-events',eventfile+'[PI='+str(PI1)+':'+str(PI2)+']',niextract_output,'timefile='+niextract_folder+'GTI'+str(i+1).zfill(6)+'.gti'])

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
    times = fits.open(eventfile)[1].data['TIME']

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

    binned_counts, bin_edges, binnumber = stats.binned_statistic(times,np.ones(len(times)),statistic='sum',bins=segment_times)
    #bin edges defined by left boundary

    niextract_folder = parent_folder + '/accelsearch_' + str(segment_length) + 's/'
    Lv2_mkdir.makedir(niextract_folder)
    for i in tqdm(range(len(segment_times)-1)):
        if binned_counts[i] != 0: #prevents from processing GTIs with no counts
            outputfile = niextract_folder + 'ni' + obsid + '_nicersoft_bary_GTI'+str(i).zfill(6)+'_'+str(segment_length).zfill(5)+'s.evt'
            if 'merged' in parent_folder:
                merged_id = str(pathlib.Path(eventfile).name)[:12]
                outputfile = niextract_folder + merged_id + '_nicersoft_bary_GTI' + str(i).zfill(6) + '_' + str(segment_length).zfill(5) + 's.evt'

            if os.path.exists(outputfile)==False:
                subprocess.run(['niextract-events',eventfile,outputfile,'timefile='+parent_folder+'/gtis/'+str(segment_length).zfill(5)+'s_GTI'+str(i).zfill(6)+'.gti'])

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
    output_file = niextract_folder + eventfile.split('/')[-1][:-4] + '_E'+str(PI1).zfill(4)+'-'+str(PI2).zfill(4)+'.evt'

    if os.path.exists(output_file)==False:
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
    times = fits.open(eventfile)[1].data['TIME']

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

#    Tobs_start = event[1].data['TIME'][0]
#    Tobs_end = event[1].data['TIME'][-1]

    segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length) #array of time values, starting
    #Jul 10: added Tobs_end + segment_length to get that final piece of data

    binned_counts, bin_edges, binnumber = stats.binned_statistic(times,np.ones(len(times)),statistic='sum',bins=segment_times)
    #bin edges defined by left boundary

    niextract_folder = parent_folder + '/accelsearch_' + str(segment_length) + 's/'
    Lv2_mkdir.makedir(niextract_folder)

    for i in tqdm(range(len(segment_times)-1)):
        if binned_counts[i] != 0: #prevents from processing GTIs with no counts
            outputfile = niextract_folder + 'ni' + obsid + '_nicersoft_bary_GTI'+str(i).zfill(6)+'_'+str(segment_length).zfill(5)+'s_' + 'E'+str(PI1).zfill(4)+'-'+str(PI2).zfill(4)+'.evt'
            if 'merged' in parent_folder:
                merged_id = str(pathlib.Path(eventfile).name)[:12]
                outputfile = niextract_folder + merged_id + '_nicersoft_bary_GTI' + str(i).zfill(6) + '_' + str(segment_length).zfill(5) + 's_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '.evt'
            if os.path.exists(outputfile)==False:
                subprocess.run(['niextract-events',eventfile+'[PI='+str(PI1)+':'+str(PI2)+']',outputfile,'timefile='+parent_folder+'/gtis/'+str(segment_length).zfill(5)+'s_GTI'+str(i).zfill(6)+'.gti'])

    return

def do_nicerfits2presto(eventfile,tbin,segment_length,mode):
    """
    Using nicerfits2presto.py to bin the data, and to convert into PRESTO-readable format.
    I can always move files to different folders to prevent repeats (especially for large files)
    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    tbin - size of the bins in time
    segment_length - length of the individual segments for combining power spectra
    mode - "all", "t", "gtis" or "E" ; basically to tell the function where to access files to run realfft for
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if mode == "all": #for non-segmented event file
        subprocess.run(['nicerfits2presto.py','--dt='+str(tbin),eventfile])

    if mode == "E":
        E_files = glob.glob(parent_folder+'/accelsearch_E/*.evt')
        if len(E_files) != 0: #if E_files is not empty
            for i in tqdm(range(len(E_files))):
                if os.path.exists(E_files[i][:-3] + 'dat'):
                    continue
                try:
                    subprocess.run(['nicerfits2presto.py','--dt='+str(tbin),E_files[i]])
                except (ValueError,subprocess.CalledProcessError):
                    pass

    if mode == 't':
        time_files = glob.glob(parent_folder+'/accelsearch_' + str(segment_length) + 's/*.evt')
        if len(time_files) != 0:
            for i in tqdm(range(len(time_files))):
                if os.path.exists(time_files[i][:-3] + 'dat'):
                    continue
                try:
                    subprocess.run(['nicerfits2presto.py','--dt='+str(tbin),time_files[i]])
                except (ValueError,subprocess.CalledProcessError):
                    pass

    if mode == "gtis":
        gti_files = glob.glob(parent_folder+'/accelsearch_GTIs/*.evt')
        if len(gti_files) != 0:
            for i in tqdm(range(len(gti_files))):
                if os.path.exists(gti_files[i][:-3] + 'dat'):
                    continue
                try:
                    subprocess.run(['nicerfits2presto.py','--dt='+str(tbin),gti_files[i]])
                except (ValueError,subprocess.CalledProcessError):
                    pass

    ##### will probably remove once I know to either NOT work in nicerpy_xrayanalysis,
    ##### and/or install my packages such that I can run the scripts from anywhere
    ##### in the terminal!

    presto_files = glob.glob('*'+obsid+'*')
    if 'merged' in eventfile:
        presto_files = glob.glob('merged*')
    for i in range(len(presto_files)):
        if mode == "all":
            subprocess.run(['mv',presto_files[i],parent_folder+'/'])
        if mode == "t":
            subprocess.run(['mv',presto_files[i],parent_folder+'/accelsearch_' + str(segment_length) + 's/'])
        if mode == "E":
            subprocess.run(['mv',presto_files[i],parent_folder+'/accelsearch_E/'])
        if mode == "gtis":
            subprocess.run(['mv',presto_files[i],parent_folder+'/accelsearch_GTIs/'])

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
    #obj_name = event_header['OBJECT']
    #obsid = event_header['OBS_ID']

    inf_files = sorted(glob.glob(parent_folder + '/accelsearch_' + str(segment_length) + 's/*GTI*' + str(segment_length).zfill(5)+'s*.inf')) #not the .evt file; some .evt files will be empty
    #inf_files = sorted(glob.glob(parent_folder + '/accelsearch_GTIs/*GTI*.inf')) #not the .evt file; some .evt files will be empty

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
    #obj_name = event_header['OBJECT']
    #obsid = event_header['OBS_ID']

    dat_files = sorted(glob.glob(parent_folder + '/accelsearch_' + str(segment_length) + 's/*GTI*' + str(segment_length).zfill(5) + 's*.dat')) #not that order matters here I think, but just in case
    #dat_files = sorted(glob.glob(parent_folder + '/accelsearch_GTIs/*GTI*.dat')) #not that order matters here I think, but just in case

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
    mode - "all", "t", "gtis", or "E" ; basically to tell the function where to access files to run realfft for
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
    elif mode == "gtis":
        dat_files = sorted(glob.glob(parent_folder+'/accelsearch_GTIs/*.dat'))
        logfile = parent_folder + '/accelsearch_GTIs/realfft_gtis.log'
    else:
        raise ValueError("mode should either of 'all', 't', 'gtis', or 'E'!")

    print("Running realfft now:")
    with open(logfile,'w') as logtextfile:
        for i in tqdm(range(len(dat_files))):
            if os.path.exists(dat_files[i][:-3] + 'fft')==False:
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
    mode - "all", "t", "gtis" or "E" ; basically to tell the function where to access files to run accelsearch for
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
    elif mode == "gtis":
        fft_files = sorted(glob.glob(parent_folder+'/accelsearch_GTIs/*.fft'))
        logfile = parent_folder + '/accelsearch_GTIs/accelsearch_gtis.log'
    else:
        raise ValueError("mode should either of 'all', 't', 'gtis', or 'E'!")

    zmax_index = flags.index('-zmax')
    zmax_num = str(flags[zmax_index+1])
    print("Running accelsearch now:")
    with open(logfile,'w') as logtextfile:
        for i in tqdm(range(len(fft_files))):
            #if os.path.exists(fft_files[i][:-4] + '_ACCEL_' + zmax_num + '.txtcand') == False:
            command = ['accelsearch'] + flags + [fft_files[i]]
            output = subprocess.run(command,capture_output=True,text=True)
            logtextfile.write(output.stdout)
            logtextfile.write('*------------------------------* \n')
            logtextfile.write(output.stderr)
        logtextfile.close()

    return

def prepfold(eventfile,segment_length,mode,zmax):
    """
    Performing PRESTO's prepfold on the pulsation candidates.
    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    segment_length - length of the individual segments
    mode - "all", "t", "gtis", or "E" ; basically to tell the function where to access files to run prepfold for
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
    elif mode == "gtis":
        ACCEL_files = sorted(glob.glob(parent_folder+'/accelsearch_GTIs/*ACCEL_'+str(zmax)))
        logfile = parent_folder + '/accelsearch_GTIs/prepfold.log'
    else:
        raise ValueError("mode should either of 'all', 't', 'gtis', or 'E'!")

    cand_files = [ACCEL_files[i] + '.cand' for i in range(len(ACCEL_files))]
    if zmax < 10:
        events_files = [cand_files[i][:-13]+'.events' for i in range(len(cand_files))]
    if (zmax >= 10) & (zmax < 100):
        events_files = [cand_files[i][:-14]+'.events' for i in range(len(cand_files))]
    if (zmax >= 100) & (zmax < 999):
        events_files = [cand_files[i][:-15]+'.events' for i in range(len(cand_files))]
    else:
        events_files = [cand_files[i][:-16]+'.events' for i in range(len(cand_files))]

    header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
    header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

    with open(logfile,'w') as logtextfile:
        for i in range(len(ACCEL_files)): #for each successful ACCEL_zmax file:
            accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
            index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
            index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
            no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates
            cand_relpath = relpath(cand_files[i],parent_folder) #relative path of .cand file ; PRESTO doesn't like using absolute paths
            events_relpath = relpath(events_files[i],parent_folder) #relative path of .events file ; PRESTO doesn't like using absolute paths

            for j in range(min(50,no_cands)): #to control the number of outputs! Sometimes there can be thousands of candidates - insane!
                command = 'cd ' + parent_folder + ' ; prepfold -double -events -noxwin -n 20 -accelcand ' + str(j+1) + ' -accelfile ' + cand_relpath + ' ' + events_relpath
                output = subprocess.run(command,shell=True,capture_output=True,text=True)
                logtextfile.write(output.stdout)
                logtextfile.write('*------------------------------* \n')
                logtextfile.write(output.stderr)
        logtextfile.close()

    return

def filter_accelsearch(eventfile,mode,min_freq,zmax):
    """
    Added on 8/31/2020. To filter out the ACCEL files by displaying a list of candidates
    that have met a frequency threshold (to avoid very low frequency candidates, really)

    The files will have been generated via prepfold already - perhaps in the future, I may
    only run prepfold on candidates I want! Maybe do this to avoid having too many candidates...

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    mode - "all", "t", "gtis", or "E" ; basically to tell the function where to access files to run prepfold for
    min_freq - minimum frequency to include in the list
    zmax - maximum acceleration
    """
    parent_folder = str(pathlib.Path(eventfile).parent)

    if mode == "all":
        ACCEL_files = sorted(glob.glob(parent_folder+'/*ACCEL_'+str(zmax)))
    elif mode == "t":
        ACCEL_files = sorted(glob.glob(parent_folder+'/accelsearch_'+str(segment_length)+'s/*ACCEL_'+str(zmax)))
    elif mode == "E":
        ACCEL_files = sorted(glob.glob(parent_folder+'/accelsearch_E/*ACCEL_'+str(zmax)))
    elif mode == "gtis":
        ACCEL_files = sorted(glob.glob(parent_folder+'/accelsearch_GTIs/*ACCEL_'+str(zmax)))
    else:
        raise ValueError("mode should either of 'all', 't', 'gtis', or 'E'!")

    ### below are just the headers in the *ACCEL_100 files ; the intention is just to get the indices for these lines later on
    header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
    header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

    output_file = open(str(pathlib.Path(ACCEL_files[0]).parent) + '/ACCEL_' + str(zmax) + '_cands_minfreq_' + str(min_freq) + '.txt','w')
    output_file.write(header1 + '\n')

    for i in range(len(ACCEL_files)):
        freqs = []
        accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
        index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
        index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
        no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates
        for j in range(no_cands):
            freq = float(accel_textfile[j+3].split()[6][:-3])
            freqs.append(freq)
        if not any(np.array(freqs)>min_freq): #if none are above the minimum frequency
            continue
        else:
            output_file.write(str(pathlib.Path(ACCEL_files[i]).name) + '\n')
            for j in range(no_cands):
                freq = float(accel_textfile[j+3].split()[6][:-3])
                if freq >= min_freq:
                    output_file.write(accel_textfile[j+3] + '\n')
    output_file.close()

def ps2pdf(eventfile,segment_length,mode):
    """
    Converting from .ps to .pdf
    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    mode - "all", "t", "gtis", or "E" ; basically to tell the function where to access files to run ps2pdf for
    """
    parent_folder = str(pathlib.Path(eventfile).parent)

    if mode == "all":
        ps_files = sorted(glob.glob(parent_folder+'/*ps'))
    elif mode == "t":
        ps_files = sorted(glob.glob(parent_folder+'/accelsearch_'+str(segment_length)+'s/*ps'))
    elif mode == "E":
        ps_files = sorted(glob.glob(parent_folder+'/accelsearch_E/*ps'))
    elif mode == "gtis":
        ps_files = sorted(glob.glob(parent_folder+'/accelsearch_GTIs/*ps'))
    else:
        raise ValueError("mode should either of 'all', 't', 'gtis', or 'E'!")

    print('Running ps2pdf now!')
    for i in tqdm(range(len(ps_files))):
        pdf_file = ps_files[i][:-2]+'pdf' #replacing .ps to .pdf
        ps2pdf_command = ['ps2pdf',ps_files[i],pdf_file]
        subprocess.run(ps2pdf_command) #using ps2pdf to convert from ps to pdf ; not just a simple change in extension

    return

if __name__ == "__main__":
    #eventfile = Lv0_dirs.NICERSOFT_DATADIR + '1034070101_pipe/ni1034070101_nicersoft_bary.evt'
    #eventfile = Lv0_dirs.NICER_DATADIR + 'xtej1739/2002131540_bary.evt'

    gap = 5
    gtifile = Lv0_dirs.NICER_DATADIR + 'xtej1739/bunched.gti'
    mode = 't'
    segment_length = 100
    PI1 = 100
    PI2 = 1000
    zmax = 100

    tbin = 0.025

    accelsearch_flags = ['-numharm','4','-zmax','100','-photon','-flo','1','-fhi','1000']

    #get_gti_file(eventfile,segment_length)
    #niextract_gti(eventfile,gap,gtifile)
    #niextract_gti_time(eventfile,segment_length)
    #niextract_gti_energy(eventfile,PI1,PI2)
    #niextract_gti_time_energy(eventfile,segment_length,PI1,PI2)
    #do_nicerfits2presto(eventfile,tbin,segment_length)
    #edit_inf(eventfile,tbin,segment_length)
    #edit_binary(eventfile,tbin,segment_length)
    #realfft(eventfile,segment_length,mode)
    #accelsearch(eventfile,segment_length,mode,accelsearch_flags)
    #prepfold(eventfile,mode,zmax)
    #ps2pdf(eventfile,mode)
    #testfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034070101_pipe/niextract/ni1034070101_nicersoft_bary_GTI0_100s.dat'
    #bins = np.fromfile(testfile,dtype='<f',count=-1)
    #print(len(bins))

    #####
    #eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/at2018cow_2020/at2018cow_filtered_highSNR.evt'
    #gap = 0
    #gtifile = 'bunched.gti'
    #niextract_gti(eventfile,gap,gtifile)

    #####
    eventfile = '/Volumes/Samsung_T5/NICER-data/at2019wey/2008041401_filt_bary.evt'
    mode = 'gtis'
    min_freq = 10
    zmax = 100

    filter_accelsearch(eventfile,mode,min_freq,zmax)
