#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 2.32pm 2019.

Performing NICER-specific HEASOFT functions. Specifically niextract-events for now.

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs,Lv0_call_eventcl,Lv0_call_ufa,Lv1_data_gtis,Lv1_data_gtis_filter
import subprocess
from astropy.io import fits
from tqdm import tqdm

Lv0_dirs.global_par()

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
        output_file = output_dir + 'GTI' + str(i) + '.gti'
        subprocess.check_call(['mkgti.py','--gtiname',output_file,str(segment_times[i]),str(segment_times[i+1])] )
        #no need to do 'startmet=', 'stopmet=', but make sure I do it in the right order!

    return

#get_gti_file('1034090108',1000)

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

        subprocess.check_call(['niextract-events',inputfile,outputfile,'timefile='+working_dir+'GTI'+str(i)+'.gti'])

    return

#niextract_gti('1034090108',1000)

"""
Updated on 5/23/19. This is old! The above incorporated advice from Paul Ray!

def niextract(obsid,bary,gap,desired_length,segment_length):

    Using the start and end time of the observation, break up the observation
    into segments of [segment_length] seconds each! Can do 1000s for example.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    gap - gap between data that would still be considered data
    -- for example, if you have data for t=1000-1040 and t=1045-1100, then if
    gap < 5, you get TWO GTIs. If gap >= 5, then you have 1 GTI from 1000-1100.
    This is different to desired_length - this is just meant to 'join up data'
    which are separated by mere seconds!
    desired_length - if GTIs > some desired length, use that data. This is used to
    weed out short GTIs.
    segment_length - length of the individual segments for combining power spectra

    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")

    shifted_gtis_desired,actual_gtis_desired = Lv1_data_gtis_filter.desired_gtis(obsid,bary,gap,desired_length)

    Tobs_start = actual_gtis_desired[0][0] #MJD for the observation start time
    Tobs_end = actual_gtis_desired[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end,segment_length) #array of time values, starting
    #from the observation start time until the observation end time, in steps of 1000s.
    #This means that you'll get, for a 4392s observation, np.array([0,1000,2000,3000,4000])!
    #We'd lose 392s worth of counts, but it can't be used in combining the power spectra anyways.

    for i in tqdm(range(len(segment_times)-1)):
        inputfile = Lv0_dirs.NICERSOFT_DATADIR+obsid+"_pipe/ni" + obsid + '_nicersoft_bary.evt'
        filter = '[TIME=' + str(segment_times[i])+':'+str(segment_times[i+1])+']'
        outputfile = Lv0_dirs.NICERSOFT_DATADIR+obsid+"_pipe/ni" + obsid + '_nicersoft_bary_GTI'+str(i+1)+'_1000s.evt'

        subprocess.check_call(['niextract-events',inputfile+filter,outputfile])

    return

#niextract('1034090111',True,50,1000,1000)

def niextract_stitch(obsid,bary,segment_length):

    Using the start and end time of the observation, break up the observation
    into segments of [segment_length] seconds each! Can do 1000s for example.
    This is done with the STITCHED data file -

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    gap - gap between data that would still be considered data
    -- for example, if you have data for t=1000-1040 and t=1045-1100, then if
    gap < 5, you get TWO GTIs. If gap >= 5, then you have 1 GTI from 1000-1100.
    This is different to desired_length - this is just meant to 'join up data'
    which are separated by mere seconds!
    desired_length - if GTIs > some desired length, use that data. This is used to
    weed out short GTIs.
    segment_length - length of the individual segments for combining power spectra

    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary_stitched.evt'
    event = fits.open(event)
    gtis = event[2].data

    Tobs_start = gtis[0][0] #MJD for the observation start time
    Tobs_end = gtis[-1][1] #MJD for the observation end time

    segment_times = np.arange(Tobs_start,Tobs_end,segment_length) #array of time values, starting
    #from the observation start time until the observation end time, in steps of 1000s.
    #This means that you'll get, for a 4392s observation, np.array([0,1000,2000,3000,4000])!
    #We'd lose 392s worth of counts, but it can't be used in combining the power spectra anyways.

    for i in tqdm(range(len(segment_times)-1)):
        inputfile = Lv0_dirs.NICERSOFT_DATADIR+obsid+"_pipe/ni" + obsid + '_nicersoft_bary_stitched.evt'
        filter = '[TIME=' + str(segment_times[i])+':'+str(segment_times[i+1])+']'
        outputfile = Lv0_dirs.NICERSOFT_DATADIR+obsid+"_pipe/ni" + obsid + '_nicersoft_bary_GTI'+str(i+1)+'_1000s_stitched.evt'

        subprocess.check_call(['niextract-events',inputfile+filter,outputfile])

    return

#niextract_stitch('1034090111',True,1000)
"""
