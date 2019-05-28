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
import subprocess
import glob

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

#get_gti_file('1060060127',100)

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

#niextract_gti_time('1060060127',100)

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

    subprocess.check_call(['niextract-events',inputfile+'[pi='+str(PI1)+':'+str(PI2)+']',outputfile])

    return

#niextract_gti_energy('1060060127',200,300)

def do_nicerfits2presto(obsid,tbin,truncations):
    #use glob to find .evt!

def edit_binary(obsid):
    #will likely be the most intense function to write, remember that we want to find the duty cycle (maybe write separate functions
    #that I'll call in Lv2_presto_main), then do a plot of % vs segment length, then edit the binary
    #[MAKE SURE I TEST THIS OUT ON TEST.PY FIRST!!!]
    #use glob to find .dat
