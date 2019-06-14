#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 1:55pm 2019

Program which incorporates all of Lv2_presto_preprocess, Lv2_presto_segments, and
Lv2_presto_all.

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv0_scp,Lv0_nicerl2,Lv2_presto_preprocess,Lv2_presto_segments,Lv2_presto_all,Lv2_mkdir
import Lv3_duty_cycle
from tqdm import tqdm
import time
import os
import subprocess
import glob

################################################################################
########################## REMEMBER TO RUN HEAINIT !!! #########################
################################################################################

Lv0_dirs.global_par()
nicersoft_dir = Lv0_dirs.NICERSOFT_DATADIR

preprocess = False #I only need to do preprocessing for each ObsID ONCE!!!
process_all = False #whether to process the ENTIRE observation or not

process_segments = True #whether to make truncations to the data
time_segments = False #truncate by time segments
energy_segments = True #truncate by energy range
time_energy_segments = False #truncate both by time segments and energy range

accelsearch = False  #really only for segments

##### From Lv2_presto_preprocess
#obsids = ['1034090111']
#obsids = ['12002501' + str(i+1).zfill(2) for i in range(26)]
#obsids = ['12002501' + str(i).zfill(2) for i in [1,2,3,6,7,8,9,14,15,16,17,18,20,21,24,26]] # Jun 5 prepfold for AT2018cow
obsids = ['12002501' + str(i).zfill(2) for i in [6,7,14,15,16,18,20,24,26]] #accelsearch for AT2018cow pulsations in short observation durations
#obsids = ['1200250108']
nicerl2_flags = ['clobber=YES']
psrpipe_flags = ['--emin','0.3','--emax','12.0'] #for psrpipe in Lv0_psrpipe
refframe = 'ICRS' #for barycorr in Lv1_barycorr
tbin = '0.00025' #time bin for PRESTO in seconds

##### From Lv2_presto_all
##### For when we analyze the entire data set (no truncations in energy or time)
accelsearch_flags = ['-numharm','8','-zmax','100','-photon','-flo','1','-fhi','2000']

##### Parameters for prepfold
prepfold = True
zmax = 100

##### From Lv2_presto_segments
segment_lengths = [200,500,1000] #desired length of segments (and 200)
## NOTE: nicerfits2presto tries to find a 'nice number' of bins to bin the data, and
## sometimes that will be longer than the segment. I'm not sure what the best remedy is,
## but currently, the best thing is just to track the segment size on the terminal, then increase
## accordingly. For example, I tried 100s segments for 0034070101 (Cen X-3), but nicerfits2presto.py
## is binning at 100.8s, so I went up to the next nice number - 102..
PI1 = [30,250,500,800]
PI2 = [250,500,800,1200]

##### For Lv3_duty_cycle:
duty_cycle_bin = 1 #for 1s bins to do duty cycle calculations
threshold = 10 #for 10%

################################################################################
################################ PRE-PROCESSING ################################
################################################################################

### This just runs Lv2_presto_preprocess.preprocess, but within it, it's actually
### running scp, gunzip, psrpipe, barycorr, and nicerfits2presto. Note that this
### is for the entire data set, if we want to truncate the data

## also remember to run heainit at first!

if preprocess == True:
    print('Preprocessing now!')
    for i in range(len(obsids)):
        ### doing it here rather than in Lv2_presto_preprocess so that I just have
        ### to type the password in once
        Lv0_scp.scp(obsids[i]) #secure copying the decrypted folder from ciri

    for i in range(len(obsids)):
        Lv2_presto_preprocess.preprocess(obsids[i],nicerl2_flags,psrpipe_flags,refframe,tbin)

################################################################################
################################## PRESTO_ALL ##################################
################################################################################

if process_all == True:

    if preprocess == True:
        for i in range(len(obsids)):
            subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),nicersoft_dir+obsids[i]+'_pipe/ni'+obsids[i]+'_nicersoft_bary.evt']) #converting to the PRESTO data format
            subprocess.check_call(['mv','ni'+obsids[i]+'_nicersoft_bary.dat',Lv0_dirs.NICERSOFT_DATADIR+obsids[i]+'_pipe/'])
            subprocess.check_call(['mv','ni'+obsids[i]+'_nicersoft_bary.events',Lv0_dirs.NICERSOFT_DATADIR+obsids[i]+'_pipe/'])
            subprocess.check_call(['mv','ni'+obsids[i]+'_nicersoft_bary.inf',Lv0_dirs.NICERSOFT_DATADIR+obsids[i]+'_pipe/'])

    print('Doing realfft/accelsearch now!')
    for i in range(len(obsids)):
        ### Running realfft and accelsearch from PRESTO
        Lv2_presto_all.realfft(obsids[i])
        Lv2_presto_all.accelsearch(obsids[i],accelsearch_flags)

        ## no_cand might be a list, if I'm processing multiple ObsIDs at once...
        if prepfold == True:
            Lv2_presto_all.prepfold(obsids[i],no_cand,zmax)
            Lv2_presto_all.ps2pdf(obsids[i])

################################################################################
############################### PRESTO_SEGMENTS ################################
################################################################################

if process_segments == True:

    if time_segments == True:
        for i in range(len(obsids)):
            for j in range(len(segment_lengths)):
                Lv2_presto_segments.get_gti_file(obsids[i],segment_lengths[j]) #make GTI files for each segment
                Lv2_presto_segments.niextract_gti_time(obsids[i],segment_lengths[j]) #performing niextract-events

                Lv2_presto_segments.do_nicerfits2presto(obsids[i],tbin,segment_lengths[j])
                Lv2_presto_segments.edit_inf(obsids[i],tbin,segment_lengths[j])
                Lv2_presto_segments.edit_binary(obsids[i],tbin,segment_lengths[j])

                evt_files = glob.glob(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/*_'+str(segment_lengths[j])+'*.evt')
                dat_files = glob.glob(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/*_'+str(segment_lengths[j])+'*.dat')
                inf_files = glob.glob(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/*_'+str(segment_lengths[j])+'*.inf')
                events_files = glob.glob(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/*_'+str(segment_lengths[j])+'*.events')
                accelsearch_dir = Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/accelsearch_' + str(segment_lengths[j]) + 's/'
                Lv2_mkdir.makedir(accelsearch_dir)
                for k in range(len(evt_files)):
                    subprocess.check_call(['mv',evt_files[k],accelsearch_dir])
                for k in range(len(dat_files)):
                    subprocess.check_call(['mv',dat_files[k],accelsearch_dir])
                    subprocess.check_call(['mv',inf_files[k],accelsearch_dir])
                    subprocess.check_call(['mv',events_files[k],accelsearch_dir])

                Lv3_duty_cycle.duty_cycle(obsids[i],tbin,segment_lengths[j],duty_cycle_bin,threshold)
                Lv3_duty_cycle.duty_cycle_dist(obsids[i],tbin,segment_lengths[j],duty_cycle_bin,threshold)

#    if energy_segments == True:
#        if len(PI1) != len(PI2):
#            raise ValueError("Make sure that the length of PI1 and PI2 are the same! Need pairs of PI values.")
#        for i in range(len(obsids)):
#            for j in range(len(PI1)):
#                Lv2_presto_segments.niextract_gti_energy(obsids[i],PI1[j],PI2[j])
#                Lv2_presto_segments.do_nicerfits2presto(obsids[i],tbin,1) #segment_length makes no sense, so 1 is a placeholder

    if time_energy_segments == True:
        for i in range(len(obsids)):
            for j in range(len(segment_lengths)):
                Lv2_presto_segments.get_gti_file(obsids[i],segment_lengths[j]) #make GTI files for each segment
                for k in range(len(PI1)):
                    Lv2_presto_segments.niextract_gti_time_energy(obsids[i],segment_lengths[j],PI1[k],PI2[k])

                Lv2_presto_segments.do_nicerfits2presto(obsids[i],tbin,segment_lengths[j])
                Lv2_presto_segments.edit_inf(obsids[i],tbin,segment_lengths[j])
                Lv2_presto_segments.edit_binary(obsids[i],tbin,segment_lengths[j])

                evt_files = glob.glob(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/*_'+str(segment_lengths[j])+'*.evt')
                dat_files = glob.glob(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/*_'+str(segment_lengths[j])+'*.dat')
                inf_files = glob.glob(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/*_'+str(segment_lengths[j])+'*.inf')
                events_files = glob.glob(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/*_'+str(segment_lengths[j])+'*.events')

                accelsearch_dir = Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/accelsearch_' + str(segment_lengths[j]) + 's/'
                Lv2_mkdir.makedir(accelsearch_dir)
                for k in range(len(evt_files)):
                    subprocess.check_call(['mv',evt_files[k],accelsearch_dir])
                for k in range(len(dat_files)):
                    subprocess.check_call(['mv',dat_files[k],accelsearch_dir])
                    subprocess.check_call(['mv',inf_files[k],accelsearch_dir])
                    subprocess.check_call(['mv',events_files[k],accelsearch_dir])

                for k in range(len(PI1)):
                    Lv3_duty_cycle.duty_cycle_tE(obsids[i],tbin,segment_lengths[j],PI1[k],PI2[k],duty_cycle_bin,threshold)
                    Lv3_duty_cycle.duty_cycle_tE_dist(obsids[i],tbin,segment_lengths[j],PI1[k],PI2[k],duty_cycle_bin,threshold)

    ### Running realfft and accelsearch from PRESTO
    if accelsearch == True:
        if time_segments == True or time_energy_segments == True:
            for i in range(len(obsids)):
                for j in range(len(segment_lengths)):
                    Lv2_presto_segments.realfft_segment(obsids[i],segment_lengths[j])
                    Lv2_presto_segments.accelsearch_segment(obsids[i],accelsearch_flags,segment_lengths[j])
        elif energy_segments == True:
            for i in range(len(obsids)):
                Lv2_presto_segments.realfft(obsids[i])
                Lv2_presto_segments.accelsearch(obsids[i],accelsearch_flags)
        else:
            "None of time_segments, time_energy_segments, or energy_segments are True!"

    ## no_cand might be a list, if I'm processing multiple ObsIDs at once...
    if prepfold == True:
        if time_segments == True or time_energy_segments == True:
            for i in range(len(obsids)):
                for j in range(len(segment_lengths)):
                    Lv2_presto_segments.prepfold_segment(obsids[i],zmax,segment_lengths[j])
        elif energy_segments == True:
            for i in range(len(obsids)):
                Lv2_presto_segments.prepfold(obsids[i],zmax)
        else:
            "None of time_segments, time_energy_segments, or energy_segments are True!"

        ### doing ps2pdf
        if time_segments == True or time_energy_segments == True:
            for i in range(len(obsids)):
                for j in range(len(segment_lengths)):
                    Lv2_presto_segments.ps2pdf_segment(obsids[i],segment_lengths[j])
        elif energy_segments == True:
            for i in range(len(obids)):
                Lv2_presto_segments.ps2pdf(obsids[i])
        else:
            "None of the time_segments, time_energy_segments, or energy_segments are True!"