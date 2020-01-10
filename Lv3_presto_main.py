#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 1:55pm 2019

Program which incorporates all of Lv2_preprocess and Lv2_presto_subroutines

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv0_scp,Lv0_nicerl2
import Lv2_preprocess,Lv2_presto_subroutines,Lv2_mkdir
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
conv_fits2presto = False
process_all = False #whether to process the ENTIRE observation or not

process_segments = True #whether to make truncations to the data
time_segments = True #truncate by time segments
energy_segments = False #truncate by energy range
time_energy_segments = False #truncate both by time segments and energy range

accelsearch = True

##### From Lv2_preprocess
obsids = ['20600603' + str(i) for i in range(51,53)]
parfile = ''
obsdirs = [Lv0_dirs.NICER_DATADIR + obsids[i] + '/' for i in range(len(obsids))]
orbitfiles = [Lv0_dirs.NICER_DATADIR + obsids[i] + '/auxil/ni' + obsids[i] + '.orb' for i in range(len(obsids))]
nicer_datafiles = [Lv0_dirs.NICER_DATADIR + obsids[i] + '/xti/event_cl/ni' + obsids[i] + '_0mpu7_cl.evt' for i in range(len(obsids))]
nicer_baryoutputs = [Lv0_dirs.NICER_DATADIR + obsids[i] + '/xti/event_cl/ni' + obsids[i] + '_0mpu7_cl_bary.evt' for i in range(len(obsids))]
nicersoft_datafiles = [Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/cleanfilt.evt' for i in range(len(obsids))]
nicersoft_baryoutputs = [Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/ni' + obsids[i] + '_nicersoft_bary.evt' for i in range(len(obsids))]
nicersoft_folder = [Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/' for i in range(len(obsids))]

nicerl2_flags = ['clobber=YES']
psrpipe_flags = ['--emin','0.3','--emax','12.0','--mask','14','34','54'] #for psrpipe in Lv0_psrpipe
refframe = 'ICRS' #for barycorr in Lv1_barycorr
tbin = '0.00025' #time bin for PRESTO in seconds

##### From Lv2_presto_all
##### For when we analyze the entire data set (no truncations in energy or time)
accelsearch_flags = ['-numharm','8','-zmax','100','-photon','-flo','1','-fhi','500']

##### Parameters for prepfold
prepfold = False
zmax = 100

##### From Lv2_presto_subroutines
segment_lengths = [1000] #desired length of segments (and 200)
## NOTE: nicerfits2presto tries to find a 'nice number' of bins to bin the data, and
## sometimes that will be longer than the segment. I'm not sure what the best remedy is,
## but currently, the best thing is just to track the segment size on the terminal, then increase
## accordingly. For example, I tried 100s segments for 0034070101 (Cen X-3), but nicerfits2presto.py
## is binning at 100.8s, so I went up to the next nice number - 102..
PI1 = [30]
PI2 = [200]

##### For Lv3_duty_cycle:
duty_cycle_bin = 1 #for 1s bins to do duty cycle calculations
threshold = 20 #for 10%

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
        Lv2_preprocess.preprocess(obsdirs[i],nicerl2_flags,psrpipe_flags,refframe,orbitfiles[i],parfile,nicer_datafiles[i],nicer_baryoutputs[i],nicersoft_datafiles[i],nicersoft_baryoutputs[i],nicersoft_folder[i])

################################################################################
################################## PRESTO_ALL ##################################
################################################################################

if process_all == True:
    if conv_fits2presto == True:
        for i in range(len(obsids)):
            Lv2_presto_subroutines.do_nicerfits2presto(nicersoft_baryoutputs[i],tbin,0)

    if accelsearch == True:
        print('Doing realfft/accelsearch now!')
        for i in range(len(obsids)):
            ### Running realfft and accelsearch from PRESTO
            Lv2_presto_subroutines.realfft(nicersoft_baryoutputs[i],0,'all')
            Lv2_presto_subroutines.accelsearch(nicersoft_baryoutputs[i],0,'all',accelsearch_flags)

    ## no_cand might be a list, if I'm processing multiple ObsIDs at once...
    if prepfold == True:
        for i in range(len(obsids)):
            Lv2_presto_all.prepfold(obsids[i],zmax)
            Lv2_presto_all.ps2pdf(obsids[i])

################################################################################
############################### PRESTO_SEGMENTS ################################
################################################################################

if process_segments == True:

    if time_segments == True:
        for i in range(len(obsids)):
            for j in range(len(segment_lengths)):
                Lv2_presto_subroutines.get_gti_file(nicersoft_baryoutputs[i],segment_lengths[j]) #make GTI files for each segment
                Lv2_presto_subroutines.niextract_gti_time(nicersoft_baryoutputs[i],segment_lengths[j]) #performing niextract-events

                Lv2_presto_subroutines.do_nicerfits2presto(nicersoft_baryoutputs[i],tbin,segment_lengths[j])
                Lv2_presto_subroutines.edit_inf(nicersoft_baryoutputs[i],tbin,segment_lengths[j])
                Lv2_presto_subroutines.edit_binary(nicersoft_baryoutputs[i],tbin,segment_lengths[j])
#
#                Lv3_duty_cycle.duty_cycle(obsids[i],tbin,segment_lengths[j],duty_cycle_bin,threshold)
#                Lv3_duty_cycle.duty_cycle_dist(obsids[i],tbin,segment_lengths[j],duty_cycle_bin,threshold)

    if energy_segments == True:
        if len(PI1) != len(PI2):
            raise ValueError("Make sure that the length of PI1 and PI2 are the same! Need pairs of PI values.")
        for i in range(len(obsids)):
            for j in range(len(PI1)):
                Lv2_presto_subroutines.niextract_gti_energy(nicersoft_baryoutputs[i],PI1[j],PI2[j])
                Lv2_presto_subroutines.do_nicerfits2presto(nicersoft_baryoutputs[i],tbin,1) #segment_length makes no sense, so 1 is a placeholder

    if time_energy_segments == True:
        for i in range(len(obsids)):
            for j in range(len(segment_lengths)):
                Lv2_presto_subroutines.get_gti_file(nicersoft_baryoutputs[i],segment_lengths[j]) #make GTI files for each segment
                for k in range(len(PI1)):
                    Lv2_presto_subroutines.niextract_gti_time_energy(nicersoft_baryoutputs[i],segment_lengths[j],PI1[k],PI2[k])

                Lv2_presto_subroutines.do_nicerfits2presto(nicersoft_baryoutputs[i],tbin,segment_lengths[j])
                Lv2_presto_subroutines.edit_inf(nicersoft_baryoutputs[i],tbin,segment_lengths[j])
                Lv2_presto_subroutines.edit_binary(nicersoft_baryoutputs[i],tbin,segment_lengths[j])

                #for k in range(len(PI1)):
                #    Lv3_duty_cycle.duty_cycle_tE(obsids[i],tbin,segment_lengths[j],PI1[k],PI2[k],duty_cycle_bin,threshold)
                #    Lv3_duty_cycle.duty_cycle_tE_dist(obsids[i],tbin,segment_lengths[j],PI1[k],PI2[k],duty_cycle_bin,threshold)

    ### Running realfft and accelsearch from PRESTO
    if accelsearch == True:
        if time_segments == True or time_energy_segments == True:
            for i in range(len(obsids)):
                for j in range(len(segment_lengths)):
                    Lv2_presto_subroutines.realfft(nicersoft_baryoutputs[i],segment_lengths[j],'t')
                    Lv2_presto_subroutines.accelsearch(nicersoft_baryoutputs[i],segment_lengths[j],'t',accelsearch_flags)
        if energy_segments == True:
            for i in range(len(obsids)):
                Lv2_presto_subroutines.realfft(nicersoft_baryoutputs[i],0,'E')
                Lv2_presto_subroutines.accelsearch(nicersoft_baryoutputs[i],0,'E',accelsearch_flags)
        else:
            "None of time_segments, time_energy_segments, or energy_segments are True!"

    ## no_cand might be a list, if I'm processing multiple ObsIDs at once...
    if prepfold == True:
        if time_segments == True or time_energy_segments == True:
            for i in range(len(obsids)):
                Lv2_presto_subroutines.prepfold(nicersoft_baryoutputs[i],'t',zmax)
        if energy_segments == True:
            for i in range(len(obsids)):
                Lv2_presto_subroutines.prepfold(nicersoft_baryoutputs[i],'E',zmax)

        ### doing ps2pdf
        if time_segments == True or time_energy_segments == True:
            for i in range(len(obsids)):
                Lv2_presto_subroutines.ps2pdf(nicersoft_baryoutputs[i],'t')
        if energy_segments == True:
            for i in range(len(obsids)):
                Lv2_presto_subroutines.ps2pdf(nicersoft_baryoutputs[i],'E')

        else:
            "None of time_segments, time_energy_segments, or energy_segments are True!"
