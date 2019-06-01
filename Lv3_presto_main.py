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
import Lv0_dirs,Lv0_scp,Lv0_nicerl2,Lv2_presto_preprocess,Lv2_presto_segments,Lv2_presto_all
from tqdm import tqdm
import os
import subprocess
import glob

################################################################################
########################## REMEMBER TO RUN HEAINIT !!! #########################
################################################################################

Lv0_dirs.global_par()

preprocess = False #I only need to do preprocessing for each ObsID ONCE!!!
process_all = False #whether to process the ENTIRE observation or not

process_segments = True #whether to make truncations to the data
time_segments = False #truncate by time segments
energy_segments = False #truncate by energy range
time_energy_segments = False #truncate both by time segments and energy range

accelsearch = False  #really only for segments

##### From Lv2_presto_preprocess
obsids = ['1034090111']
#obsids = ['12002501' + str(i+1).zfill(2) for i in range(26)]
nicerl2_flags = ['elv=20','br_earth=30','clobber=YES']
psrpipe_flags = ['--emin','0.3','--emax','12.0','--shrinkelvcut'] #for psrpipe in Lv0_psrpipe
refframe = 'ICRS' #for barycorr in Lv1_barycorr
tbin = '0.00025' #time bin for PRESTO in seconds

##### From Lv2_presto_all
##### For when we analyze the entire data set (no truncations in energy or time)
accelsearch_flags = ['-numharm','8','-zmax','200','-photon','-flo','1','-fhi','2000']

##### Parameters for prepfold
prepfold = True
no_cands = [2,1,10,1,3,3,4,7,1,3,6,2,1,3,3,16,5,16,2,2,4,19,5,2,10]
zmax = 200

##### From Lv2_presto_segments
segment_length = 1000 #desired length of segments
## NOTE: nicerfits2presto tries to find a 'nice number' of bins to bin the data, and
## sometimes that will be longer than the segment. I'm not sure what the best remedy is,
## but currently, the best thing is just to track the segment size on the terminal, then increase
## accordingly. For example, I tried 100s segments for 0034070101 (Cen X-3), but nicerfits2presto.py
## is binning at 100.8s, so I went up to the next nice number - 102..
PI1 = [300]
PI2 = [800]

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
        subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),nicersoft_dir+'ni'+obsids[i]+'_nicersoft_bary.evt']) #converting to the PRESTO data format

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
            Lv2_presto_segments.get_gti_file(obsids[i],segment_length) #make GTI files for each segment
            Lv2_presto_segments.niextract_gti_time(obsids[i],segment_length) #performing niextract-events

            Lv2_presto_segments.do_nicerfits2presto(obsids[i],tbin)
            Lv2_presto_segments.edit_inf(obsids[i],tbin,segment_length)
            Lv2_presto_segments.edit_binary(obsids[i],tbin,segment_length)

    if energy_segments == True:
        if len(PI1) != len(PI2):
            raise ValueError("Make sure that the length of PI1 and PI2 are the same! Need pairs of PI values.")
        for i in range(len(obsids)):
            for j in range(len(PI1)):
                Lv2_presto_segments.niextract_gti_energy(obsids[i],PI1[j],PI2[j])

    if time_energy_segments == True:
        for i in range(len(obsids)):
            Lv2_presto_segments.get_gti_file(obsids[i],segment_length) #make GTI files for each segment
            for j in range(len(PI1)):
                Lv2_presto_segments.niextract_gti_time_energy(obsids[i],segment_length,PI1[j],PI2[j])

            Lv2_presto_segments.do_nicerfits2presto(obsids[i],tbin)
            Lv2_presto_segments.edit_inf(obsids[i],tbin,segment_length)
            Lv2_presto_segments.edit_binary(obsids[i],tbin,segment_length)

    ### Running realfft and accelsearch from PRESTO
    for i in range(len(obsids)):
        if accelsearch == True:
            Lv2_presto_segments.realfft(obsids[i])
            Lv2_presto_segments.accelsearch(obsids[i],accelsearch_flags)

    ## no_cand might be a list, if I'm processing multiple ObsIDs at once...
        if prepfold == True:
            Lv2_presto_segments.prepfold(obsids[i],no_cands,zmax)
            Lv2_presto_segments.ps2pdf(obsids[i])
