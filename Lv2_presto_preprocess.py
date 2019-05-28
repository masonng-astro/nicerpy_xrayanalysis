#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 2:22pm 2019

Program for doing Lv0_scp, Lv0_gunzip, Lv0_psrpipe, Lv1_barycorr, that is,
this is pre-processing data before running them through to PRESTO!

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv0_scp,Lv0_gunzip,Lv0_psrpipe,Lv1_barycorr
import os
import subprocess
import glob

Lv0_dirs.global_par()

def preprocess(obsid,psrpipe_flags,refframe,tbin):
    """
    Preprocessing the NICER data for use in PRESTO, so running scp, gunzip, psrpipe,
    barycorr, and nicerfits2presto.

    obsid - Observation ID of the object of interest (10-digit str)
    psrpipe_flags - a LIST of input flags for psrpipe
    refframe - reference frame for barycenter corrections (usually ICRS)
    tbin - time bin for PRESTO data in seconds
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(psrpipe_flags) != list:
        raise TypeError("flags should be a list! Not even an array.")
    if refframe != 'ICRS' and refframe != 'FK5':
        raise ValueError("refframe should either be ICRS or FK5! Otherwise, update Lv1_barycorr.py if there are options I was unaware of.")

    logfile = obsid+'_preprocess.log' # name of the log file
    nicersoft_dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/' #absolute path for the $OBSID_pipe folder

    Lv0_scp.scp(obsid) #secure copying the decrypted folder from ciri
    Lv0_gunzip.unzip_all(obsid) #unzipping the contents within the observation
    Lv0_psrpipe.psrpipe(obsid,psrpipe_flags) #applying custom cuts (though no need --shrinkelv after HEASOFT 6.26)
    Lv1_barycorr.barycorr(obsid,refframe) #applying barycenter corrections to the cut data
    subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),nicersoft_dir+'ni'+obsid+'_nicersoft_bary.evt']) #converting to the PRESTO data format

    subprocess.check_call(['mv',obsid+'_psrpipe.log',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/']) #copying the psrpipe log file to the $OBSID_pipe folder
    subprocess.check_call(['mv','ni'+obsid+'_nicersoft_bary.dat',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/']) #copying the psrpipe log file to the $OBSID_pipe folder
    subprocess.check_call(['mv','ni'+obsid+'_nicersoft_bary.events',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/']) #copying the psrpipe log file to the $OBSID_pipe folder
    subprocess.check_call(['mv','ni'+obsid+'_nicersoft_bary.inf',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/']) #copying the psrpipe log file to the $OBSID_pipe folder

obsid = '1060060127'
psrpipe_flags = ['--emin','0.3','--emax','12.0','--shrinkelvcut'] #for psrpipe in Lv0_psrpipe
refframe = 'ICRS' #for barycorr in Lv1_barycorr
tbin = '0.00025' #time bin for PRESTO in seconds

#preprocess(obsid,psrpipe_flags,refframe,tbin)
