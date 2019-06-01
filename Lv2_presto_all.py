#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 2:22pm 2019

Program for doing realfft, accelsearch, prepfold, and ps2pdf.
This is for when we're using data from the WHOLE time series though.
Use Lv2_presto_segments.py if you want to look at segments instead!

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs
import os
import subprocess
import glob

Lv0_dirs.global_par()

def realfft(obsid):
    """
    Performing PRESTO's realfft on the binned data (.dat)

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    dat_file = nicersoft_output_folder + 'ni' + obsid + '_nicersoft_bary.dat'
    logfile = nicersoft_output_folder + 'realfft.log'

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['realfft',dat_file]))
        logtextfile.close()

    return

def accelsearch(obsid,flags):
    """
    Performing PRESTO's accelsearch on the FFT data (.fft)

    obsid - Observation ID of the object of interest (10-digit str)
    flags - a LIST of input flags for psrpipe
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    fft_file = nicersoft_output_folder + 'ni' + obsid + '_nicersoft_bary.fft'
    logfile = nicersoft_output_folder + 'accelsearch.log'

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['accelsearch']+flags+[fft_file]))
        logtextfile.close()

    return

#accelsearch_flags = ['-numharm','8','-zmax','200','-photon','-flo','1','-fhi','1000']

def prepfold(obsid,no_cand,zmax):
    """
    Performing PRESTO's prepfold on the pulsation candidates.

    obsid - Observation ID of the object of interest (10-digit str)
    no_cand - number of candidates. I haven't yet thought of a way to automate this,
    so I'll have to do all of the above first BEFORE I do prepfold. It's fine though,
    since this is the only 'big' manual step.
    zmax - maximum acceleration
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    cand_file = nicersoft_output_folder + 'ni' + obsid + '_nicersoft_bary_ACCEL_' + str(zmax) + '.cand'
    events_file = nicersoft_output_folder + 'ni' + obsid + '_nicersoft_bary.events'
    logfile = nicersoft_output_folder + 'accelsearch.log'

    if os.path.isfile(logfile): #basically to overwrite the previous accelsearch.log
        os.remove(logfile)

    with open(logfile,'a+') as logtextfile:
        for i in range(1,no_cand+1):
            subprocess.check_call(['prepfold','-double','-events','-noxwin','-accelcand',str(i),'-accelfile',cand_file,events_file])
            logtextfile.write(subprocess.check_output(['accelsearch']+flags+[fft_file]))
            logtextfile.close()

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
        pdf_file = ps_files[i].replace('.ps','.pdf') #replacing .ps to .pdf
        subprocess.check_call(['ps2pdf',ps_files[i],pdf_file]) #using ps2pdf to convert from ps to pdf ; not just a simple change in extension

    return
