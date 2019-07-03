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
from os.path import relpath
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
    logfile = nicersoft_output_folder + 'accelsearch_all.log'

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['accelsearch']+flags+[fft_file]))
        logtextfile.close()

    return

#accelsearch_flags = ['-numharm','8','-zmax','200','-photon','-flo','1','-fhi','1000']

def prepfold(obsid,zmax):
    """
    Performing PRESTO's prepfold on the pulsation candidates.

    obsid - Observation ID of the object of interest (10-digit str)
    zmax - maximum acceleration
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    ACCEL_files = sorted(glob.glob(nicersoft_output_folder+'ni'+obsid+'_nicersoft_bary_ACCEL_'+str(zmax)))
    cand_files = [ACCEL_files[i] + '.cand' for i in range(len(ACCEL_files))]
    events_files = [cand_files[i][:-15]+'.events' for i in range(len(cand_files))]
    logfile = nicersoft_output_folder + 'prepfold_all.log'
    log = open(logfile,'a')

    if os.path.isfile(logfile): #basically to overwrite the previous accelsearch.log
        os.remove(logfile)

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
        pdf_file = ps_files[i].replace('.ps','.pdf') #replacing .ps to .pdf
        subprocess.check_call(['ps2pdf',ps_files[i],pdf_file]) #using ps2pdf to convert from ps to pdf ; not just a simple change in extension

    return

if __name__ == "__main__":
    print('hi') #placeholder
