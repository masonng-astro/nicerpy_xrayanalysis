#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs May 30 4:03pm 2019

Program which plots the histogram of number of powers above some certain value, vs
the power...

Should just need .fft itself? Can incorporate with my own program in the future too...

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv2_mkdir
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import glob
import subprocess

Lv0_dirs.global_par()

def histogram(obsid,tbin,segment_length):
    """
    Plot a histogram of N(>P) vs P, directly from the .fft file! They should not be
    particularly large files, so for now (5/30), let's just 'glob' in all the .fft files,
    and save them to PDFs!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsdir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    fft_files = sorted(glob.glob(obsdir+'*.fft')) #grab the .fft files

    hist_dir = obsdir + 'powers_hist/' #create the directory where we'll store the PDFs of the histogram plots
    Lv2_mkdir.makedir(hist_dir)

    for i in tqdm(range(len(fft_files))):
        binned_data = np.fromfile(fft_files[i],dtype='<f',count=-1)
#        print(fft_files[i],len(binned_data),min(binned_data[int(len(binned_data)/2):-1]),max(binned_data[int(len(binned_data)/2):-1]))
#        binned_data = binned_data[binned_data>0]
#        print(len(binned_data))
        f,(ax1,ax2) = plt.subplots(2,1)
        freqs = np.linspace(-1/(tbin*2),1/(tbin*2),len(binned_data)+1)
        ax1.hist(np.log(np.abs(binned_data[freqs[:-1]>1])),bins=50,log=True) #log number of bars!
        ax2.semilogy(freqs[:-1],np.abs(binned_data))
        ax2.set_xlim([1,2000])
        plt.show()
    #    print(np.mean(np.abs(binned_data[2500000:])))
        #plt.savefig(fft_files[i][:-3]+'pdf',dpi=900,format='pdf')
        #plt.close()

    fft_pdfs = glob.glob(obsdir+'*'+str(segment_length)+'s.pdf')
    for i in range(len(fft_pdfs)):
        subprocess.check_call(['mv',fft_pdfs[i],obsdir+'powers_hist/'])

    return

#histogram('1034090111',0.00025,1000)
