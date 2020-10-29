#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Created on Wed Aug 5 5:54pm 2020

Script that visualizes the candidates from PRESTO's acceleration search algorithm.

"""
from __future__ import division, print_function
import mplcursors
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import numpy as np

def visual_accelsearch(timetype,min_freq,min_sigma):

    ACCEL_files = glob.glob('/Volumes/Samsung_T5/NICER-data/at2019wey/accelsearch_GTIs/*E0200-1200_ACCEL_100') #gets the path to all the *ACCEL_100 files
    eventfiles_ACCEL = [str(ACCEL_files[i][:-10]) + '.evt' for i in range(len(ACCEL_files))] #gets the path to all the corresponding event files where there exists a *ACCEL_100 file
    gti_no = [float(ACCEL_files[i][-25:-21]) for i in range(len(ACCEL_files))] #gets the GTI number where there exists a *ACCEL_100 file

    ### below are just the headers in the *ACCEL_100 files ; the intention is just to get the indices for these lines later on
    header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
    header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

    for i in range(len(ACCEL_files)): #for each successful ACCEL_zmax file
        freqs = []
        sigmas = []
        accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
        index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
        index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
        no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates
        for j in range(no_cands):
            freq = accel_textfile[j+3].split()[6][:-3]
            sigma = float(accel_textfile[j+3].split()[1])
            if float(freq) > min_freq and float(sigma) > min_sigma:
                freqs.append(float(freq))
                sigmas.append(float(sigma))
        if timetype == 'seconds':
            times = np.ones(len(freqs))*(fits.open(eventfiles_ACCEL[i])[1].header['TSTART']+fits.open(eventfiles_ACCEL[i])[1].header['TSTOP'])/2
            plt.xlabel('Time (s)')
        elif timetype == 'GTIs':
            times = np.ones(len(freqs))*gti_no[i]
            plt.xlabel('GTI number',fontsize=12)
        #print(times,freqs,sigmas)
        plt.scatter(x=times,y=freqs,c=sigmas,marker='o',cmap='gist_heat',vmin=3,vmax=4,edgecolors='k')

    #plt.yscale('log')
    plt.colorbar().set_label('Significance (sigma)')
    plt.title('Energy range: 2-12 keV',fontsize=12)

    mplcursors.cursor(hover=True)
    plt.ylabel('Frequency',fontsize=12)
    plt.show()

if __name__ == "__main__":
    #visual_accelsearch('seconds',10)
    visual_accelsearch('GTIs',10,3)
