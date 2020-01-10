#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import glob
import subprocess
import os
from os.path import relpath
import Lv0_dirs

Lv0_dirs.global_par()

##### INITIAL LOOK AT THE HISTOGRAM

obsids = ['12002501' + str(i+1).zfill(2) for i in range(26)]
cands_file = open('candidates.txt','w')
for k in range(len(obsids)):
    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR
    ACCEL_files = sorted(glob.glob(nicersoft_output_folder+obsids[k]+'_pipe/accelsearch_64s/ni'+obsids[k]+'_nicersoft_bary_*ACCEL_100')) #getting absolute paths for the *ACCEL_100 files
    header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
    header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

    for i in range(len(ACCEL_files)): #for each successful ACCEL_zmax file:
        accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
        index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
        index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
        no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates

        candidates = np.genfromtxt(ACCEL_files[i],dtype='str',skip_header=3,usecols=6,unpack=True,max_rows=no_cands) #to obtain frequency values
        for j in range(candidates.size):
            if candidates.size == 1: #if there's only one candidate...
                cands_file.write(ACCEL_files[i] + ' ' + str(candidates)[:-3] + '\n') #[:-3] takes away the uncertainty
            else:
                cands_file.write(ACCEL_files[i] + ' ' + candidates[j][:-3] + '\n')

cands_file.close()

candidates_file = 'candidates.txt'
check_cands = np.genfromtxt(candidates_file,dtype='float',skip_footer=0) #get the frequency values
filtered_check_cands = check_cands[check_cands>1]
print(len(check_cands))
print(len(filtered_check_cands))

plt.figure(1)
plt.hist(filtered_check_cands,bins=100)

##### NEXT, DOING THE SHUFFLING OF THE LIGHT CURVES!
"""
obsids = ['12002501' + str(i+1).zfill(2) for i in range(26)]
at2018cow_noise = '/Volumes/Samsung_T5/NICERsoft_outputs/at2018cow_noise/'

for i in range(len(obsids)):
    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR
    accelsearch_folder = nicersoft_output_folder + obsids[i] + '_pipe/accelsearch_64s/'
    ACCEL_files = sorted(glob.glob(nicersoft_output_folder+obsids[i]+'_pipe/accelsearch_64s/ni'+obsids[i]+'_nicersoft_bary_*ACCEL_100')) #getting absolute paths for the *ACCEL_100 files
    dat_files = [ACCEL_files[j][:-10] + '.dat' for j in range(len(ACCEL_files))] #get absolute paths of binary .dat files
    inf_files = [ACCEL_files[j][:-10] + '.inf' for j in range(len(ACCEL_files))] #get absolute paths of inf files (information files for PRESTO)

    for k in range(len(dat_files)): #for every .dat file
        counts = np.fromfile(dat_files[k],dtype='<f',count=-1) #read the binary data into Python
        np.random.shuffle(counts) #shuffle the counts in each 1ms bin
        filename = relpath(dat_files[k],accelsearch_folder) #get the relative path of the binary data file
        counts.tofile(at2018cow_noise + filename) #save the NEW file into the folder 'at2018cow_noise'

    for k in range(len(inf_files)): #to save the inf files into the same 'at2018cow_noise' folder as the data files. Need this for accelsearch
        inf_relpath = relpath(inf_files[k],accelsearch_folder)
        inf_newpath = at2018cow_noise + inf_relpath
        subprocess.check_output(['mv',inf_files[k],inf_newpath])

all_dat_files = sorted(glob.glob(at2018cow_noise+'*.dat'))
for i in range(len(all_dat_files)):
    subprocess.check_output(['realfft',all_dat_files[i]]) #need to do FFT on the binary data file

all_fft_files = sorted(glob.glob(at2018cow_noise+'*.fft'))
logfile = at2018cow_noise + 'accelsearch.log'
accelsearch_flags = ['-numharm','8','-zmax','100','-photon','-flo','1','-fhi','500'] #accelsearch flags
with open(logfile,'w') as logtextfile:
    for i in range(len(all_fft_files)):
        logtextfile.write(subprocess.check_output(['accelsearch']+accelsearch_flags+[all_fft_files[i]]))
    logtextfile.close()
"""
##### GET THE HISTOGRAMS AGAIN
at2018cow_noise = '/Volumes/Samsung_T5/NICERsoft_outputs/at2018cow_noise/'

cands_file = open('candidates_noise.txt','w')
ACCEL_files = sorted(glob.glob(at2018cow_noise +'*ACCEL_100'))
header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

for i in range(len(ACCEL_files)): #for each successful ACCEL_zmax file:
    accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
    index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
    index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
    no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates

    candidates = np.genfromtxt(ACCEL_files[i],dtype='str',skip_header=3,usecols=6,unpack=True,max_rows=no_cands)
    for j in range(candidates.size):
        if candidates.size == 1:
            cands_file.write(ACCEL_files[i] + ' ' + str(candidates)[:-3] + '\n')
        else:
            cands_file.write(ACCEL_files[i] + ' ' + candidates[j][:-3] + '\n')

cands_file.close()

noise_candidates_file = 'candidates_noise.txt'
noise_check_cands = np.genfromtxt(noise_candidates_file,dtype='float',skip_footer=0,usecols=1,unpack=True)
noise_filtered_check_cands = noise_check_cands[noise_check_cands>1]
print(len(noise_check_cands))
print(len(noise_filtered_check_cands))

#plt.figure(2)
plt.hist(noise_filtered_check_cands,bins=100)

plt.legend(('Original search','Shuffled data'),loc='best')

plt.axvline(x=260,lw=0.5,alpha=0.5)
plt.axvline(x=200,lw=0.5,alpha=0.5)
plt.axvline(x=100,lw=0.5,alpha=0.5)
plt.axvline(x=130,lw=0.5,alpha=0.5)
plt.axvline(x=50,lw=0.5,alpha=0.5)
plt.axvline(x=65,lw=0.5,alpha=0.5)

plt.xlabel('Frequency (Hz)',fontsize=12)
plt.ylabel('Number',fontsize=12)

plt.show()
