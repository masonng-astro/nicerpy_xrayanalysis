#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28th 4:57pm 2019

Temporary script to obtain spectra, from background-subtracted spectra, that have
the response and arf matrices folded in.

The script will generate the .xcm file needed to be read into XSPEC!
Will also break up the big all_spectra.txt file into individual text files! 

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess
import glob
import Lv0_dirs
from matplotlib.backends.backend_pdf import PdfPages

Lv0_dirs.global_par()

"""
response = '/Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf'
arf = '/Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf'

spectra = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/bgsub_cl50*.pha'))

data_list = []
response_list = []
arf_list = []

for i in range(len(spectra)):
    data_list.append(str(i+1)+':'+str(i+1))
    data_list.append(spectra[i])

    response_list.append(str(i+1))
    response_list.append(response)

    arf_list.append(str(i+1))
    arf_list.append(arf)

data_string = 'data ' + ' '.join([data_list[i] for i in range(len(data_list))])
response_string = 'response ' + ' '.join(response_list[i] for i in range(len(response_list)))
arf_string = 'arf ' + ' '.join(arf_list[i] for i in range(len(arf_list)))
ignore_string = 'ignore **:0.0-0.285,12.01-**'

xcm_file = open('get_spectra.xcm','w')
xcm_file.write(data_string + '\n')
xcm_file.write(response_string + '\n')
xcm_file.write(arf_string + '\n')
xcm_file.write(ignore_string + '\n')
xcm_file.close()
"""

##### Get spectra, and save them into individual files. Remember that I got them
##### simultaneously in XSPEC!

spectra_bgsub = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/bgsub_cl50*.pha'))
spectra_all = '/Volumes/Samsung_T5/NGC300_ULX/all_spectra.txt'
spectra = np.array(open(spectra_all,'r').read().split('\n'))

separator = list(np.where(spectra=='NO NO NO NO')[0])
separator.insert(0,2)
separator.insert(len(separator),len(spectra)-1)

for i in tqdm(range(len(spectra_bgsub))):
    if i != 397:
        spectrum_file = spectra_bgsub[i][:-3] + 'txt'
        spectrum = open(spectrum_file,'w')

        if i < 397:
            for j in range(separator[i]+1,separator[i+1]):
                spectrum.write(spectra[j] + '\n')
        else:
            for j in range(separator[i-1]+1,separator[i]):
                spectrum.write(spectra[j] + '\n')

        spectrum.close()


#data 1:1 bgsub_cl50_100.pha 2:2 bgsub_cl50_101.pha 3:3 bgsub_cl50_102.pha
#response 1 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf 2 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf 3 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf
#arf 1 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf 2 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf 3 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf
#ignore **:0.0-0.285,12.01-**
