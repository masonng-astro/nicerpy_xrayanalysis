#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 7th 2:28pm 2019

Binning the NGC300 data in a custom way

"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from astropy.io import fits
import subprocess
import pathlib
import glob
import Lv0_dirs

Lv0_dirs.global_par()

binsize_s = 5
bin_size = str(binsize_s).zfill(2) + 'd' #bins of 1 day!
bgsub_type = 'xsbgsub'

#RGcms = Lv0_dirs.NGC300_2020 + 'n300_ulx.bgsub_cl50_g2020.fffphot' NOT USING for 2020. No fffphot file for this anyways
norm_fffphot = Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020norm.fffphot'

### extracting photometry/count rates for soft_1 (0.2-0.3 keV), soft_2 (0.3-0.4 keV),
### A (0.4-1 keV), B (1-2 keV), C (2-4 keV), D (4-12) keV, in-band (A+B+C+D = 0.4-12 keV) bands, plus MJD
time = np.genfromtxt(norm_fffphot,dtype='float64',usecols=(11),unpack=True)
filenames = np.genfromtxt(norm_fffphot,dtype='str',usecols=(0),unpack=True)

time_bins = np.arange(58239,58606,binsize_s)
time = np.floor(time) #can't use int ; basically round the time values down to the integer value

def binned_text(bin_size):
    """
    Given the MJDs, binned counts, and associated uncertainties, put them into a text file

    No arguments because I'll put all the bands in here
    """
    E_bins_low = np.array([20-1,30-1,40-1,100-1,200-1,400-1,40-1,1300-1])
    E_bins_high = np.array([30-1,40-1,100-1,200-1,400-1,1200-1,1200-1,1501-1])

    bgsub_files = sorted(glob.glob(Lv0_dirs.NGC300_2020 + 'spectra_' + bin_size + '/58*_' + bgsub_type + '*_cl50.pha'))

    output_text = open(Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020norm_' + str(binsize_s).zfill(2) + 'd.fffphot','w')
    output_text_err = open(Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020err_norm_' + str(binsize_s).zfill(2) + 'd.fffphot','w')

    for i in tqdm(range(len(bgsub_files))): #for each averaged spectrum
        mjd = str(pathlib.Path(bgsub_files[i]).name)[:5]
        files_in_interval = filenames[(time>=float(mjd))&(time<float(mjd)+binsize_s)]

        if bgsub_type == 'bgsub':
            filtered_files = [Lv0_dirs.NGC300_2020 + 'bgsub_cl50/pha/' + files_in_interval[a] for a in range(len(files_in_interval))]
        elif bgsub_type == 'xbgsub':
            filtered_files = [Lv0_dirs.NGC300_2020 + 'bgsub_cl50/xpha/' + files_in_interval[a] for a in range(len(files_in_interval))]
        elif bgsub_type == 'xsbgsub':
            filtered_files = [Lv0_dirs.NGC300_2020 + 'bgsub_cl50/xspha/' + files_in_interval[a] for a in range(len(files_in_interval))]

        rates_all = []
        errs_all = []
        for j in range(len(E_bins_low)): #for each energy band
            rates_all.append(sum(fits.open(bgsub_files[i])[1].data['RATE'][E_bins_low[j]:E_bins_high[j]]))
            errs = fits.open(bgsub_files[i])[1].data['STAT_ERR'][E_bins_low[j]:E_bins_high[j]]
            errs_all.append( np.sqrt(np.mean(errs**2) * len(errs)))

        rate_1 = mjd #for the MJD
        rate_2 = '' #for the count rates in each band
        rate_3 = '' #for the list of files that goes into each bin
        for j in range(len(rates_all)):
            rate_2 += str(round(rates_all[j],4)) + ' '
        for j in range(len(filtered_files)):
            if j != len(filtered_files)-1:
                rate_3 += filtered_files[j] + ','
            else:
                rate_3 += filtered_files[j]
        output_text.write(rate_1 + ' ' + rate_2 + rate_3 + '\n')

        errs_1 = mjd #for the MJD
        errs_2 = '' #for the count rates in each band
        errs_3 = '' #for the list of files that goes into each bin
        for j in range(len(errs_all)):
            errs_2 += str(round(errs_all[j],4)) + ' '
        for j in range(len(filtered_files)):
            if j != len(filtered_files)-1:
                errs_3 += filtered_files[j] + ','
            else:
                errs_3 += filtered_files[j]
        output_text_err.write(errs_1 + ' ' + errs_2 + errs_3 + '\n')

    output_text.close()
    output_text_err.close()


if __name__ == "__main__":
    binned_text(bin_size)
