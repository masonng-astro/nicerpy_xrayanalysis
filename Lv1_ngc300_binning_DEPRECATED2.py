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
import glob
import Lv0_dirs

Lv0_dirs.global_par()

binsize_s = 10
bin_size = str(binsize_s).zfill(2) + 'd' #bins of 1 day!
bgsub_type = 'xbgsub'

#RGcms = Lv0_dirs.NGC300_2020 + 'n300_ulx.bgsub_cl50_g2020.fffphot' NOT USING for 2020. No fffphot file for this anyways
norm_fffphot = Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020norm.fffphot'

### extracting photometry/count rates for soft_1 (0.2-0.3 keV), soft_2 (0.3-0.4 keV),
### A (0.4-1 keV), B (1-2 keV), C (2-4 keV), D (4-12) keV, in-band (A+B+C+D = 0.4-12 keV) bands, plus MJD
time,norm = np.genfromtxt(norm_fffphot,dtype='float64',usecols=(11,12),unpack=True)
filenames = np.genfromtxt(norm_fffphot,dtype='str',usecols=(0),unpack=True)

time_bins = np.arange(58239,58606,binsize_s)
time = np.floor(time) #can't use int ; basically round the time values down to the integer value

def binned_text():
    """
    Given the MJDs, binned counts, and associated uncertainties, put them into a text file

    No arguments because I'll put all the bands in here
    """
    E_bins_low = np.array([20-1,30-1,40-1,100-1,200-1,400-1,40-1,1300-1])
    E_bins_high = np.array([30-1,40-1,100-1,200-1,400-1,1200-1,1200-1,1501-1])

    mjds_used = []
    rates_text = [] #row of rate values to put into the text file (each line = each time bin)
    errs_text = [] #row of error values to put into the text file (each line = each time bin)
    files_text = []
    for i in tqdm(range(len(time_bins))):
        ## so for each time bin, gather up the relevant bgsub.pha spectra, and extract rates/errors and combine them appropriately
        files_in_interval = filenames[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
        norms = norm[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]

        if len(files_in_interval) == 0:
            continue
        if bgsub_type == 'bgsub':
            filtered_files = [Lv0_dirs.NGC300_2020 + 'bgsub_cl50/pha/' + files_in_interval[a] for a in range(len(files_in_interval))]
        elif bgsub_type == 'xbgsub':
            filtered_files = [Lv0_dirs.NGC300_2020 + 'bgsub_cl50/xpha/' + files_in_interval[a] for a in range(len(files_in_interval))]

        files_text.append(filtered_files)
        mjds_used.append(time_bins[i])

        rates_row = []
        errs_row = []
        for j in range(len(E_bins_low)): #for each energy band
            exp_all = []
            rates_all = []
            errs_all = []

            for k in range(len(filtered_files)): #for each spectrum
                exp_all.append(fits.open(filtered_files[k])[1].header['EXPOSURE'])
                rates_all.append(sum(fits.open(filtered_files[k])[1].data['RATE'][E_bins_low[j]:E_bins_high[j]])*norms[k])
                errs = fits.open(filtered_files[k])[1].data['STAT_ERR'][E_bins_low[j]:E_bins_high[j]]
                errs_all.append( np.sqrt(np.mean(errs)**2 * len(errs)) * norms[k] )

            exp_all = np.array(exp_all)
            rates_all = np.array(rates_all)
            errs_all = np.array(errs_all)

            rates_row.append(sum(exp_all*rates_all)/sum(exp_all))
            errs_row.append(sum(exp_all*errs_all)/sum(exp_all))

        rates_text.append(rates_row)
        errs_text.append(errs_row)

    output_text = open(Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020norm_' + str(binsize_s).zfill(2) + 'd.fffphot','w')
    output_text_err = open(Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020err_norm_' + str(binsize_s).zfill(2) + 'd.fffphot','w')

    for i in range(len(rates_text)):
        rate_1 = str(mjds_used[i])
        rate_2 = ''
        for j in range(len(rates_text[i])):
            rate_2 += str(round(rates_text[i][j],4)) + ' '
        rate_3 = ''
        for j in range(len(files_text[i])):
            if j != len(files_text[i])-1:
                rate_3 += files_text[i][j] + ','
            else:
                rate_3 += files_text[i][j]
        output_text.write(rate_1 + ' ' + rate_2 + rate_3 + '\n')

    for i in range(len(errs_text)):
        errs_1 = str(mjds_used[i])
        errs_2 = ''
        for j in range(len(errs_text[i])):
            errs_2 += str(round(errs_text[i][j],4)) + ' '
        errs_3 = ''
        for j in range(len(files_text[i])):
            if j != len(files_text[i])-1:
                errs_3 += files_text[i][j] + ','
            else:
                errs_3 += files_text[i][j]
        output_text_err.write(errs_1 + ' ' + errs_2 + errs_3 + '\n')

    output_text.close()
    output_text_err.close()


if __name__ == "__main__":
    binned_text()
