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
import subprocess
import glob
import Lv0_dirs

Lv0_dirs.global_par()

binsize_s = 5
bin_size = str(binsize_s).zfill(2) + 'd' #bins of 1 day!
bgsub_type = 'xbgsub'

#RGcms = Lv0_dirs.NGC300_2020 + 'n300_ulx.bgsub_cl50_g2020.fffphot' NOT USING for 2020. No fffphot file for this anyways
norm = Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020norm.fffphot'
error = Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020err_norm.fffphot'

### extracting photometry/count rates for soft_1 (0.2-0.3 keV), soft_2 (0.3-0.4 keV),
### A (0.4-1 keV), B (1-2 keV), C (2-4 keV), D (4-12) keV, in-band (A+B+C+D = 0.4-12 keV) bands, plus MJD
soft1,soft2,A_band,B_band,C_band,D_band,inband,time = np.genfromtxt(norm,usecols=(3,4,5,6,7,8,9,11),unpack=True)
soft1_err,soft2_err,A_err,B_err,C_err,D_err,inband_err,time_error = np.genfromtxt(error,usecols=(3,4,5,6,7,8,9,11),unpack=True)
# time_error does not mean error in the time value, I just mean the time value found in the error text file

### extracting file names
data_filename = np.genfromtxt(norm,dtype='str',usecols=(0),unpack=True)
data_filename_err = np.genfromtxt(error,dtype='str',usecols=(0),unpack=True)

if bgsub_type == 'bgsub':
    filename = np.array([Lv0_dirs.NGC300_2020 + 'bgsub_cl50/pha/' + data_filename[i] for i in range(len(data_filename))])
    filename_err = np.array([Lv0_dirs.NGC300_2020 + 'bgsub_cl50/pha/' + data_filename_err[i] for i in range(len(data_filename_err))])
elif bgsub_type == 'xbgsub':
    filename = np.array([Lv0_dirs.NGC300_2020 + 'bgsub_cl50/xpha/' + data_filename[i] for i in range(len(data_filename))])
    filename_err = np.array([Lv0_dirs.NGC300_2020 + 'bgsub_cl50/xpha/' + data_filename_err[i] for i in range(len(data_filename_err))])

time_bins = np.arange(58239,58606,binsize_s)
filename_dict,soft1_dict,soft2_dict,A_dict,B_dict,C_dict,D_dict,inband_dict = {},{},{},{},{},{},{},{}
filename_err_dict,soft1_err_dict,soft2_err_dict,A_err_dict,B_err_dict,C_err_dict,D_err_dict,inband_err_dict = {},{},{},{},{},{},{},{}
time = np.floor(time) #can't use int ; basically round the time values down to the integer value

#### bin up the counts and associated errors into dictionaries!
for i in tqdm(range(len(time_bins))):
    ## so for each time bin, gather up the count rates (and errors) in a dictionary
    filename_dict[time_bins[i]] = filename[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    soft1_dict[time_bins[i]] = soft1[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    soft2_dict[time_bins[i]] = soft2[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    A_dict[time_bins[i]] = A_band[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    B_dict[time_bins[i]] = B_band[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    C_dict[time_bins[i]] = C_band[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    D_dict[time_bins[i]] = D_band[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    inband_dict[time_bins[i]] = inband[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]

    filename_err_dict[time_bins[i]] = filename_err[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    soft1_err_dict[time_bins[i]] = soft1_err[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    soft2_err_dict[time_bins[i]] = soft2_err[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    A_err_dict[time_bins[i]] = A_err[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    B_err_dict[time_bins[i]] = B_err[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    C_err_dict[time_bins[i]] = C_err[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    D_err_dict[time_bins[i]] = D_err[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]
    inband_err_dict[time_bins[i]] = inband_err[(time>=time_bins[i])&(time<time_bins[i]+binsize_s)]

### doing the calculations for weighted averages and associated uncertainties
def get_binned_data(counts_dict,err_dict):
    """
    Given the dictionaries where the data (counts and uncertainty) are already
    binned, put them into lists!

    counts_dict - dictionary where the keys are MJD values, and the entries correspond
    to the counts/rate for a given MJD
    err_dict - dictionary where the keys are MJD values, and the entries correspond
    to the UNCERTAINTY in the counts/rate for a given MJD
    """
    if type(counts_dict) != dict and type(err_dict) != dict:
        raise TypeError("Make sure that counts_dict and err_dict are actually dictionaries!")

    binned_MJD = []
    binned_counts = []
    binned_uncertainty = []
    binned_pha_files = []

    dict_keys = sorted(counts_dict.keys())
    for i in range(len(dict_keys)):
        counts = counts_dict[dict_keys[i]] #get count rates corresponding to a given MJD
        if len(counts) != 0:
            pha_string = ''
            pha_list = filename_dict[dict_keys[i]]
            for j in range(len(pha_list)):
                if j != len(pha_list)-1:
                    pha_string += pha_list[j] + ','
                else:
                    pha_string += pha_list[j]
            binned_pha_files.append(pha_string)

            unc = err_dict[dict_keys[i]] #get errors corresponding to a given MJD
            weights = 1/unc**2 #calculating the weights
            weighted_ave = sum(weights*counts)/sum(weights) #weighted average
            weighted_ave_unc = 1/np.sqrt(sum(weights)) #uncertainty for the weighted average
            binned_MJD.append(dict_keys[i])
            binned_counts.append(weighted_ave)
            binned_uncertainty.append(weighted_ave_unc)

    return binned_pha_files,binned_MJD,binned_counts,binned_uncertainty

def binned_text():
    """
    Given the MJDs, binned counts, and associated uncertainties, put them into a text file

    No arguments because I'll put all the bands in here
    """
    binned_pha_files,binned_MJD,binned_counts_soft1,binned_unc_soft1 = get_binned_data(soft1_dict,soft1_err_dict)
    binned_pha_files,binned_MJD,binned_counts_soft2,binned_unc_soft2 = get_binned_data(soft2_dict,soft2_err_dict)
    binned_pha_files,binned_MJD,binned_counts_A,binned_unc_A = get_binned_data(A_dict,A_err_dict)
    binned_pha_files,binned_MJD,binned_counts_B,binned_unc_B = get_binned_data(B_dict,B_err_dict)
    binned_pha_files,binned_MJD,binned_counts_C,binned_unc_C = get_binned_data(C_dict,C_err_dict)
    binned_pha_files,binned_MJD,binned_counts_D,binned_unc_D = get_binned_data(D_dict,D_err_dict)
    binned_pha_files,binned_MJD,binned_counts_inband,binned_unc_inband = get_binned_data(inband_dict,inband_err_dict)

    counts_file = Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020norm_' + bin_size + '.fffphot'
    output_file = open(counts_file,'w')

    ### get MJD (int), soft1, soft2, A, B, C, D, inband, all associated pha files
    for i in range(len(binned_MJD)):
        output_file.write(str(binned_MJD[i]) + ' ' + str(round(binned_counts_soft1[i],4)) + ' ' + str(round(binned_counts_soft2[i],4)) + ' ' + str(round(binned_counts_A[i],4)) + ' ' + str(round(binned_counts_B[i],4)) + ' ' + str(round(binned_counts_C[i],4)) + ' ' + str(round(binned_counts_D[i],4)) + ' ' + str(round(binned_counts_inband[i],4)) + ' ' + binned_pha_files[i] + '\n')
    output_file.close()

    unc_file = Lv0_dirs.NGC300_2020 + 'n300_ulx.' + bgsub_type + '_cl50_g2020err_norm_' + bin_size + '.fffphot'
    output_file = open(unc_file,'w')

    for i in range(len(binned_MJD)):
        output_file.write(str(binned_MJD[i]) + ' ' + str(round(binned_unc_soft1[i],4)) + ' ' + str(round(binned_unc_soft2[i],4)) + ' ' + str(round(binned_unc_A[i],4)) + ' ' + str(round(binned_unc_B[i],4)) + ' ' + str(round(binned_unc_C[i],4)) + ' ' + str(round(binned_unc_D[i],4)) + ' ' + str(round(binned_unc_inband[i],4)) + ' ' + binned_pha_files[i] + '\n')
    output_file.close()


if __name__ == "__main__":
    binned_text()
