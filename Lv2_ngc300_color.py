#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thurs Aug 8th 10:32am 2019

Calculating colors from the customized binned data with NGC300 ULX

"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import Lv0_dirs

Lv0_dirs.global_par()

bin_size = '03d' #bins of 1 day!

def get_color(bin_size,band1,band2):
    """
    Obtain colors and the corresponding uncertainties. Will NOT use values where
    either the counts/rate from band1 OR band2 are negative! Allowed band values are
    "soft1, soft2, A, B, C, D, and inband."

    bin_size - binning size desired (1 day or 10 days, for example)
    band1 - energy band 1
    band2 - energy band 2
    """
    if band1 != 'soft1' and band1 != 'soft2' and band1 != 'A' and band1 != 'B' and band1 != 'C' and band1 != 'D' and band1 != 'inband':
        raise ValueError("Make sure band1 is either of soft1, soft2, A, B, C, D, or inband!")
    if band2 != 'soft1' and band2 != 'soft2' and band2 != 'A' and band2 != 'B' and band2 != 'C' and band2 != 'D' and band2 != 'inband':
        raise ValueError("Make sure band2 is either of soft1, soft2, A, B, C, D, or inband!")

    binned_counts_file = Lv0_dirs.NGC300 + 'n300_ulx.bgsub_cl50_RGnorm_' + bin_size + '.ffphot'
    binned_unc_file = Lv0_dirs.NGC300 + 'n300_ulx.bgsub_cl50_RGerr_' + bin_size + '.ffphot'

    mjds = np.genfromtxt(binned_counts_file,usecols=(0),unpack=True)
    if band1 == 'soft1':
        counts_band1 = np.genfromtxt(binned_counts_file,usecols=(1),unpack=True)
        unc_band1 = np.genfromtxt(binned_unc_file,usecols=(1),unpack=True)
    if band1 == 'soft2':
        counts_band1 = np.genfromtxt(binned_counts_file,usecols=(2),unpack=True)
        unc_band1 = np.genfromtxt(binned_unc_file,usecols=(2),unpack=True)
    if band1 == 'A':
        counts_band1 = np.genfromtxt(binned_counts_file,usecols=(3),unpack=True)
        unc_band1 = np.genfromtxt(binned_unc_file,usecols=(3),unpack=True)
    if band1 == 'B':
        counts_band1 = np.genfromtxt(binned_counts_file,usecols=(4),unpack=True)
        unc_band1 = np.genfromtxt(binned_unc_file,usecols=(4),unpack=True)
    if band1 == 'C':
        counts_band1 = np.genfromtxt(binned_counts_file,usecols=(5),unpack=True)
        unc_band1 = np.genfromtxt(binned_unc_file,usecols=(5),unpack=True)
    if band1 == 'D':
        counts_band1 = np.genfromtxt(binned_counts_file,usecols=(6),unpack=True)
        unc_band1 = np.genfromtxt(binned_unc_file,usecols=(6),unpack=True)
    if band1 == 'inband':
        counts_band1 = np.genfromtxt(binned_counts_file,usecols=(7),unpack=True)
        unc_band1 = np.genfromtxt(binned_unc_file,usecols=(7),unpack=True)

    if band2 == 'soft1':
        counts_band2 = np.genfromtxt(binned_counts_file,usecols=(1),unpack=True)
        unc_band2 = np.genfromtxt(binned_unc_file,usecols=(1),unpack=True)
    if band2 == 'soft2':
        counts_band2 = np.genfromtxt(binned_counts_file,usecols=(2),unpack=True)
        unc_band2 = np.genfromtxt(binned_unc_file,usecols=(2),unpack=True)
    if band2 == 'A':
        counts_band2 = np.genfromtxt(binned_counts_file,usecols=(3),unpack=True)
        unc_band2 = np.genfromtxt(binned_unc_file,usecols=(3),unpack=True)
    if band2 == 'B':
        counts_band2 = np.genfromtxt(binned_counts_file,usecols=(4),unpack=True)
        unc_band2 = np.genfromtxt(binned_unc_file,usecols=(4),unpack=True)
    if band2 == 'C':
        counts_band2 = np.genfromtxt(binned_counts_file,usecols=(5),unpack=True)
        unc_band2 = np.genfromtxt(binned_unc_file,usecols=(5),unpack=True)
    if band2 == 'D':
        counts_band2 = np.genfromtxt(binned_counts_file,usecols=(6),unpack=True)
        unc_band2 = np.genfromtxt(binned_unc_file,usecols=(6),unpack=True)
    if band2 == 'inband':
        counts_band2 = np.genfromtxt(binned_counts_file,usecols=(7),unpack=True)
        unc_band2 = np.genfromtxt(binned_unc_file,usecols=(7),unpack=True)

    counts_band1_pos = counts_band1[(counts_band1>0)&(counts_band2>0)]
    counts_band2_pos = counts_band2[(counts_band1>0)&(counts_band2>0)]
    unc_band1_pos = unc_band1[(counts_band1>0)&(counts_band2>0)]
    unc_band2_pos = unc_band2[(counts_band1>0)&(counts_band2>0)]
    mjds_pos = mjds[(counts_band1>0)&(counts_band2>0)]

    color = counts_band2_pos/counts_band1_pos
    color_unc = np.sqrt( (unc_band2_pos/counts_band1_pos)**2 + (counts_band2_pos*unc_band1_pos/counts_band1_pos**2)**2 )

    return mjds_pos, color, color_unc

if __name__ == "__main__":
    binned_text()
