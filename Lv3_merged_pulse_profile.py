#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 5:18pm 2019

Script that takes in an input (barycentered) merged time series from multiple
ObsIDs, and outputs a pulse profile! Methods script.

"""
from __future__ import division, print_function
import numpy as np
import subprocess
from tqdm import tqdm
from astropy.io import fits

import matplotlib.pyplot as plt
import Lv0_dirs,Lv2_presto_segments,Lv2_average_merge_ps_methods,Lv2_merged_pulse_methods

preprocessing = False
merging = True
merged_id = '000008'
#obsids = ['1030180' + str(i) for i in range(101,131)] + ['1030180' + str(i) for i in range(149,188)]
#obsids = ['1030180' + str(i) for i in range(101,188)]
#obsids.remove('1030180113')
#obsids.remove('1030180164')

"""
obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113'] + [str(i) for i in range(1060060101,1060060200)] + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
obsids.remove('0060060113')
obsids.remove('1060060116')
obsids.remove('1060060117')
obsids.remove('1060060118')
obsids.remove('1060060119')
obsids.remove('1060060120')
obsids.remove('1060060121')
obsids.remove('1060060122')
obsids.remove('1060060125')
obsids.remove('1060060126')
obsids.remove('1060060232')
obsids.remove('1060060233')
obsids.remove('1060060254')
obsids.remove('1060060258')
obsids.remove('1060060274')
"""

#obsids = [str(i) for i in range(1050230101,1050230108)]

par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1231-1411.par' #parameter file for demodulation
#par_file = Lv0_dirs.NICERSOFT_DATADIR + 'B1957+20.par'
#par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1757-2508.par'

E_trunc = True
PI1 = 35 #lower bound for PI
PI2 = 150 #upper bound for PI
#pulse_pars = [182.06580377,0,0] #J1756.9-2508
#pulse_pars = [29.639575,-3.77535E-10,1.1147E-20] #Crab
pulse_pars = [271.453019624388,-1.66705E-15,0] #J1231
#pulse_pars = [622.12202499673376738,-6.5114520166830696246e-15,0] #B1957
no_phase_bins = 50

if preprocessing == True:
    Lv2_average_merge_ps_methods.merging(obsids)
    Lv2_average_merge_ps_methods.merging_GTIs(obsids,merged_id)

if merging == True:
    Lv2_merged_pulse_methods.niextract_gti_energy(merging,merged_id,PI1,PI2)
    Lv2_merged_pulse_methods.do_demodulate(merging,merged_id,par_file,E_trunc,PI1,PI2)
    Lv2_merged_pulse_methods.plot_pf(merging,merged_id,E_trunc,PI1,PI2,pulse_pars,no_phase_bins)

average_profile = np.zeros(no_phase_bins*2)
if merging == False:
    for i in tqdm(range(len(obsids))):
        Lv2_merged_pulse_methods.niextract_gti_energy(merging,obsids[i],PI1,PI2)
        Lv2_merged_pulse_methods.do_demodulate(merging,obsids[i],par_file,E_trunc,PI1,PI2)
        phase_bins,summed_profile = Lv2_merged_pulse_methods.pulse_profile(merging,obsids[i],E_trunc,PI1,PI2,pulse_pars,no_phase_bins)
        average_profile += summed_profile

    plt.figure()
    plt.step(phase_bins[:-1],(average_profile/len(obsids))/np.mean(average_profile/len(obsids)))
    #plt.plot(phase_bins[:-1],(average_profile/len(obsids))/np.mean(average_profile/len(obsids)),'r+')
    plt.xlabel('Phase', fontsize=12)
    plt.ylabel('Count/s',fontsize=12)
    plt.show()

#################################################################################
