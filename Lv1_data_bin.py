#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 8 2:11pm 2019

Extracting the GTIs from the FITS files. Use the event_cl files.

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv0_call_eventcl
from scipy import stats
import matplotlib.pyplot as plt

def binning(obsid,bary,par_list,tbin_size):
    """
    Binning routine. Got to make sure I have TIME and PI called!

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(tbin_size) != int and type(tbin_size) != np.float:
        raise TypeError("tbin_size should be a float or integer!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if 'PI' and 'TIME' not in par_list:
        raise ValueError("You should have BOTH 'PI' and 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    data_dict = Lv0_call_eventcl.get_eventcl(obsid,bary,par_list)

    times = data_dict['TIME']
    PI_data = data_dict['PI']

    shifted_t = times-times[0]
    counts_data = np.ones(len(PI_data))

    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])*1/tbin_size+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts_data,statistic='sum',bins=t_bins)

    return t_bins, summed_data
