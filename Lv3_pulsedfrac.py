#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 8 11:20am 2019

Obtaining the pulsed fraction in terms of fractional rms amplitude.

"""
from __future__ import division, print_function
import numpy as np
import Lv2_phase
from scipy import stats
from PyAstronomy.pyasl import foldAt
import matplotlib.pyplot as plt
import os

def pf_all(obsid,bary,par_list,tbin_size,f_pulse,shift,no_phase_bins,mode):
    """
    Obtaining the pulsed fraction from the entire raw pulse profile without
    any cuts to the data.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    f_pulse - the frequency of the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how the pulse profile is 'displaced'.
    no_phase_bins - number of phase bins desired
    mode - whether we want to show or save the plot.
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    phase,counts = Lv2_phase.whole(obsid,bary,par_list,tbin_size,f_pulse,shift,no_phase_bins,mode)

    phase_ave_mean = np.mean(counts)
    variance = np.var(counts) #gives the same result as me writing out the formula manually! np.var is ok.

    return np.sqrt(variance)/phase_ave_mean

def pf_t(obsid,bary,par_list,tbin_size,f_pulse,shift,no_phase_bins,t1,t2,mode):
    """
    Obtain the pulsed fraction from the pulse profile for a desired time interval.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    f_pulse - the frequency of the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how it is presented.
    no_phase_bins - number of phase bins desired
    t1 - lower time boundary
    t2 - upper time boundary
    mode - whether we want to show or save the plot
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if t2<t1:
        raise ValueError("t2 should be greater than t1!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    phase,counts = Lv2_phase.partial_t(obsid,bary,par_list,tbin_size,f_pulse,shift,no_phase_bins,t1,t2,mode)

    phase_ave_mean = np.mean(counts)
    variance = np.var(counts)

    return np.sqrt(variance)/phase_ave_mean

def pf_E(obsid,bary,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,E1,E2,mode):
    """
    Obtain the pulsed fraction from the pulse profile for a desired energy range.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means by in 0.05s
    Ebin_size - the size of the energy bins (in keV!)
    >> e.g., Ebin_size = 0.1 means bin by 0.1keV
    >> e.g., Ebin_size = 0.01 means bin by 0.01keV!
    f_pulse - the frequency of the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how it is presented.
    no_phase_bins - number of phase bins desired
    E1 - lower energy boundary
    E2 - upper energy boundary
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if E2<E1:
        raise ValueError("E2 should be greater than E1!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    phase,counts = Lv2_phase.partial_E(obsid,bary,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,E1,E2,mode)

    phase_ave_mean = np.mean(counts)
    variance = np.var(counts)

    return np.sqrt(variance)/phase_ave_mean

def pf_tE(obsid,bary,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,t1,t2,E1,E2,mode):
    """
    Obtain the pulsed fraction from a desired time interval and desired energy range.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means by in 0.05s
    Ebin_size - the size of the energy bins (in keV!)
    >> e.g., Ebin_size = 0.1 means bin by 0.1keV
    >> e.g., Ebin_size = 0.01 means bin by 0.01keV!
    f_pulse - the frequency of the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how it is presented.
    no_phase_bins - number of phase bins desired
    t1 - lower time boundary
    t2 - upper time boundary
    E1 - lower energy boundary
    E2 - upper energy boundary
    mode - whether we want to show or save the plot
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if E2<E1:
        raise ValueError("E2 should be greater than E1!")
    if t2<t1:
        raise ValueError("t2 should be greater than t1!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    phase,counts = Lv2_phase.partial_tE(obsid,bary,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,t1,t2,E1,E2,mode)

    phase_ave_mean = np.mean(counts)
    variance = np.var(counts)

    return np.sqrt(variance)/phase_ave_mean

print(pf_tE('0034070101',True,['TIME','PI','PI_FAST'],1,0.05,0.20846118761251825,0.4,50,0,500,0.3,12,'show'))
