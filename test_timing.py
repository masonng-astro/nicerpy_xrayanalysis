#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 7:29am 2019

Program for timing

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv1_data_gtis,Lv2_mkdir
from PyAstronomy.pyasl import foldAt
from scipy import stats
from scipy.optimize import curve_fit
from tqdm import tqdm
import matplotlib.pyplot as plt
import time
import pathlib

Lv0_dirs.global_par()

def linear_f(t,phi0,f0,T0):
    return (phi0 + f0*(t-T0))%1

def get_phases(eventfile,f_pulse,shift):
    """
    Folding the time series by the pulse frequency to obtain phase values

    eventfile - path to the event file.
    f_pulse - the frequency of the pulse
    shift - how much to shift the pulse by in the phase axis.
    """
    raw_times = fits.open(eventfile)[1].data['TIME']
    times = raw_times - raw_times[0]

    period = 1/f_pulse
    phases = foldAt(times,period,T0=shift*period)

    return phases

def fit_to_linear(eventfile,f_pulse,shift,T0):
    """
    Fitting the phases to a linear phase model

    eventfile - path to the event file.
    f_pulse - the frequency of the pulse
    shift - how much to shift the pulse by in the phase axis.
    T0 - some reference T0 time
    """
    raw_times = fits.open(eventfile)[1].data['TIME']
    times = raw_times - raw_times[0]

    phases = get_phases(eventfile,f_pulse,shift)
    popt,pcov = curve_fit(linear_f,times,phases,p0=[shift,f_pulse,T0],bounds=([0.3,0.2085,40],[0.5,0.2095,60]))

    return popt,np.sqrt(np.diag(pcov))

if __name__ == "__main__":
    eventfile = Lv0_dirs.NICERSOFT_DATADIR + '0034070101_pipe/ni0034070101_nicersoft_bary.evt'
    raw_times = fits.open(eventfile)[1].data['TIME']
    times = raw_times - raw_times[0]

    phases = get_phases(eventfile,0.20879991,0.35705417)

    #phase_bins = np.linspace(0,1,21)
    #summed_profile, bin_edges, binnumber = stats.binned_statistic(phases,np.ones(len(phases)),statistic='sum',bins=phase_bins)

    vars,unc = fit_to_linear(eventfile,0.2088,0.366,50)

    print(vars)
    #plt.plot(times,linear_f(times,vars[0],vars[1],vars[2]),'rx')
    resid_phases = phases - linear_f(times,vars[0],vars[1],vars[2])

    resid_phases_plot = resid_phases[resid_phases>=-0.5]
    times_plot = times[resid_phases>=-0.5]
    plt.plot(times_plot,resid_phases_plot)

    #plt.plot(phase_bins[:-1],summed_profile,'rx')
    plt.show()
