#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 7 11:01am 2019

Automating the process of determining detection levels! Can also do calculations
pertaining to "z" (no. of Fourier bins signal drifts over) and acceleration!

Average acceleration of a pulsar in a binary system is a = z*c/(T^2 * f_0)

Based on van der Klis 1988 or 1989!

"""
from __future__ import division, print_function
import numpy as np
from scipy import stats, signal, special
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits

def max_acc(zmax,T,f0):
    """
    To obtain the maximum acceleration 'detectable' by PRESTO.

    zmax - (expected) maximum number of Fourier bins that the pulsar frequency
    f0 drifts
    T - observation duration (s)
    f0 - pulsar's frequency (Hz)
    """
    c = 299792458 #speed of light in m/s

    return zmax*c/(T**2*f0)

#print(max_acc(200,200,230))

def N_trials(tbin,T):
    """
    To obtain the number of trials used in the FFT. Divided by two to get number
    of trials f >= 0!

    tbin- size of the bins in time
    T - observation duration (s) or segment length (s)
    """

    return 1/2 * T/tbin

def single_trial_prob(significance,N):
    """
    To obtain the single trial probability required for a statistically significant
    "significance" detection, with N trials.

    significance - the number of 'sigmas' desired for detection
    N - number of trials
    """
    prob = 1-special.erf(significance/np.sqrt(2))
    single_trial = 1 - (1 - prob)**(1/N)
    single_trial_signif = special.erfinv(1-single_trial)*np.sqrt(2)

    print('The single trial probability required for a ' + str(significance) + ' sigma detection with ' + str(int(N)) + ' trials is ' + str(single_trial) + ', i.e., a significance of ' + str(single_trial_signif))

    return single_trial, single_trial_signif

#single_trial_prob(3,1e3)

def signal_significance(M,W,Pthreshold):
    """
    Calculating the significance of a particular signal in the power spectrum,
    given M (number of segments), W (number of consecutive bins summed), and
    Pthreshold (the power [Leahy-normalized] of the signal).

    M - number of segments
    W - number of consecutive bins summed
    Pthreshold - the power of the signal (Leahy-normalized)
    """
    chi2 = M*W*Pthreshold
    dof = 2*M*W

    Q_chi2_dof = 1-stats.chi2.cdf(chi2,dof) #stats.chi2 is from 0 to chi2, so do the complement if we want chi2 to infinity
    ## Q(M*W*Pthreshold|2*M*W) ; where Q(chi^2|nu) = 1/[2^(nu/2)*Gamma(nu/2)] * \int_{chi^2}^\infty t^{nu/2 - 1} e^{-t/2} dt
    significance = special.erfinv(1-Q_chi2_dof)*np.sqrt(2)
    print('The signal has a significance of ' + str(significance) + ' sigma.')
    #confidence_level = special.erf(significance/np.sqrt(2))*100
    return significance

sig_sig = signal_significance(45,1,3.90978)
single_trial_prob(3,2e6)