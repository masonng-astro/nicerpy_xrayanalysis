#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Friday Jan 11 11:22am 2019

Getting the power spectra from the periodogram or manual method

"""
from __future__ import division, print_function
import numpy as np
from scipy import stats, signal
import matplotlib.pyplot as plt

def padding(counts):
    """
    For use in the function manual. Recall: The optimal number of bins is 2^n,
    n being some natural number. We pad 0s onto the original data set, where
    the number of 0s to pad is determined by the difference between the optimal
    number of bins and the length of the data set (where the former should be
    greater than the latter!)

    counts - array of counts from the binned data
    """
    if type(counts) != list and type(counts) != np.ndarray:
        raise TypeError("counts should either be a list or an array!")

    data_size = len(counts)
    diff = [np.abs(data_size-2**n) for n in range(0,30)]
    min_diff_index = np.argmin(diff)

    optimal_bins = 2**(min_diff_index)
    if optimal_bins < data_size: #if the number of optimal bins for FFT is less
    #than that of the number of data bins in the data set, *2 to get 2^N bins
    #that will yield it
        optimal_bins *= 2

    ### padding the data set
    len_diff = optimal_bins - data_size
    if len_diff < 0:
        raise ValueError("For some reason, the number of optimal bins is less than that of the data size, CHECK program!")

    padded_counts = list(counts) + list(np.zeros(len_diff))
    padded_counts = np.array(padded_counts)

    return padded_counts

def oversample(factor,counts):
    """
    Perform oversampling on the data. Return the padded array of counts.

    factor - N-times oversampling; factor = 5 means 5x oversampling
    counts - array of counts from the binned data
    """
    if type(factor) != int:
        raise TypeError("Make sure the second entry in the array is an integer!")
    pad_zeros = np.zeros(len(counts)*(factor-1))
    oversampled_counts = np.array(list(counts) + list(pad_zeros))
    padded_counts = padding(oversampled_counts)

    return padded_counts

def pdgm(times,counts,xlims,vlines,toplot,oversampling):
    """
    Generating the power spectrum through the signal.periodogram method.

    times - array of binned times
    counts - array of counts from the binned data
    xlims - a list or array: first entry = True/False as to whether to impose an
    xlim; second and third entry correspond to the desired x-limits of the plot
    vlines - a list or array: first entry = True/False as to whether to draw
    a vertical line in the plot; second entry is the equation for the vertical line
    toplot - whether to show the plot or not
    oversampling - whether to perform oversampling. Array will consist of
    [True/False, oversampling factor]
    """
    if type(times) != list and type(times) != np.ndarray:
        raise TypeError("times should either be a list or an array!")
    if type(counts) != list and type(counts) != np.ndarray:
        raise TypeError("counts should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")
    if toplot != True and toplot != False:
        raise ValueError("toplot should either be True or False!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")

    T = times[-1] - times[0]
    dt = T/len(times)
    freq = 1/dt
    f,pxx = signal.periodogram(counts,fs=freq)

    if oversampling[0] == True:
        os_counts = oversample(oversampling[1],counts)
        dt = T/len(times)
        freq = 1/dt
        f,pxx = signal.periodogram(os_counts,fs=freq)

    if toplot == True:
        plt.semilogy(f[1:],pxx[1:])#/np.mean(pxx[1:]))
        plt.xlabel('Hz',fontsize=12)
        plt.ylabel('Normalized power spectrum',fontsize=12)
        if xlims[0] == True:
            plt.xlim([xlims[1],xlims[2]])
        if vlines[0] == True:
            plt.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)

    return f[1:],pxx[1:]

def manual(times,counts,xlims,vlines,toplot,oversampling):
    """
    Generating the power spectrum through the manual FFT method.

    times - array of binned times
    counts - array of counts from the binned data
    xlims - a list or array: first entry = True/False as to whether to impose an
    xlim; second and third entry correspond to the desired x-limits of the plot
    vlines - a list or array: first entry = True/False as to whether to draw
    a vertical line in the plot; second entry is the equation for the vertical line
    toplot - whether to show the plot or not
    oversampling - whether to perform oversampling. Array will consist of
    [True/False, oversampling factor]
    """
    if type(times) != list and type(times) != np.ndarray:
        raise TypeError("times should either be a list or an array!")
    if type(counts) != list and type(counts) != np.ndarray:
        raise TypeError("counts should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")
    if toplot != True and toplot != False:
        raise ValueError("toplot should either be True or False!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")

    T = times[-1] - times[0]
    dt = T/len(times)

    padded_counts = padding(counts)
    mean_corrected = padded_counts-np.mean(counts)
    power_spec = 2.0/sum(counts)*np.abs(np.fft.fft(mean_corrected))**2 #sum(counts) = N_photons
    freqs = np.fft.fftfreq(padded_counts.size,dt)

    if oversampling[0] == True:
        os_counts = oversample(oversampling[1],counts)
        dt = T/len(times)

        padded_counts = padding(os_counts)
        mean_corrected = padded_counts-np.mean(counts)
        power_spec = 2.0/sum(counts)*np.abs(np.fft.fft(mean_corrected))**2 #sum(counts) = N_photons
        freqs = np.fft.fftfreq(padded_counts.size,dt)

    N = len(freqs)
#    checkmean = power_spec[(freqs>=10)&(freqs<=50)]
#    print(np.mean(checkmean))
    ####freqs = np.linspace(0,1/(2.0*dt),int(len(test_times)/2)) #correct method! Same as that of fftfreq
    if toplot == True:
        plt.semilogy(freqs[1:int(N/2)],power_spec[1:int(N/2)],'rx-')#/np.mean(power_spec[1:int(N/2)]),'rx-')
        plt.xlabel('Hz')
        plt.ylabel('Normalized power spectrum')
        if xlims[0] == True:
            plt.xlim([xlims[1],xlims[2]])
        if vlines[0] == True:
            plt.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)

    return freqs[1:int(N/2)], power_spec[1:int(N/2)]
