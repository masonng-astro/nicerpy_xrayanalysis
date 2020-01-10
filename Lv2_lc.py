#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thurs Jan 10 3:02pm 2019

Updated on Mon Jun 3 - Added name_par_list for NICERsoft segments

Plotting light curves

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv0_fits2dict,Lv1_data_bin
from scipy import stats
import pathlib
import matplotlib.pyplot as plt
import os

Lv0_dirs.global_par() #obtaining the global parameters

def whole(eventfile,par_list,tbin_size,mode):
    """
    Plot the entire raw time series without any cuts to the data.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    mode - whether we want to show or save the plot.
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,par_list)
    times = data_dict['TIME']
    counts = np.ones(len(times))

    shifted_t = times-times[0]
    t_bins = np.linspace(0,np.ceil(shifted_t[-1]),np.ceil(shifted_t[-1])*1/tbin_size+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    plt.plot(t_bins[:-1],summed_data)
    plt.title('Light curve for ' + obj_name + ', ObsID ' + str(obsid),fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)
    if mode == 'show':
        plt.show()
    elif mode == 'save':
        filename = 'lc_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        plt.savefig(parent_folder+filename,dpi=900)
        plt.close()

def partial_t(eventfile,par_list,tbin_size,t1,t2,mode):
    """
    Plot the time series for a desired time interval.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    t1 - lower time boundary
    t2 - upper time boundary
    mode - whether we want to show or save the plot
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if t2<t1:
        raise ValueError("t2 should be greater than t1!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    truncated_t, truncated_counts = Lv1_data_bin.binning_t(eventfile,par_list,tbin_size,t1,t2)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    plt.plot(truncated_t[:-1], truncated_counts)
    plt.title('Light curve for ' + obj_name + ', ObsID ' + str(obsid) + '\n Time interval: ' + str(t1) + 's - ' + str(t2) + 's',fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        filename = 'lc_' + obsid + '_bin' + str(tbin_size) + 's_' + str(t1) + 's-' + str(t2) + 's.pdf'
        plt.savefig(parent_folder+filename,dpi=900)
        plt.close()

def partial_E(eventfile,par_list,tbin_size,E1,E2,mode):
    """
    Plot the time series for a desired energy range.
    [Though I don't think this will be used much. Count/s vs energy is pointless,
    since we're not folding in response matrix information here to get the flux.
    So we're just doing a count/s vs time with an energy cut to the data.]

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means by in 0.05s
    E1 - lower energy boundary
    E2 - upper energy boundary
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if E2<E1:
        raise ValueError("E2 should be greater than E1!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_E(eventfile,par_list,tbin_size,0.05,E1,E2)
    #put in Ebin_size of 0.05 keV; I'm plotting light curves in this script, so truncated_E and truncated_E_counts are not needed

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    plt.figure()
    plt.plot(truncated_t[:-1], truncated_t_counts)
    plt.title('Light curve for ' + obj_name + ', ObsID ' + str(obsid)+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        filename = 'lc_' + obsid + '_bin' + str(tbin_size) + 's_' + str(E1) + 'keV-' + str(E2) + 'keV.pdf'
        plt.savefig(parent_folder+filename,dpi=900)
        plt.close()

def partial_tE(eventfile,par_list,tbin_size,t1,t2,E1,E2,mode):
    """
    Plot the time series for a desired time interval and desired energy range.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means by in 0.05s
    t1 - lower time boundary
    t2 - upper time boundary
    E1 - lower energy boundary
    E2 - upper energy boundary
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
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

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_tE(eventfile,par_list,tbin_size,0.05,t1,t2,E1,E2)
    #put in Ebin_size of 0.05 keV; I'm plotting light curves in this script, so truncated_E and truncated_E_counts are not needed

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    plt.figure()
    plt.plot(truncated_t[:-1], truncated_t_counts)
    plt.title('Light curve for ' + obj_name + ', ObsID ' + str(obsid)+ '\n Time interval: ' + str(t1) + 's - ' + str(t2) + 's'+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        filename = 'lc_' + obsid + '_bin' + str(tbin_size) + 's_' + str(t1) + 's-' + str(t2) + 's_' + str(E1) + 'keV-' + str(E2) + 'keV.pdf'
        plt.savefig(parent_folder+filename,dpi=900)
        plt.close()

if __name__ == "__main__":
    eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034070101_pipe/cleanfilt.evt'

    whole(eventfile,['TIME','PI','PI_FAST'],1,'save')
    partial_t(eventfile,['TIME','PI','PI_FAST'],1,0,400,'save')
    partial_E(eventfile,['TIME','PI','PI_FAST'],1,0.3,6,'save')
    partial_tE(eventfile,['TIME','PI','PI_FAST'],1,0,200,0.3,6,'save')
