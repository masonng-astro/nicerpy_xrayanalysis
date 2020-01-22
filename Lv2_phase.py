#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 9:56am 2019

Plotting phase curves

Updated on Mon Jun 3 - Added name_par_list for NICERsoft segments

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv0_fits2dict,Lv1_data_bin,Lv2_mkdir
from scipy import stats
from PyAstronomy.pyasl import foldAt
import pathlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import os

Lv0_dirs.global_par() #obtaining the global parameters

def pulse_profile(f_pulse,times,counts,shift,no_phase_bins):
    """
    Calculating the pulse profile for the observation. Goes from 0 to 2!
    Thoughts on 1/14/2020: I wonder if the count rate is calculated from times[-1]-times[0]?
    If so, this is WRONG! I should be using the total from the GTIs!

    f_pulse - the frequency of the pulse
    times - the array of time values
    counts - the array of counts values
    shift - how much to shift the pulse by in the phase axis.
    It only affects how it is presented.
    no_phase_bins - number of phase bins desired
    """
    period = 1/f_pulse
    phases = foldAt(times,period,T0=shift*period)

    index_sort = np.argsort(phases)
    phases = list(phases[index_sort]) + list(phases[index_sort]+1)
    counts = list(counts[index_sort])*2

    phase_bins = np.linspace(0,2,no_phase_bins*2+1)
    summed_profile, bin_edges, binnumber = stats.binned_statistic(phases,counts,statistic='sum',bins=phase_bins)

    return phases, phase_bins, summed_profile

def pulse_folding(t,T,T0,f,fdot,fdotdot,no_phase_bins):
    """
    Calculating the pulse profile by also incorporating \dot{f} corrections!
    Goes from 0 to 2.

    t - array of time values
    T - sum of all the GTIs
    T0 - reference epoch in MJD
    f - pulse/folding Frequency
    fdot - frequency derivative
    fdotdot - second derivative of frequency
    no_phase_bins - number of phase bins desired (recommended 20!)

    Returns the pulse profile in counts/s/phase bin vs phase. The number of counts
    is divided by the exposure time (calculated through total sum of the GTIs)

    Also added a "TIMEZERO" manually in the script since it'd be inconvenient to call the eventfile here.
    """
    MJDREFI = 56658
    MJDREFF = 0.000777592592592593
    TIMEZERO = -1
    t_MJDs =  MJDREFI + MJDREFF + (TIMEZERO+t)/86400

    tau = (t_MJDs-T0)*86400
    phase = (f*tau + fdot/2 *tau**2 + fdotdot/6*tau**3)%1

    counts = np.ones(len(phase))
    phase_bins = np.linspace(0,1,no_phase_bins+1)

    summed_profile,bin_edges,binnumber = stats.binned_statistic(phase,counts,statistic='sum',bins=phase_bins)

    phase_bins_pad = np.linspace(1,2,no_phase_bins+1)
    summed_profile_pad = summed_profile

    phase_bins_total = np.array(list(phase_bins[:-1]) + list(phase_bins_pad))
    summed_profile_total = np.array(list(summed_profile) + list(summed_profile_pad))

    return phase_bins_total, summed_profile_total/T

def whole(eventfile,par_list,tbin_size,pulse_pars,shift,no_phase_bins,mode):
    """
    Plot the entire raw pulse profile without any cuts to the data.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    pulse_pars - parameters corresponding to the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how the pulse profile is 'displaced'.
    no_phase_bins - number of phase bins desired
    mode - whether we want to show or save the plot.

    pulse_pars will have [f,fdot,fdotdot]
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(pulse_pars) != list and type(pulse_pars) != np.ndarray:
        raise TypeError("pulse_pars should either be a list or an array!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,par_list)
    gtis = Lv0_fits2dict.fits2dict(eventfile,2,['START','STOP'])
    T = sum([ (gtis['STOP'][i]-gtis['START'][i]) for i in range(len(gtis['START'])) ])

    times = data_dict['TIME']

    if pulse_pars[1] == 0 and pulse_pars[2] == 0: #i.e., if both \dot{f} and \ddot{f} are zero; that is, if we have no frequency derivatives
        counts = np.ones(len(times))
        shifted_t = times-times[0]
        t_bins = np.linspace(0,np.ceil(shifted_t[-1]),np.ceil(shifted_t[-1])*1/tbin_size+1)
        summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data
        phases, phase_bins, summed_profile = pulse_profile(pulse_pars[0],t_bins[:-1],summed_data,shift,no_phase_bins)
    else:
        phase_bins, summed_profile = pulse_folding(times,T,times[0],pulse_pars[0],pulse_pars[1],pulse_pars[2],no_phase_bins)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if mode != 'overlap':
        plt.figure()
        plt.title('Pulse profile for ' + obj_name + ', ObsID ' + str(obsid),fontsize=12)
#    plt.plot(phase_bins[:-1],summed_profile)
    plt.step(phase_bins[:-1],summed_profile)

    plt.xlabel('Phase', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)
    if mode == 'show':
        plt.show()
    elif mode == 'save':
        filename = 'pp_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        plt.savefig(parent_folder+'/'+filename,dpi=900)
        plt.close()

    return phase_bins[:-1],summed_profile

def partial_t(eventfile,par_list,tbin_size,pulse_pars,shift,no_phase_bins,t1,t2,mode):
    """
    Plot the pulse profile for a desired time interval.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    pulse_pars - parameters corresponding to the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how it is presented.
    no_phase_bins - number of phase bins desired
    t1 - lower time boundary
    t2 - upper time boundary
    mode - whether we want to show or save the plot

    pulse_pars will have [f,fdot,fdotdot]
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(pulse_pars) != list and type(pulse_pars) != np.ndarray:
        raise TypeError("pulse_pars should either be a list or an array!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if t2<t1:
        raise ValueError("t2 should be greater than t1!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,par_list)
    gtis = Lv0_fits2dict.fits2dict(eventfile,2,['START','STOP'])
    T = sum([ (gtis['STOP'][i]-gtis['START'][i]) for i in range(len(gtis['START'])) ])

    if pulse_pars[1] == 0 and pulse_pars[2] == 0:
        truncated_t, truncated_counts = Lv1_data_bin.binning_t(eventfile,par_list,tbin_size,t1,t2)
        phases, phase_bins, summed_profile = pulse_profile(pulse_pars[0],truncated_t[:-1],truncated_counts,shift,no_phase_bins)
    else:
        truncated_t, truncated_counts = Lv1_data_bin.binning_t(eventfile,par_list,tbin_size,t1,t2)
        phase_bins, summed_profile = pulse_folding(truncated_t,T,0,pulse_pars[0],pulse_pars[1],pulse_pars[2],no_phase_bins)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if mode != 'overlap':
        plt.figure()
        plt.title('Pulse profile for ' + obj_name + ', ObsID ' + str(obsid) + '\n Time interval: ' + str(t1) + 's - ' + str(t2) + 's',fontsize=12)
#    plt.plot(phase_bins[:-1], summed_profile)
    plt.step(phase_bins[:-1],summed_profile)
    plt.xlabel('Phase', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        filename = 'pp_' + obsid + '_bin' + str(tbin_size) + 's_' + str(t1) + 's-' + str(t2) + 's.pdf'
        plt.savefig(parent_folder+'/'+filename,dpi=900)
        plt.close()

    return phase_bins[:-1],summed_profile

def partial_E(eventfile,par_list,tbin_size,Ebin_size,pulse_pars,shift,no_phase_bins,E1,E2,mode):
    """
    Plot the pulse profile for a desired energy range.
    [Though I don't think this will be used much. Count/s vs energy is pointless,
    since we're not folding in response matrix information here to get the flux.
    So we're just doing a count/s vs time with an energy cut to the data.]
    INTERJECTION: This caveat is for the spectrum, NOT the pulse profile!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means by in 0.05s
    Ebin_size - the size of the energy bins (in keV!)
    >> e.g., Ebin_size = 0.1 means bin by 0.1keV
    >> e.g., Ebin_size = 0.01 means bin by 0.01keV!
    pulse_pars - parameters corresponding to the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how it is presented.
    no_phase_bins - number of phase bins desired
    E1 - lower energy boundary
    E2 - upper energy boundary

    pulse_pars will have [f,fdot,fdotdot]
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(pulse_pars) != list and type(pulse_pars) != np.ndarray:
        raise TypeError("pulse_pars should either be a list or an array!")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if E2<E1:
        raise ValueError("E2 should be greater than E1!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,par_list)
    gtis = Lv0_fits2dict.fits2dict(eventfile,2,['START','STOP'])
    T = sum([ (gtis['STOP'][i]-gtis['START'][i]) for i in range(len(gtis['START'])) ])

    if pulse_pars[1] == 0 and pulse_pars[2] == 0:
        truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_E(eventfile,par_list,tbin_size,Ebin_size,E1,E2)
        phases, phase_bins, summed_profile = pulse_profile(pulse_pars[0],truncated_t[:-1],truncated_t_counts,shift,no_phase_bins)
    else:
        phase_bins, summed_profile = pulse_folding(truncated_t,T,0,pulse_pars[0],pulse_pars[1],pulse_pars[2],no_phase_bins)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if mode != 'overlap':
        plt.figure()
#    plt.plot(phase_bins[:-1], summed_profile,'-')
#    print(sum(summed_profile)/truncated_t[-1])
    plt.step(phase_bins[:-1],summed_profile)
    plt.xlabel('Phase', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode != 'overlap':
        plt.title('Pulse profile for ' + obj_name + ', ObsID ' + str(obsid)+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        filename = 'pp_' + obsid + '_bin' + str(tbin_size) + 's_' + str(E1) + 'keV-' + str(E2) + 'keV.pdf'
        plt.savefig(parent_folder+'/'+filename,dpi=900)
        plt.close()

    return phase_bins[:-1],summed_profile

def partial_tE(eventfile,par_list,tbin_size,Ebin_size,pulse_pars,shift,no_phase_bins,t1,t2,E1,E2,mode):
    """
    Plot the pulse profile for a desired time interval and desired energy range.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means by in 0.05s
    Ebin_size - the size of the energy bins (in keV!)
    >> e.g., Ebin_size = 0.1 means bin by 0.1keV
    >> e.g., Ebin_size = 0.01 means bin by 0.01keV!
    pulse_pars - parameters corresponding to the pulse
    shift - how much to shift the pulse by in the phase axis.
    It only affects how it is presented.
    no_phase_bins - number of phase bins desired
    t1 - lower time boundary
    t2 - upper time boundary
    E1 - lower energy boundary
    E2 - upper energy boundary
    mode - whether we want to show or save the plot

    pulse_pars will have [f,fdot,fdotdot]
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(pulse_pars) != list and type(pulse_pars) != np.ndarray:
        raise TypeError("pulse_pars should either be a list or an array!")
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

    parent_folder = str(pathlib.Path(eventfile).parent)

    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,par_list)
    gtis = Lv0_fits2dict.fits2dict(eventfile,2,['START','STOP'])
    T = sum([ (gtis['STOP'][i]-gtis['START'][i]) for i in range(len(gtis['START'])) ])

    if pulse_pars[1] == 0 and pulse_pars[2] == 0:
        truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_tE(eventfile,par_list,tbin_size,Ebin_size,t1,t2,E1,E2)
        phases, phase_bins, summed_profile = pulse_profile(pulse_pars[0],truncated_t[:-1],truncated_t_counts,shift,no_phase_bins)
    else:
        phase_bins, summed_profile = pulse_folding(truncated_t,T,0,pulse_pars[0],pulse_pars[1],pulse_pars[2],no_phase_bins)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if mode != 'overlap':
        plt.figure()
        plt.title('Pulse profile for ' + obj_name + ', ObsID ' + str(obsid)+ '\n Time interval: ' + str(t1) + 's - ' + str(t2) + 's'+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)
#    plt.plot(phase_bins[:-1], summed_profile)
    plt.step(phase_bins[:-1],summed_profile)
    plt.xlabel('Phase', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode == 'show':
        plt.show()

    elif mode == 'save':
        filename = 'pp_' + obsid + '_bin' + str(tbin_size) + 's_' + str(t1) + 's-' + str(t2) + 's_' + str(E1) + 'keV-' + str(E2) + 'keV.pdf'
        plt.savefig(parent_folder+'/'+filename,dpi=900)
        plt.close()

    return phase_bins[:-1],summed_profile

################################################################################
### SUBPLOTS

def partial_subplots_E(eventfile,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,subplot_Es,E1,E2,mode):
    """
    Plot the pulse profile for a desired energy range.
    [Though I don't think this will be used much. Count/s vs energy is pointless,
    since we're not folding in response matrix information here to get the flux.
    So we're just doing a count/s vs time with an energy cut to the data.]
    INTERJECTION: This caveat is for the spectrum, NOT the pulse profile!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
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
    subplot_Es - list of tuples defining energy boundaries for pulse profiles
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
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    #should find a way to generalize calling p10,p20,etc in the future..!
    fig,(p10,p20,p30,p40,p50,p60) = plt.subplots(6,1)
    gs = gridspec.GridSpec(6,1)

    for i in range(len(subplot_Es)): #for each tuple of energy boundaries
        truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_E(eventfile,par_list,tbin_size,Ebin_size,subplot_Es[i][0],subplot_Es[i][1])
        phases, phase_bins, summed_profile = pulse_profile(pulse_pars[0],truncated_t[:-1],truncated_t_counts,shift,no_phase_bins)
        plt.subplot(gs[i]).plot(phase_bins[:-1],summed_profile,'-')

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    fig.suptitle(str(obsid),fontsize=12)

    if mode != 'overlap':
        plt.figure()
    plt.xlabel('Phase', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode != 'overlap':
        plt.title('Pulse profile for ' + obj_name + ', ObsID ' + str(obsid)+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        filename = 'pp_subplots_' + obsid + '_bin' + str(tbin_size) + 's_' + str(E1) + 'keV-' + str(E2) + 'keV.pdf'
        plt.savefig(parent_folder+'/'+filename,dpi=900)
        plt.close()

    return phase_bins[:-1],summed_profile

if __name__ == "__main__":
    eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034070101_pipe/ni1034070101_nicersoft_bary.evt'

    #whole(eventfile,['TIME','PI','PI_FAST'],0.01,[0.20801275,0,0],0.4,21,'show')
    #partial_t(eventfile,['TIME','PI','PI_FAST'],1,[0.2081,0,0],0.4,21,0,400,'show')
    #partial_E(eventfile,['TIME','PI','PI_FAST'],1,0.05,[0.2081,0,0],0.4,21,0.3,12,'show')
    #partial_tE(eventfile,['TIME','PI','PI_FAST'],1,0.05,[0.2081,0,0],0.4,21,0,400,0.3,12,'show')
