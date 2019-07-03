#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 9:56am 2019

Plotting phase curves

Updated on Mon Jun 3 - Added name_par_list for NICERsoft segments

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv0_call_eventcl,Lv0_call_nicersoft_eventcl,Lv1_data_bin,Lv2_sources,Lv2_mkdir
from scipy import stats
from PyAstronomy.pyasl import foldAt
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import os

Lv0_dirs.global_par() #obtaining the global parameters

def pulse_profile(f_pulse,times,counts,shift,no_phase_bins):
    """
    Calculating the pulse profile for the observation. Goes from 0 to 2!

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

    phase_bins = np.linspace(0,2,no_phase_bins)
    summed_profile, bin_edges, binnumber = stats.binned_statistic(phases,counts,statistic='sum',bins=phase_bins)

    return phases, phase_bins, summed_profile

def whole(obsid,bary,name_par_list,par_list,tbin_size,f_pulse,shift,no_phase_bins,mode):
    """
    Plot the entire raw pulse profile without any cuts to the data.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
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
    if type(name_par_list) != list and type(name_par_list) != np.ndarray:
        raise TypeError("name_par_list should either be a list or an array!")
    if len(name_par_list) != 6:
        raise ValueError("There seems to be fewer or more values in the list/array than there should be! You should have [GTI_true, E_true, GTIno, segment length, PI1, PI2]")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    if all(name_par_list[i] == '' for i in range(len(name_par_list))):
        data_dict = Lv0_call_eventcl.get_eventcl(obsid,bary,par_list)
    else:
        data_dict = Lv0_call_nicersoft_eventcl.get_eventcl(obsid,name_par_list,par_list)

    times = data_dict['TIME']
    counts = np.ones(len(times))

    shifted_t = times-times[0]
    t_bins = np.linspace(0,np.ceil(shifted_t[-1]),np.ceil(shifted_t[-1])*1/tbin_size+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

    phases, phase_bins, summed_profile = pulse_profile(f_pulse,t_bins[:-1],summed_data,shift,no_phase_bins)

    obj_name = Lv2_sources.obsid_to_obj(obsid)
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
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/pp/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/pp/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()

    return phase_bins[:-1],summed_profile

def partial_t(obsid,bary,name_par_list,par_list,tbin_size,f_pulse,shift,no_phase_bins,t1,t2,mode):
    """
    Plot the pulse profile for a desired time interval.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
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

    name_par_list should be [GTI_true,E_true,GTIno,segment_length,PI1,PI2]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if type(name_par_list) != list and type(name_par_list) != np.ndarray:
        raise TypeError("name_par_list should either be a list or an array!")
    if len(name_par_list) != 6:
        raise ValueError("There seems to be fewer or more values in the list/array than there should be! You should have [GTI_true, E_true, GTIno, segment length, PI1, PI2]")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if t2<t1:
        raise ValueError("t2 should be greater than t1!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    truncated_t, truncated_counts = Lv1_data_bin.binning_t(obsid,bary,name_par_list,par_list,tbin_size,t1,t2)
    phases, phase_bins, summed_profile = pulse_profile(f_pulse,truncated_t[:-1],truncated_counts,shift,no_phase_bins)

    obj_name = Lv2_sources.obsid_to_obj(obsid)
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
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/pp/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/pp/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()

    return phase_bins[:-1],summed_profile

def partial_E(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,E1,E2,mode):
    """
    Plot the pulse profile for a desired energy range.
    [Though I don't think this will be used much. Count/s vs energy is pointless,
    since we're not folding in response matrix information here to get the flux.
    So we're just doing a count/s vs time with an energy cut to the data.]
    INTERJECTION: This caveat is for the spectrum, NOT the pulse profile!

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
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

    name_par_list should be [GTI_true,E_true,GTIno,segment_length,PI1,PI2]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if type(name_par_list) != list and type(name_par_list) != np.ndarray:
        raise TypeError("name_par_list should either be a list or an array!")
    if len(name_par_list) != 6:
        raise ValueError("There seems to be fewer or more values in the list/array than there should be! You should have [GTI_true, E_true, GTIno, segment length, PI1, PI2]")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if E2<E1:
        raise ValueError("E2 should be greater than E1!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_E(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,E1,E2)
    phases, phase_bins, summed_profile = pulse_profile(f_pulse,truncated_t[:-1],truncated_t_counts,shift,no_phase_bins)

    obj_name = Lv2_sources.obsid_to_obj(obsid)
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
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/pp/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/pp/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()

    return phase_bins[:-1],summed_profile

def partial_tE(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,t1,t2,E1,E2,mode):
    """
    Plot the pulse profile for a desired time interval and desired energy range.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
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

    name_par_list should be [GTI_true,E_true,GTIno,segment_length,PI1,PI2]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if type(name_par_list) != list and type(name_par_list) != np.ndarray:
        raise TypeError("name_par_list should either be a list or an array!")
    if len(name_par_list) != 6:
        raise ValueError("There seems to be fewer or more values in the list/array than there should be! You should have [GTI_true, E_true, GTIno, segment length, PI1, PI2]")
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

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_tE(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,t1,t2,E1,E2)
    phases, phase_bins, summed_profile = pulse_profile(f_pulse,truncated_t[:-1],truncated_t_counts,shift,no_phase_bins)

    obj_name = Lv2_sources.obsid_to_obj(obsid)
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
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/pp/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
                Lv2_mkdir.makedir(dir)
                plt.savefig(filename,dpi=900)
                plt.close()
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/pp/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()

    return phase_bins[:-1],summed_profile

################################################################################
### SUBPLOTS

def partial_subplots_E(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,subplot_Es,E1,E2,mode):
    """
    Plot the pulse profile for a desired energy range.
    [Though I don't think this will be used much. Count/s vs energy is pointless,
    since we're not folding in response matrix information here to get the flux.
    So we're just doing a count/s vs time with an energy cut to the data.]
    INTERJECTION: This caveat is for the spectrum, NOT the pulse profile!

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
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

    name_par_list should be [GTI_true,E_true,GTIno,segment_length,PI1,PI2]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if type(name_par_list) != list and type(name_par_list) != np.ndarray:
        raise TypeError("name_par_list should either be a list or an array!")
    if len(name_par_list) != 6:
        raise ValueError("There seems to be fewer or more values in the list/array than there should be! You should have [GTI_true, E_true, GTIno, segment length, PI1, PI2]")
    if 'TIME' not in par_list:
        raise ValueError("You should have 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if E2<E1:
        raise ValueError("E2 should be greater than E1!")
    if mode != 'show' and mode != 'save' and mode != 'overlap':
        raise ValueError("Mode should either be 'show' or 'save' or 'overlap'!")

    #should find a way to generalize calling p10,p20,etc in the future..!
    fig,(p10,p20,p30,p40,p50,p60) = plt.subplots(6,1)
    gs = gridspec.GridSpec(6,1)

    for i in range(len(subplot_Es)): #for each tuple of energy boundaries
        truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_E(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,subplot_Es[i][0],subplot_Es[i][1])
        phases, phase_bins, summed_profile = pulse_profile(f_pulse,truncated_t[:-1],truncated_t_counts,shift,no_phase_bins)
        plt.subplot(gs[i]).plot(phase_bins[:-1],summed_profile,'-')
    fig.suptitle(str(obsid),fontsize=12)
    obj_name = Lv2_sources.obsid_to_obj(obsid)

    if mode != 'overlap':
        plt.figure()
    plt.xlabel('Phase', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode != 'overlap':
        plt.title('Pulse profile for ' + obj_name + ', ObsID ' + str(obsid)+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/pp/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/pp/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
            plt.close()

    return phase_bins[:-1],summed_profile

if __name__ == "__main__":
    whole('0034070101',True,[True,True,1,100,300,800],['TIME','PI','PI_FAST'],1,0.2081,0.4,51,'show')

    partial_t('0034070101',True,[True,True,1,100,300,800],['TIME','PI','PI_FAST'],1,0.2081,0.4,51,0,100,'show')

    partial_E('0034070101',True,[True,True,1,100,300,800],['TIME','PI','PI_FAST'],1,0.05,0.2081,0.4,51,3,8,'show')

    partial_tE('0034070101',True,[True,True,1,100,300,800],['TIME','PI','PI_FAST'],1,0.05,0.2081,0.4,51,0,100,3,8,'show')
