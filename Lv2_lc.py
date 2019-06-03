#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs Jan 10 3:02pm 2019

Updated on Mon Jun 3 - Added name_par_list for NICERsoft segments

Plotting light curves

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv0_call_eventcl,Lv0_call_nicersoft_eventcl,Lv1_data_bin,Lv2_sources,Lv2_mkdir
from scipy import stats
import matplotlib.pyplot as plt
import os

Lv0_dirs.global_par() #obtaining the global parameters

def whole(obsid,bary,name_par_list,par_list,tbin_size,mode):
    """
    Plot the entire raw time series without any cuts to the data.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    mode - whether we want to show or save the plot.

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
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    if all(name_par_list[i] == '' for i in range(len(name_par_list))):
        data_dict = Lv0_call_eventcl.get_eventcl(obsid,bary,par_list)
    else:
        data_dict = Lv0_call_nicersoft_eventcl.get_eventcl(obsid,name_par_list,par_list)

    times = data_dict['TIME']
    counts = np.ones(len(times))

    shifted_t = times-times[0]
    t_bins = np.linspace(0,np.ceil(shifted_t[-1]),np.ceil(shifted_t[-1])*1/tbin_size+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

    obj_name = Lv2_sources.obsid_to_obj(obsid)
    plt.plot(t_bins[:-1],summed_data)
    plt.title('Light curve for ' + obj_name + ', ObsID ' + str(obsid),fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)
    if mode == 'show':
        plt.show()
    elif mode == 'save':
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/lc/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/lc/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

#whole('1034090109',True,['TIME','PI','PI_FAST'],1,'show')
#whole('0034070101',True,[True,False,4,100,300,800],['TIME','PI','PI_FAST'],1,'show')

def partial_t(obsid,bary,name_par_list,par_list,tbin_size,t1,t2,mode):
    """
    Plot the time series for a desired time interval.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
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
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    truncated_t, truncated_counts = Lv1_data_bin.binning_t(obsid,bary,name_par_list,par_list,tbin_size,t1,t2)

    obj_name = Lv2_sources.obsid_to_obj(obsid)
    plt.plot(truncated_t[:-1], truncated_counts)
    plt.title('Light curve for ' + obj_name + ', ObsID ' + str(obsid) + '\n Time interval: ' + str(t1) + 's - ' + str(t2) + 's',fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/lc/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/lc/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

#partial_t('0034070101',True,[True,False,1,100,300,800],['PI','TIME','PI_FAST'],1,0,100,'show')

def partial_E(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,E1,E2,mode):
    """
    Plot the time series for a desired energy range.
    [Though I don't think this will be used much. Count/s vs energy is pointless,
    since we're not folding in response matrix information here to get the flux.
    So we're just doing a count/s vs time with an energy cut to the data.]

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
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_E(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,E1,E2)

    obj_name = Lv2_sources.obsid_to_obj(obsid)
    plt.figure()
    plt.plot(truncated_t[:-1], truncated_t_counts)
    plt.title('Light curve for ' + obj_name + ', ObsID ' + str(obsid)+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/lc/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/lc/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

#partial_E('1034090111',True,['TIME','PI','PI_FAST'],1,0.05,8,12,'show')
#partial_E('0034070101',True,[True,False,1,100,300,800],['PI','TIME','PI_FAST'],1,0.05,0.3,12,'show')

def partial_tE(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode):
    """
    Plot the time series for a desired time interval and desired energy range.

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
    t1 - lower time boundary
    t2 - upper time boundary
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
    if t2<t1:
        raise ValueError("t2 should be greater than t1!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_tE(obsid,bary,name_par_list,par_list,tbin_size,Ebin_size,t1,t2,E1,E2)

    obj_name = Lv2_sources.obsid_to_obj(obsid)
    plt.figure()
    plt.plot(truncated_t[:-1], truncated_t_counts)
    plt.title('Light curve for ' + obj_name + ', ObsID ' + str(obsid)+ '\n Time interval: ' + str(t1) + 's - ' + str(t2) + 's'+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Count/' + str(tbin_size) + 's',fontsize=12)

    if mode == 'show':
        plt.show()
    elif mode == 'save':
        if all(name_par_list[i] == '' for i in range(len(name_par_list))):
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/lc/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
                Lv2_mkdir.makedir(dir)
                plt.savefig(filename,dpi=900)
        else:
            dir = Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/outputs/lc/'
            filename = dir + obsid + '_nicersoft_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)
