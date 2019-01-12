#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Friday Jan 11 10:17am 2019

Plotting power spectra/power spectral densities

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv0_call_eventcl,Lv1_data_bin,Lv2_sources,Lv2_mkdir,Lv2_ps_method
from scipy import stats, signal
import matplotlib.pyplot as plt
import os

Lv0_dirs.global_par() #obtaining the global parameters

def whole(obsid,bary,par_list,tbin_size,mode,ps_type,oversampling,xlims,vlines):
    """
    Plot the entire power spectrum without any cuts to the data.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    mode - whether we want to show or save the plot.
    ps_type - obtain power spectrum through the periodogram method ('period') or
    the manual FFT way ('manual') or both ('both')
    oversampling - whether to perform oversampling. Array will consist of
    [True/False, oversampling factor]
    xlims - a list or array: first entry = True/False as to whether to impose an
    xlim; second and third entry correspond to the desired x-limits of the plot
    vlines - a list or array: first entry = True/False as to whether to draw
    a vertical line in the plot; second entry is the equation for the vertical line
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
    if ps_type != 'period' and ps_type != 'manual' and ps_type != 'both':
        raise ValueError("ps_type should either be 'period' or 'show' or 'save'!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")

    data_dict = Lv0_call_eventcl.get_eventcl(obsid,bary,par_list)

    times = data_dict['TIME']
    counts = np.ones(len(times))

    shifted_t = times-times[0]
    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])*1/tbin_size+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

    obj_name = Lv2_sources.obsid_to_obj(obsid)

    if ps_type == 'period':
        plt.figure()
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(t_bins,summed_data,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Periodogram method' + '\n Includes whole time interval and energy range',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_pdgm.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_pdgm.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

    if ps_type == 'manual':
        plt.figure()
        manual_f,manual_ps = Lv2_ps_method.manual(t_bins,summed_data,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Manual FFT method' + '\n Includes whole time interval and energy range',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_manual.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_manual.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

    if ps_type == 'both':
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(t_bins,summed_data,xlims,vlines,False,oversampling)
        manual_f,manual_ps = Lv2_ps_method.manual(t_bins,summed_data,xlims,vlines,False,oversampling)
        fig, (ax1,ax2) = plt.subplots(2,1)
        fig.suptitle('Power spectra for ' + obj_name + ', ObsID ' + str(obsid) + '\n both periodogram and manual FFT method' + '\n Includes whole time interval and energy range' , fontsize=12)

        ax1.semilogy(pdgm_f,pdgm_ps/np.mean(pdgm_ps),'b-') #periodogram; arrays already truncated!
        ax1.set_xlabel('Hz',fontsize=12)
        ax1.set_ylabel('Normalized power spectrum',fontsize=12)

        ax2.semilogy(manual_f,manual_ps/np.mean(manual_ps),'r-') #manual FFT; arrays already truncated!
        ax2.set_xlabel('Hz',fontsize=12)
        ax2.set_ylabel('Normalized power spectrum',fontsize=12)

        if xlims[0] == True:
            ax1.set_xlim([xlims[1],xlims[2]])
            ax2.set_xlim([xlims[1],xlims[2]])
        if vlines[0] == True:
            ax1.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
            ax2.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
            ax2.axhline(y=2,color='k',alpha=0.3,lw=0.3)

        plt.subplots_adjust(hspace=0.2)

        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_both.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_both.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

def partial_t(obsid,bary,par_list,tbin_size,t1,t2,mode,ps_type,oversampling,xlims,vlines):
    """
    Plot the power spectrum for a desired time interval.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    t1 - lower time boundary
    t2 - upper time boundary
    mode - whether we want to show or save the plot.
    ps_type - obtain power spectrum through the periodogram method ('period') or
    the manual FFT way ('manual') or both ('both')
    oversampling - whether to perform oversampling. Array will consist of
    [True/False, oversampling factor]
    xlims - a list or array: first entry = True/False as to whether to impose an
    xlim; second and third entry correspond to the desired x-limits of the plot
    vlines - a list or array: first entry = True/False as to whether to draw
    a vertical line in the plot; second entry is the equation for the vertical line
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
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")
    if ps_type != 'period' and ps_type != 'manual' and ps_type != 'both':
        raise ValueError("ps_type should either be 'period' or 'show' or 'save'!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")

    truncated_t, truncated_counts = Lv1_data_bin.binning_t(obsid,bary,par_list,tbin_size,t1,t2)

    obj_name = Lv2_sources.obsid_to_obj(obsid)

    if ps_type == 'period':
        plt.figure()
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Periodogram method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s' + '\n Whole energy range',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_pdgm_'+str(t1)+'s-'+str(t2)+'s.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_pdgm_'+str(t1)+'s-'+str(t2)+'s.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

    if ps_type == 'manual':
        plt.figure()
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Manual FFT method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s' + '\n Whole energy range',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_manual_'+str(t1)+'s-'+str(t2)+'s.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_manual_'+str(t1)+'s-'+str(t2)+'s.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

    if ps_type == 'both':
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_counts,xlims,vlines,False,oversampling)
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_counts,xlims,vlines,False,oversampling)
        fig, (ax1,ax2) = plt.subplots(2,1)
        fig.suptitle('Power spectra for ' + obj_name + ', ObsID ' + str(obsid) + '\n both periodogram and manual FFT method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s' + '\n Whole energy range',fontsize=12)

        ax1.semilogy(pdgm_f,pdgm_ps/np.mean(pdgm_ps),'b-') #periodogram; arrays already truncated!
        ax1.set_xlabel('Hz',fontsize=12)
        ax1.set_ylabel('Normalized power spectrum',fontsize=12)

        ax2.semilogy(manual_f,manual_ps/np.mean(manual_ps),'r-') #manual FFT; arrays already truncated!
        ax2.set_xlabel('Hz',fontsize=12)
        ax2.set_ylabel('Normalized power spectrum',fontsize=12)

        if xlims[0] == True:
            ax1.set_xlim([xlims[1],xlims[2]])
            ax2.set_xlim([xlims[1],xlims[2]])
        if vlines[0] == True:
            ax1.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
            ax2.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
            ax2.axhline(y=2,color='k',alpha=0.3,lw=0.3)

        plt.subplots_adjust(hspace=0.2)

        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_both_'+str(t1)+'s-'+str(t2)+'s.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_both_'+str(t1)+'s-'+str(t2)+'s.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

def partial_E(obsid,bary,par_list,tbin_size,Ebin_size,E1,E2,mode,ps_type,oversampling,xlims,vlines):
    """
    Plot the time series for a desired energy range.
    [Though I don't think this will be used much. Count/s vs energy is pointless,
    since we're not folding in response matrix information here to get the flux.
    So we're just doing a count/s vs time with an energy cut to the data.]

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
    E1 - lower energy boundary
    E2 - upper energy boundary
    mode - whether we want to show or save the plot.
    ps_type - obtain power spectrum through the periodogram method ('period') or
    the manual FFT way ('manual') or both ('both')
    oversampling - whether to perform oversampling. Array will consist of
    [True/False, oversampling factor]
    xlims - a list or array: first entry = True/False as to whether to impose an
    xlim; second and third entry correspond to the desired x-limits of the plot
    vlines - a list or array: first entry = True/False as to whether to draw
    a vertical line in the plot; second entry is the equation for the vertical line
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
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")
    if ps_type != 'period' and ps_type != 'manual' and ps_type != 'both':
        raise ValueError("ps_type should either be 'period' or 'show' or 'save'!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_E(obsid,bary,par_list,tbin_size,Ebin_size,E1,E2)

    obj_name = Lv2_sources.obsid_to_obj(obsid)

    if ps_type == 'period':
        plt.figure()
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_t_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Periodogram method'+ '\n Whole time interval'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_pdgm_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_pdgm_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

    if ps_type == 'manual':
        plt.figure()
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_t_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Manual FFT method'+ '\n Whole time interval'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_manual_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_manual_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

    if ps_type == 'both':
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_t_counts,xlims,vlines,False,oversampling)
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_t_counts,xlims,vlines,False,oversampling)
        fig, (ax1,ax2) = plt.subplots(2,1)
        fig.suptitle('Power spectra for ' + obj_name + ', ObsID ' + str(obsid) + '\n both periodogram and manual FFT method'+ '\n Whole time interval'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)

        ax1.semilogy(pdgm_f,pdgm_ps/np.mean(pdgm_ps),'b-') #periodogram; arrays already truncated!
        ax1.set_xlabel('Hz',fontsize=12)
        ax1.set_ylabel('Normalized power spectrum',fontsize=12)

        ax2.semilogy(manual_f,manual_ps/np.mean(manual_ps),'r-') #manual FFT; arrays already truncated!
        ax2.set_xlabel('Hz',fontsize=12)
        ax2.set_ylabel('Normalized power spectrum',fontsize=12)

        if xlims[0] == True:
            ax1.set_xlim([xlims[1],xlims[2]])
            ax2.set_xlim([xlims[1],xlims[2]])
        if vlines[0] == True:
            ax1.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
            ax2.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
            ax2.axhline(y=2,color='k',alpha=0.3,lw=0.3)

        plt.subplots_adjust(hspace=0.2)

        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_both_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_both_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)


def partial_tE(obsid,bary,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode,ps_type,oversampling,xlims,vlines):
    """
    Plot the time series for a desired time interval and desired energy range.

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
    t1 - lower time boundary
    t2 - upper time boundary
    E1 - lower energy boundary
    E2 - upper energy boundary
    mode - whether we want to show or save the plot.
    ps_type - obtain power spectrum through the periodogram method ('period') or
    the manual FFT way ('manual') or both ('both')
    oversampling - whether to perform oversampling. Array will consist of
    [True/False, oversampling factor]
    xlims - a list or array: first entry = True/False as to whether to impose an
    xlim; second and third entry correspond to the desired x-limits of the plot
    vlines - a list or array: first entry = True/False as to whether to draw
    a vertical line in the plot; second entry is the equation for the vertical line
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
    if ps_type != 'period' and ps_type != 'manual' and ps_type != 'both':
        raise ValueError("ps_type should either be 'period' or 'show' or 'save'!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_tE(obsid,bary,par_list,tbin_size,Ebin_size,t1,t2,E1,E2)

    obj_name = Lv2_sources.obsid_to_obj(obsid)

    if ps_type == 'period':
        plt.figure()
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_t_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Periodogram method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

    if ps_type == 'manual':
        plt.figure()
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_t_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Manual FFT method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

    if ps_type == 'both':
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_t_counts,xlims,vlines,False,oversampling)
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_t_counts,xlims,vlines,False,oversampling)
        fig, (ax1,ax2) = plt.subplots(2,1)
        fig.suptitle('Power spectra for ' + obj_name + ', ObsID ' + str(obsid) + '\n both periodogram and manual FFT method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)

        ax1.semilogy(pdgm_f,pdgm_ps/np.mean(pdgm_ps),'b-') #periodogram; arrays already truncated!
        ax1.set_xlabel('Hz',fontsize=12)
        ax1.set_ylabel('Normalized power spectrum',fontsize=12)

        ax2.semilogy(manual_f,manual_ps/np.mean(manual_ps),'r-') #manual FFT; arrays already truncated!
        ax2.set_xlabel('Hz',fontsize=12)
        ax2.set_ylabel('Normalized power spectrum',fontsize=12)

        if xlims[0] == True:
            ax1.set_xlim([xlims[1],xlims[2]])
            ax2.set_xlim([xlims[1],xlims[2]])
        if vlines[0] == True:
            ax1.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
            ax2.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
            ax2.axhline(y=2,color='k',alpha=0.3,lw=0.3)

        plt.subplots_adjust(hspace=0.2)

        if mode == 'show':
            plt.show()
        elif mode == 'save':
            dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/ps/'
            if bary == True:
                filename = dir + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            elif bary == False:
                filename = dir + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            Lv2_mkdir.makedir(dir)
            plt.savefig(filename,dpi=900)

partial_tE('1034070104',True,['TIME','PI','PI_FAST'],0.1,0.05,11113,11945,0.3,2.7,'save','both',[True,5],[True,0,1],[True,0.2078])
partial_tE('1034070104',True,['TIME','PI','PI_FAST'],0.1,0.05,11113,11945,2.7,12,'save','both',[True,5],[True,0,1],[True,0.2078])
