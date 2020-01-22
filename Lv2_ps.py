#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Friday Jan 11 10:17am 2019

Plotting power spectra/power spectral densities

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import pathlib
import Lv0_dirs,Lv0_fits2dict,Lv1_data_bin,Lv2_mkdir,Lv2_ps_method,Lv3_detection_level
from scipy import stats, signal
import matplotlib.pyplot as plt
import os

Lv0_dirs.global_par() #obtaining the global parameters

def whole(eventfile,par_list,tbin_size,mode,ps_type,oversampling,xlims,vlines):
    """
    Plot the entire power spectrum without any cuts to the data.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
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
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
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

    if ps_type == 'period':
        plt.figure()
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(t_bins,summed_data,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID: ' + str(obsid) + '\n Periodogram method' + '\n Includes whole time interval and energy range',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_pdgm.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return pdgm_f, pdgm_ps

    if ps_type == 'manual':
        plt.figure()
        manual_f,manual_ps = Lv2_ps_method.manual(t_bins,summed_data,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Manual FFT method' + '\n Includes whole time interval and energy range',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_manual.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return manual_f, manual_ps

    if ps_type == 'both':
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(t_bins,summed_data,xlims,vlines,False,oversampling)
        manual_f,manual_ps = Lv2_ps_method.manual(t_bins,summed_data,xlims,vlines,False,oversampling)
        fig, (ax1,ax2) = plt.subplots(2,1)
        fig.suptitle('Power spectra for ' + obj_name + ', ObsID ' + str(obsid) + '\n both periodogram and manual FFT method' + '\n Includes whole time interval and energy range' , fontsize=12)

        ax1.semilogy(pdgm_f,pdgm_ps,'b-')#/np.mean(pdgm_ps),'b-') #periodogram; arrays already truncated!
        ax1.set_xlabel('Hz',fontsize=12)
        ax1.set_ylabel('Normalized power spectrum',fontsize=10)

        ax2.semilogy(manual_f,manual_ps,'r-')#/np.mean(manual_ps),'r-') #manual FFT; arrays already truncated!
        ax2.set_xlabel('Hz',fontsize=12)
        ax2.set_ylabel('Normalized power spectrum',fontsize=10)

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
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_both.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return pdgm_f, pdgm_ps, manual_f, manual_ps

def partial_t(eventfile,par_list,tbin_size,t1,t2,mode,ps_type,oversampling,xlims,vlines):
    """
    Plot the power spectrum for a desired time interval.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
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
    if ps_type != 'period' and ps_type != 'manual' and ps_type != 'both':
        raise ValueError("ps_type should either be 'period' or 'show' or 'save'!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    truncated_t, truncated_counts = Lv1_data_bin.binning_t(eventfile,par_list,tbin_size,t1,t2)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if ps_type == 'period':
        plt.figure()
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Periodogram method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s' + '\n Whole energy range',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_pdgm_'+str(t1)+'s-'+str(t2)+'s.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return pdgm_f, pdgm_ps

    if ps_type == 'manual':
        plt.figure()
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Manual FFT method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s' + '\n Whole energy range',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_manual_'+str(t1)+'s-'+str(t2)+'s.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return manual_f, manual_ps

    if ps_type == 'both':
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_counts,xlims,vlines,False,oversampling)
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_counts,xlims,vlines,False,oversampling)
        fig, (ax1,ax2) = plt.subplots(2,1)
        fig.suptitle('Power spectra for ' + obj_name + ', ObsID ' + str(obsid) + '\n both periodogram and manual FFT method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s' + '\n Whole energy range',fontsize=12)

        ax1.semilogy(pdgm_f,pdgm_ps,'b-')#/np.mean(pdgm_ps),'b-') #periodogram; arrays already truncated!
        ax1.set_xlabel('Hz',fontsize=12)
        ax1.set_ylabel('Normalized power spectrum',fontsize=10)

        ax2.semilogy(manual_f,manual_ps,'r-')#,'r-')#/np.mean(manual_ps),'r-') #manual FFT; arrays already truncated!
        ax2.set_xlabel('Hz',fontsize=12)
        ax2.set_ylabel('Normalized power spectrum',fontsize=10)

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
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_both_'+str(t1)+'s-'+str(t2)+'s.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return pdgm_f, pdgm_ps, manual_f, manual_ps

def partial_E(eventfile,par_list,tbin_size,Ebin_size,E1,E2,mode,ps_type,oversampling,xlims,vlines):
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
    if ps_type != 'period' and ps_type != 'manual' and ps_type != 'both':
        raise ValueError("ps_type should either be 'period' or 'show' or 'save'!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_E(eventfile,par_list,tbin_size,Ebin_size,E1,E2)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if ps_type == 'period':
        plt.figure()
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_t_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Periodogram method'+ '\n Whole time interval'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_pdgm_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return pdgm_f, pdgm_ps

    if ps_type == 'manual':
        plt.figure()
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_t_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Manual FFT method'+ '\n Whole time interval'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_manual_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return manual_f, manual_ps

    if ps_type == 'both':
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_t_counts,xlims,vlines,False,oversampling)
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_t_counts,xlims,vlines,False,oversampling)
        fig, (ax1,ax2) = plt.subplots(2,1)
        fig.suptitle('Power spectra for ' + obj_name + ', ObsID ' + str(obsid) + '\n both periodogram and manual FFT method'+ '\n Whole time interval'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)

        ax1.semilogy(pdgm_f,pdgm_ps,'b-')#/np.mean(pdgm_ps),'b-') #periodogram; arrays already truncated!
        ax1.set_xlabel('Hz',fontsize=12)
        ax1.set_ylabel('Normalized power spectrum',fontsize=10)

        ax2.semilogy(manual_f,manual_ps,'r-')#/np.mean(manual_ps),'r-') #manual FFT; arrays already truncated!
        ax2.set_xlabel('Hz',fontsize=12)
        ax2.set_ylabel('Normalized power spectrum',fontsize=10)

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
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_both_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return pdgm_f, pdgm_ps, manual_f, manual_ps

def partial_tE(eventfile,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode,ps_type,oversampling,xlims,vlines):
    """
    Plot the time series for a desired time interval and desired energy range.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
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
    if ps_type != 'period' and ps_type != 'manual' and ps_type != 'both':
        raise ValueError("ps_type should either be 'period' or 'show' or 'save'!")
    if type(oversampling) != list and type(oversampling) != np.ndarray:
        raise TypeError("oversampling should either be a list or an array!")
    if type(xlims) != list and type(xlims) != np.ndarray:
        raise TypeError("xlims should either be a list or an array!")
    if type(vlines) != list and type(vlines) != np.ndarray:
        raise TypeError("vlines should either be a list or an array!")

    parent_folder = str(pathlib.Path(eventfile).parent)

    truncated_t, truncated_t_counts, truncated_E, truncated_E_counts = Lv1_data_bin.binning_tE(eventfile,par_list,tbin_size,Ebin_size,t1,t2,E1,E2)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    if ps_type == 'period':
        plt.figure()
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_t_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Periodogram method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_pdgm'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return pdgm_f, pdgm_ps

    if ps_type == 'manual':
        plt.figure()
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_t_counts,xlims,vlines,True,oversampling)
        plt.title('Power spectrum for ' + obj_name + ', ObsID ' + str(obsid) + '\n Manual FFT method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)
        if mode == 'show':
            plt.show()
        elif mode == 'save':
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_manual'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return manual_f, manual_ps

    if ps_type == 'both':
        pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(truncated_t[:-1],truncated_t_counts,xlims,vlines,False,oversampling)
        manual_f,manual_ps = Lv2_ps_method.manual(truncated_t[:-1],truncated_t_counts,xlims,vlines,False,oversampling)
        fig, (ax1,ax2) = plt.subplots(2,1)
        fig.suptitle('Power spectra for ' + obj_name + ', ObsID ' + str(obsid) + '\n both periodogram and manual FFT method'+ '\n Time interval: '+str(t1)+'s-'+str(t2)+'s'+'\n Energy range: '+str(E1)+'keV-'+str(E2)+'keV',fontsize=12)

        ax1.semilogy(pdgm_f,pdgm_ps,'b-')#/np.mean(pdgm_ps),'b-') #periodogram; arrays already truncated!
        ax1.set_xlabel('Hz',fontsize=12)
        ax1.set_ylabel('Normalized power spectrum',fontsize=10)

        ax2.semilogy(manual_f,manual_ps,'r-')#/np.mean(manual_ps),'r-') #manual FFT; arrays already truncated!
        ax2.set_xlabel('Hz',fontsize=12)
        ax2.set_ylabel('Normalized power spectrum',fontsize=10)

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
            filename = 'ps_' + obsid + '_bin' + str(tbin_size) + 's_both'+str(t1)+'s-'+str(t2)+'s_'+str(E1)+'keV-'+str(E2)+'keV.pdf'
            plt.savefig(parent_folder+'/'+filename,dpi=900)
            plt.close()

        return pdgm_f, pdgm_ps, manual_f, manual_ps

if __name__ == "__main__":
    eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034070101_pipe/ni1034070101_nicersoft_bary.evt'
    whole(eventfile,['TIME','PI','PI_FAST'],0.1,'show','both',[True,5],[False,0,1],[True,0.2081])
    #partial_t(eventfile,['TIME','PI','PI_FAST'],0.1,0,100,'show','both',[False,5],[False,0,1],[True,0.2081])
    #partial_E(eventfile,['TIME','PI','PI_FAST'],0.1,0.05,3,8,'show','both',[False,5],[False,0,1],[True,0.2081])
    #partial_tE(eventfile,['TIME','PI','PI_FAST'],0.1,0.05,0,100,3,8,'show','both',[False,5],[False,0,1],[True,0.2081])
