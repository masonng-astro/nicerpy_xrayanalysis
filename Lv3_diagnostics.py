#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 10:21am 2019

Getting diagnostic plots - so say, how does angular offset change over time
for some desired time interval and/or energy range.

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv0_fits2dict,Lv0_nicer_housekeeping,Lv1_data_bin,Lv2_mkdir
import Lv3_diagnostics_display
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import pathlib
import matplotlib.pyplot as plt

Lv0_dirs.global_par() #obtaining the global parameters

def diag_all(eventfile,par_list,tbin_size,mode,diag_vars):
    """
    Get the diagnostic plots for a desired time interval.
    [Likely too large a range in time (and energy) to be sufficiently useful for
    diagnosis.]

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    mode - whether we want to show or save the plot.
    diag_vars - a dictionary where each key = 'att','mkf','hk', or 'cl', and
    diag_vars[key] provides the list of variables to loop over.
    """
    if type(tbin_size) != int and type(tbin_size) != np.float:
        raise TypeError("tbin_size should be a float or integer!")
    if 'PI' and 'TIME' not in par_list:
        raise ValueError("You should have BOTH 'PI' and 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    parent_folder = str(pathlib.Path(eventfile).parent)
    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    #get the binned light curve
    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,par_list)

    times = data_dict['TIME']
    counts = np.ones(len(times))

    shifted_t = times-times[0]

    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1]*1/tbin_size)+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

    binned_t = t_bins
    binned_counts = summed_data
    #define the variables that we'd like to compare their behavior with the light curve
    att_var = diag_vars['att']
    mkf_var = diag_vars['mkf']
    hk_var = diag_vars['hk']

    ### FOR ATTITUDE
    dict_att = Lv0_nicer_housekeeping.get_att(eventfile,att_var)
    times_att = dict_att['TIME']
    shifted_t = times_att - times_att[0]
    for i in range(1,len(att_var)): #as in, don't compare time with time...
        filtered_att = dict_att[att_var[i]]
        if len(shifted_t) != len(filtered_att):
            raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(att_var[i]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_att)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_all(eventfile,att_var[i],binned_t,binned_counts,shifted_t,filtered_att,'.att')
            plt.show()

    if mode == 'save':
        filename = parent_folder + '/diag_att_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        with PdfPages(filename) as pdf:
            for i in range(1,len(att_var)):
                filtered_att = dict_att[att_var[i]]
                Lv3_diagnostics_display.display_all(eventfile,att_var[i],binned_t,binned_counts,shifted_t,filtered_att,'.att')
                pdf.savefig()
                plt.close()

    ### FOR FILTER
    dict_mkf = Lv0_nicer_housekeeping.get_mkf(eventfile,mkf_var)
    times_mkf = dict_mkf['TIME']
    shifted_t = times_mkf - times_mkf[0]
    for i in range(1,len(mkf_var)): #as in, don't compare time with time...
        filtered_mkf = dict_mkf[mkf_var[i]]
        if len(shifted_t) != len(filtered_mkf):
            raise ValueError("The lengths of arrays shifted t and filtered mkf for variable " + str(mkf_var[i]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_mkf)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_all(eventfile,mkf_var[i],binned_t,binned_counts,shifted_t,filtered_mkf,'.mkf')
            plt.show()

    if mode == 'save':
        filename = parent_folder + '/diag_mkf_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        with PdfPages(filename) as pdf:
            for i in range(1,len(mkf_var)):
                filtered_mkf = dict_mkf[mkf_var[i]]
                Lv3_diagnostics_display.display_all(eventfile,mkf_var[i],binned_t,binned_counts,shifted_t,filtered_mkf,'.mkf')
                pdf.savefig()
                plt.close()

    ### FOR HK
    if mode == 'show':
        for i in range(7):
            dict_hk = Lv0_nicer_housekeeping.get_hk(eventfile,str(i),hk_var)
            times_hk = dict_hk['TIME']
            shifted_t = times_hk - times_hk[0]
            for j in range(1,len(hk_var)): #as in, don't compare time with time...
                filtered_hk = dict_hk[hk_var[j]]
                if len(shifted_t) != len(filtered_hk):
                    raise ValueError("The lengths of arrays shifted t and filtered att for variable " + str(hk_var[j]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_hk)) + ' respectively. This is for HK MPU=' + str(i))
                Lv3_diagnostics_display.display_all(eventfile,hk_var[j],binned_t,binned_counts,shifted_t,filtered_hk,['.hk',str(i)])
                plt.show()

    if mode == 'save':
        filename = parent_folder + '/diag_hk_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        with PdfPages(filename) as pdf:
            for i in range(7):
                dict_hk = Lv0_nicer_housekeeping.get_hk(eventfile,str(i),hk_var)
                times_hk = dict_hk['TIME']
                shifted_t = times_hk - times_hk[0]
                for j in range(1,len(hk_var)): #as in, don't compare time with time...
                    filtered_hk = dict_hk[hk_var[j]]
                    if len(shifted_t) != len(filtered_hk):
                        raise ValueError("The lengths of arrays shifted t and filtered att for variable " + str(hk_var[j]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_hk)) + ' respectively. This is for HK MPU=' + str(i))
                    Lv3_diagnostics_display.display_all(eventfile,hk_var[j],binned_t,binned_counts,shifted_t,filtered_hk,['.hk',str(i)])
                    pdf.savefig()
                    plt.close()

    ### FOR EVENT_CL (BARY)
    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,par_list)
    times_cl = data_dict['TIME']
    shifted_t = times_cl - times_cl[0]

    for i in range(1,len(par_list)): #as in, don't compare time with time...
        filtered_cl = data_dict[par_list[i]]
        if len(shifted_t) != len(filtered_cl):
            raise ValueError("The lengths of arrays shifted t and filtered cl for variable " + str(eventcl_var[i]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_cl)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_all(eventfile,par_list[i],binned_t,binned_counts,shifted_t,filtered_cl,'.cl')
            plt.show()

    if mode == 'save':
        filename = parent_folder + '/diag_cl_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        with PdfPages(filename) as pdf:
            for i in range(1,len(par_list)):
                filtered_cl = data_dict[par_list[i]]
                Lv3_diagnostics_display.display_all(eventfile,eventcl_var[i],binned_t,binned_counts,shifted_t,filtered_cl,'.cl')
                pdf.savefig()
                plt.close()

def diag_t(eventfile,par_list,tbin_size,t1,t2,mode,diag_vars):
    """
    Get the diagnostic plots for a desired time interval.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    t1 - lower time boundary
    t2 - upper time boundary
    mode - whether we want to show or save the plot.
    diag_vars - a dictionary where each key = 'att','mkf','hk', or 'cl', and
    diag_vars[key] provides the list of variables to loop over.
    """
    if type(tbin_size) != int and type(tbin_size) != np.float:
        raise TypeError("tbin_size should be a float or integer!")
    if 'PI' and 'TIME' not in par_list:
        raise ValueError("You should have BOTH 'PI' and 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    parent_folder = str(pathlib.Path(eventfile).parent)
    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    #get the binned light curve
    binned_t, binned_counts = Lv1_data_bin.binning_t(eventfile,par_list,tbin_size,t1,t2)

    #define the variables that we'd like to compare their behavior with the light curve
    att_var = diag_vars['att']
    mkf_var = diag_vars['mkf']
    hk_var = diag_vars['hk']

    ### FOR ATTITUDE
    dict_att = Lv0_nicer_housekeeping.get_att(eventfile,att_var)
    times_att = dict_att['TIME']
    shifted_t_att = times_att - times_att[0]
    filtered_t = shifted_t_att[(shifted_t_att>=t1)&(shifted_t_att<=t2)]
    for i in range(1,len(att_var)): #as in, don't compare time with time...
        filtered_att = dict_att[att_var[i]][(shifted_t_att>=t1)&(shifted_t_att<=t2)]
        if len(filtered_t) != len(filtered_att):
            raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(att_var[i]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_att)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_t(eventfile,att_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_att,'.att')
            plt.show()

    if mode == 'save':
        filename = parent_folder + '/diag_att_' + obsid + '_bin' + str(tbin_size) + 's_' + str(t1) + 's-' + str(t2) + 's.pdf'
        with PdfPages(filename) as pdf:
            for i in range(1,len(att_var)):
                filtered_att = dict_att[att_var[i]][(shifted_t_att>=t1)&(shifted_t_att<=t2)]
                Lv3_diagnostics_display.display_t(eventfile,att_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_att,'.att')
                pdf.savefig()
                plt.close()

    ### FOR FILTER
    dict_mkf = Lv0_nicer_housekeeping.get_mkf(eventfile,mkf_var)
    times_mkf = dict_mkf['TIME']
    shifted_t_mkf = times_mkf - times_mkf[0]
    filtered_t = shifted_t_mkf[(shifted_t_mkf>=t1)&(shifted_t_mkf<=t2)]
    for i in range(1,len(mkf_var)): #as in, don't compare time with time...
        filtered_mkf = dict_mkf[mkf_var[i]][(shifted_t_mkf>=t1)&(shifted_t_mkf<=t2)]
        if len(filtered_t) != len(filtered_mkf):
            raise ValueError("The lengths of arrays filtered t and filtered mkf for variable " + str(mkf_var[i]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_mkf)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_t(eventfile,mkf_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_mkf,'.mkf')
            plt.show()

    if mode == 'save':
        filename = parent_folder + '/diag_mkf_' + obsid + '_bin' + str(tbin_size) + 's_' + str(t1) + 's-' + str(t2) + 's.pdf'
        with PdfPages(filename) as pdf:
            for i in range(1,len(mkf_var)):
                filtered_mkf = dict_mkf[mkf_var[i]][(shifted_t_mkf>=t1)&(shifted_t_mkf<=t2)]
                Lv3_diagnostics_display.display_t(eventfile,mkf_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_mkf,'.mkf')
                pdf.savefig()
                plt.close()

    ### FOR HK
    if mode == 'show':
        for i in range(7):
            dict_hk = Lv0_nicer_housekeeping.get_hk(eventfile,str(i),hk_var)
            times_hk = dict_hk['TIME']
            shifted_t_hk = times_hk - times_hk[0]
            filtered_t = shifted_t_hk[(shifted_t_hk>=t1)&(shifted_t_hk<=t2)]
            for j in range(1,len(hk_var)): #as in, don't compare time with time...
                filtered_hk = dict_hk[hk_var[j]][(shifted_t_hk>=t1)&(shifted_t_hk<=t2)]
                if len(filtered_t) != len(filtered_hk):
                    raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(hk_var[j]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_hk)) + ' respectively. This is for HK MPU=' + str(i))
                Lv3_diagnostics_display.display_t(obsid,hk_var[j],t1,t2,binned_t,binned_counts,filtered_t,filtered_hk,['.hk',str(i)])
                plt.show()

    if mode == 'save':
        filename = parent_folder + '/diag_hk_' + obsid + '_bin' + str(tbin_size) + 's_' + str(t1) + 's-' + str(t2) + 's.pdf'
        with PdfPages(filename) as pdf:
            for i in range(7):
                dict_hk = Lv0_nicer_housekeeping.get_hk(eventfile,str(i),hk_var)
                times_hk = dict_hk['TIME']
                shifted_t_hk = times_hk - times_hk[0]
                filtered_t = shifted_t_hk[(shifted_t_hk>=t1)&(shifted_t_hk<=t2)]
                for j in range(1,len(hk_var)): #as in, don't compare time with time...
                    filtered_hk = dict_hk[hk_var[j]][(shifted_t_hk>=t1)&(shifted_t_hk<=t2)]
                    if len(filtered_t) != len(filtered_hk):
                        raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(hk_var[j]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_hk)) + ' respectively. This is for HK MPU=' + str(i))
                    Lv3_diagnostics_display.display_t(eventfile,hk_var[j],t1,t2,binned_t,binned_counts,filtered_t,filtered_hk,['.hk',str(i)])
                    pdf.savefig()
                    plt.close()

    ### FOR EVENT_CL (BARY)
    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,par_list)
    times_cl = data_dict['TIME']
    shifted_t_cl = times_cl - times_cl[0]
    filtered_t = shifted_t_cl[(shifted_t_cl>=t1)&(shifted_t_cl<=t2)]
    for i in range(1,len(par_list)): #as in, don't compare time with time...
        filtered_cl = data_dict[par_list[i]][(shifted_t_cl>=t1)&(shifted_t_cl<=t2)]
        if len(filtered_t) != len(filtered_cl):
            raise ValueError("The lengths of arrays filtered t and filtered cl for variable " + str(eventcl_var[i]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_cl)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_t(eventfile,par_list[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_cl,'.cl')
            plt.show()

    if mode == 'save':
        filename = parent_folder + '/diag_cl_' + obsid + '_bin' + str(tbin_size) + 's_' + str(t1) + 's-' + str(t2) + 's.pdf'
        with PdfPages(filename) as pdf:
            for i in range(1,len(par_list)):
                filtered_cl = data_dict[par_list[i]][(shifted_t_cl>=t1)&(shifted_t_cl<=t2)]
                Lv3_diagnostics_display.display_t(eventfile,par_file[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_cl,'.cl')
                pdf.savefig()
                plt.close()

if __name__ == "__main__":
    #eventfile = '/Volumes/Samsung_T5/NICER-data/1030180113/'
    #diag_all(eventfile,['TIME','ANG_DIST'],1,'save',{})
    print('hi')
