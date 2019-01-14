#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 10:21am 2019

Getting diagnostic plots - so say, how does angular offset change over time
for some desired time interval and/or energy range.

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs,Lv1_data_bin,Lv2_sources,Lv2_mkdir
import Lv0_call_eventcl,Lv0_call_att,Lv0_call_hk
import Lv0_call_mkf,Lv0_call_orb,Lv0_call_uf,Lv0_call_ufa
import Lv3_diagnostics_display
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import matplotlib.pyplot as plt

Lv0_dirs.global_par() #obtaining the global parameters

def diag_all(obsid,bary,par_list,tbin_size,mode,diag_vars):
    """
    Get the diagnostic plots for a desired time interval.
    [Likely too large a range in time (and energy) to be sufficiently useful for
    diagnosis.]

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    tbin_size - the size of the time bins (in seconds!)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    mode - whether we want to show or save the plot.
    diag_vars - a dictionary where each key = 'att','mkf','hk', or 'cl', and
    diag_vars[key] provides the list of variables to loop over.
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(tbin_size) != int and type(tbin_size) != np.float:
        raise TypeError("tbin_size should be a float or integer!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if 'PI' and 'TIME' not in par_list:
        raise ValueError("You should have BOTH 'PI' and 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    #get the binned light curve
    data_dict = Lv0_call_eventcl.get_eventcl(obsid,bary,par_list)

    times = data_dict['TIME']
    counts = np.ones(len(times))

    shifted_t = times-times[0]
    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])*1/tbin_size+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

    binned_t = t_bins
    binned_counts = summed_data
    #define the variables that we'd like to compare their behavior with the light curve
    att_var = diag_vars['att']
    mkf_var = diag_vars['mkf']
    hk_var = diag_vars['hk']
    eventcl_var = diag_vars['cl']

    ### FOR ATTITUDE
    dict_att = Lv0_call_att.get_att(obsid,att_var)
    times_att = dict_att['TIME']
    shifted_t = times_att - times_att[0]
    for i in range(1,len(att_var)): #as in, don't compare time with time...
        filtered_att = dict_att[att_var[i]]
        if len(shifted_t) != len(filtered_att):
            raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(att_var[i]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_att)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_all(obsid,att_var[i],binned_t,binned_counts,shifted_t,filtered_att,'.att')
            plt.show()

    if mode == 'save':
        dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/diagnostics/'
        if bary == True:
            filename = dir + 'att_' + obsid + '_bary_bin' + str(tbin_size) + 's.pdf'
        elif bary == False:
            filename = dir + 'att_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        Lv2_mkdir.makedir(dir)
        with PdfPages(filename) as pdf:
            for i in range(1,len(att_var)):
                filtered_att = dict_att[att_var[i]]
                Lv3_diagnostics_display.display_all(obsid,att_var[i],binned_t,binned_counts,shifted_t,filtered_att,'.att')
                pdf.savefig()
                plt.close()

    ### FOR FILTER
    dict_mkf = Lv0_call_mkf.get_mkf(obsid,mkf_var)
    times_mkf = dict_mkf['TIME']
    shifted_t = times_mkf - times_mkf[0]
    for i in range(1,len(mkf_var)): #as in, don't compare time with time...
        filtered_mkf = dict_mkf[mkf_var[i]]
        if len(shifted_t) != len(filtered_mkf):
            raise ValueError("The lengths of arrays filtered t and filtered mkf for variable " + str(mkf_var[i]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_mkf)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_all(obsid,mkf_var[i],binned_t,binned_counts,shifted_t,filtered_mkf,'.mkf')
            plt.show()

    if mode == 'save':
        dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/diagnostics/'
        if bary == True:
            filename = dir + 'mkf_' + obsid + '_bary_bin' + str(tbin_size) + 's.pdf'
        elif bary == False:
            filename = dir + 'mkf_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        Lv2_mkdir.makedir(dir)
        with PdfPages(filename) as pdf:
            for i in range(1,len(mkf_var)):
                filtered_mkf = dict_mkf[mkf_var[i]]
                Lv3_diagnostics_display.display_all(obsid,mkf_var[i],binned_t,binned_counts,shifted_t,filtered_mkf,'.mkf')
                pdf.savefig()
                plt.close()

    ### FOR HK
    if mode == 'show':
        for i in range(7):
            dict_hk = Lv0_call_hk.get_hk(obsid,str(i),hk_var)
            times_hk = dict_hk['TIME']
            shifted_t = times_hk - times_hk[0]
            for j in range(1,len(hk_var)): #as in, don't compare time with time...
                filtered_hk = dict_hk[hk_var[j]]
                if len(shifted_t) != len(filtered_hk):
                    raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(hk_var[j]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_hk)) + ' respectively. This is for HK MPU=' + str(i))
                Lv3_diagnostics_display.display_all(obsid,hk_var[j],binned_t,binned_counts,shifted_t,filtered_hk,['.hk',str(i)])
                plt.show()

    if mode == 'save':
        dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/diagnostics/'
        if bary == True:
            filename = dir + 'hk_' + obsid + '_bary_bin' + str(tbin_size) + 's.pdf'
        elif bary == False:
            filename = dir + 'hk_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        Lv2_mkdir.makedir(dir)
        with PdfPages(filename) as pdf:
            for i in range(7):
                dict_hk = Lv0_call_hk.get_hk(obsid,str(i),hk_var)
                times_hk = dict_hk['TIME']
                shifted_t = times_hk - times_hk[0]
                for j in range(1,len(hk_var)): #as in, don't compare time with time...
                    filtered_hk = dict_hk
                    if len(shifted_t) != len(filtered_hk):
                        raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(hk_var[j]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_hk)) + ' respectively. This is for HK MPU=' + str(i))
                    Lv3_diagnostics_display.display_all(obsid,hk_var[j],binned_t,binned_counts,shifted_t,filtered_hk,['.hk',str(i)])
                    pdf.savefig()
                    plt.close()

    ### FOR EVENT_CL (BARY)
    dict_eventcl = Lv0_call_eventcl.get_eventcl(obsid,bary,eventcl_var)
    times_cl = dict_eventcl['TIME']
    shifted_t = times_cl - times_cl[0]

    for i in range(1,len(eventcl_var)): #as in, don't compare time with time...
        filtered_cl = dict_eventcl[eventcl_var[i]]
        if len(shifted_t) != len(filtered_cl):
            raise ValueError("The lengths of arrays filtered t and filtered cl for variable " + str(eventcl_var[i]) + ' are different, with ' + str(len(shifted_t)) + ' and ' + str(len(filtered_cl)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_all(obsid,eventcl_var[i],binned_t,binned_counts,shifted_t,filtered_cl,'.cl')
            plt.show()

    if mode == 'save':
        dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/diagnostics/'
        if bary == True:
            filename = dir + 'cl_' + obsid + '_bary_bin' + str(tbin_size) + 's.pdf'
        elif bary == False:
            filename = dir + 'cl_' + obsid + '_bin' + str(tbin_size) + 's.pdf'
        Lv2_mkdir.makedir(dir)
        with PdfPages(filename) as pdf:
            for i in range(1,len(eventcl_var)):
                filtered_cl = dict_eventcl[eventcl_var[i]]
                Lv3_diagnostics_display.display_all(obsid,eventcl_var[i],binned_t,binned_counts,shifted_t,filtered_cl,'.cl')
                pdf.savefig()
                plt.close()

def diag_t(obsid,bary,par_list,tbin_size,t1,t2,mode,diag_vars):
    """
    Get the diagnostic plots for a desired time interval.

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
    diag_vars - a dictionary where each key = 'att','mkf','hk', or 'cl', and
    diag_vars[key] provides the list of variables to loop over.
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(tbin_size) != int and type(tbin_size) != np.float:
        raise TypeError("tbin_size should be a float or integer!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if 'PI' and 'TIME' not in par_list:
        raise ValueError("You should have BOTH 'PI' and 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")

    #get the binned light curve
    binned_t, binned_counts = Lv1_data_bin.binning_t(obsid,bary,par_list,tbin_size,t1,t2)

    #define the variables that we'd like to compare their behavior with the light curve
    att_var = diag_vars['att']
    mkf_var = diag_vars['mkf']
    hk_var = diag_vars['hk']
    eventcl_var = diag_vars['cl']

    ### FOR ATTITUDE
    dict_att = Lv0_call_att.get_att(obsid,att_var)
    times_att = dict_att['TIME']
    shifted_t_att = times_att - times_att[0]
    filtered_t = shifted_t_att[(shifted_t_att>=t1)&(shifted_t_att<=t2)]
    for i in range(1,len(att_var)): #as in, don't compare time with time...
        filtered_att = dict_att[att_var[i]][(shifted_t_att>=t1)&(shifted_t_att<=t2)]
        if len(filtered_t) != len(filtered_att):
            raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(att_var[i]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_att)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_t(obsid,att_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_att,'.att')
            plt.show()

    if mode == 'save':
        dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/diagnostics/'
        if bary == True:
            filename = dir + 'att_' + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
        elif bary == False:
            filename = dir + 'att_' + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
        Lv2_mkdir.makedir(dir)
        with PdfPages(filename) as pdf:
            for i in range(1,len(att_var)):
                filtered_att = dict_att[att_var[i]][(shifted_t_att>=t1)&(shifted_t_att<=t2)]
                Lv3_diagnostics_display.display_t(obsid,att_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_att,'.att')
                pdf.savefig()
                plt.close()

    ### FOR FILTER
    dict_mkf = Lv0_call_mkf.get_mkf(obsid,mkf_var)
    times_mkf = dict_mkf['TIME']
    shifted_t_mkf = times_mkf - times_mkf[0]
    filtered_t = shifted_t_mkf[(shifted_t_mkf>=t1)&(shifted_t_mkf<=t2)]
    for i in range(1,len(mkf_var)): #as in, don't compare time with time...
        filtered_mkf = dict_mkf[mkf_var[i]][(shifted_t_mkf>=t1)&(shifted_t_mkf<=t2)]
        if len(filtered_t) != len(filtered_mkf):
            raise ValueError("The lengths of arrays filtered t and filtered mkf for variable " + str(mkf_var[i]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_mkf)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_t(obsid,mkf_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_mkf,'.mkf')
            plt.show()

    if mode == 'save':
        dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/diagnostics/'
        if bary == True:
            filename = dir + 'mkf_' + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
        elif bary == False:
            filename = dir + 'mkf_' + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
        Lv2_mkdir.makedir(dir)
        with PdfPages(filename) as pdf:
            for i in range(1,len(mkf_var)):
                filtered_mkf = dict_mkf[mkf_var[i]][(shifted_t_mkf>=t1)&(shifted_t_mkf<=t2)]
                Lv3_diagnostics_display.display_t(obsid,mkf_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_mkf,'.mkf')
                pdf.savefig()
                plt.close()

    ### FOR HK
    if mode == 'show':
        for i in range(7):
            dict_hk = Lv0_call_hk.get_hk(obsid,str(i),hk_var)
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
        dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/diagnostics/'
        if bary == True:
            filename = dir + 'hk_' + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
        elif bary == False:
            filename = dir + 'hk_' + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
        Lv2_mkdir.makedir(dir)
        with PdfPages(filename) as pdf:
            for i in range(7):
                dict_hk = Lv0_call_hk.get_hk(obsid,str(i),hk_var)
                times_hk = dict_hk['TIME']
                shifted_t_hk = times_hk - times_hk[0]
                filtered_t = shifted_t_hk[(shifted_t_hk>=t1)&(shifted_t_hk<=t2)]
                for j in range(1,len(hk_var)): #as in, don't compare time with time...
                    filtered_hk = dict_hk[hk_var[j]][(shifted_t_hk>=t1)&(shifted_t_hk<=t2)]
                    if len(filtered_t) != len(filtered_hk):
                        raise ValueError("The lengths of arrays filtered t and filtered att for variable " + str(hk_var[j]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_hk)) + ' respectively. This is for HK MPU=' + str(i))
                    Lv3_diagnostics_display.display_t(obsid,hk_var[j],t1,t2,binned_t,binned_counts,filtered_t,filtered_hk,['.hk',str(i)])
                    pdf.savefig()
                    plt.close()

    ### FOR EVENT_CL (BARY)
    dict_eventcl = Lv0_call_eventcl.get_eventcl(obsid,bary,eventcl_var)
    times_cl = dict_eventcl['TIME']
    shifted_t_cl = times_cl - times_cl[0]
    filtered_t = shifted_t_cl[(shifted_t_cl>=t1)&(shifted_t_cl<=t2)]
    for i in range(1,len(eventcl_var)): #as in, don't compare time with time...
        filtered_cl = dict_eventcl[eventcl_var[i]][(shifted_t_cl>=t1)&(shifted_t_cl<=t2)]
        if len(filtered_t) != len(filtered_cl):
            raise ValueError("The lengths of arrays filtered t and filtered cl for variable " + str(eventcl_var[i]) + ' are different, with ' + str(len(filtered_t)) + ' and ' + str(len(filtered_cl)) + ' respectively.')

        if mode == 'show':
            Lv3_diagnostics_display.display_t(obsid,eventcl_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_cl,'.cl')
            plt.show()

    if mode == 'save':
        dir = Lv0_dirs.BASE_DIR+'outputs/' + obsid + '/diagnostics/'
        if bary == True:
            filename = dir + 'cl_' + obsid + '_bary_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
        elif bary == False:
            filename = dir + 'cl_' + obsid + '_bin' + str(tbin_size) + 's_'+str(t1)+'s-'+str(t2)+'s.pdf'
        Lv2_mkdir.makedir(dir)
        with PdfPages(filename) as pdf:
            for i in range(1,len(eventcl_var)):
                filtered_cl = dict_eventcl[eventcl_var[i]][(shifted_t_cl>=t1)&(shifted_t_cl<=t2)]
                Lv3_diagnostics_display.display_t(obsid,eventcl_var[i],t1,t2,binned_t,binned_counts,filtered_t,filtered_cl,'.cl')
                pdf.savefig()
                plt.close()
