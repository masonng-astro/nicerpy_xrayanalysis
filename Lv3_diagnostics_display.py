#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 2:54pm 2019

Saving/showing (aka displaying) the diagnostic plots. This is simply to reduce
the clutter in Lv3_diagnostics.

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs,Lv1_data_bin,Lv2_sources
import Lv0_call_eventcl,Lv0_call_att,Lv0_call_hk
import Lv0_call_mkf,Lv0_call_orb,Lv0_call_uf,Lv0_call_ufa
from scipy import stats
import matplotlib.pyplot as plt

def display_all(obsid,diag_var,lc_t,lc_counts,diag_t,diag_counts,filetype):
    """
    To display the plots for desired time interval. Whether to save or show the
    plots is determined in Lv3_diagnostics.

    obsid - Observation ID of the object of interest (10-digit str)
    diag_var - the diagnostic variable we are looking at
    lc_t - array corresponding to time values for the light curve
    lc_counts - array corresponding to counts for the light curve
    diag_t - array corresponding to times for the diagnostic variable
    diag_counts - array corresponding to counts for the diagnostic variable
    filetype = '.att', '.mkf', '.cl' or ['.hk',mpuno]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(diag_var) != str:
        raise TypeError("diag_var should be a string!")
    if filetype not in ['.att','.mkf','.cl'] and type(filetype) != list and type(filetype) != np.ndarray:
        raise ValueError("filetype should be one of '.att','.mkf','.hk','.eventcl'! Or filetype = ['.hk',mpuno]")
    obj_name = Lv2_sources.obsid_to_obj(obsid)

    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(10,8))
    if filetype == '.att' or filetype == '.mkf' or filetype == '.cl':
        fig.suptitle('Diagnostic plots for ' + obj_name + ', ObsID ' + str(obsid) + '\n Comparing binned light curve and ' + diag_var + ' from ' + filetype + '\n for whole time interval and energy range', fontsize=12)
    elif len(filetype) == 2:
        fig.suptitle('Diagnostic plots for ' + obj_name + ', ObsID ' + str(obsid) + '\n MPU='+filetype[1] + ' - Comparing binned light curve and ' + diag_var + ' from ' + filetype[0] + '\n for whole time interval and energy range', fontsize=12)

    ax1.plot(lc_t[:-1],lc_counts,'b')
    ax1.set_xlabel('Time (s)',fontsize=12)
    ax1.set_ylabel('Counts',fontsize=12)

    ax2.plot(diag_t,diag_counts,'rx-')
    ax2.set_xlabel('Time (s)',fontsize=12)
    ax2.set_ylabel(diag_var)

    plt.subplots_adjust(hspace=0.2)

def display_t(obsid,diag_var,t1,t2,lc_t,lc_counts,diag_t,diag_counts,filetype):
    """
    To display the plots for desired time interval. Whether to save or show the
    plots is determined in Lv3_diagnostics.

    obsid - Observation ID of the object of interest (10-digit str)
    diag_var - the diagnostic variable we are looking at
    t1 - lower time boundary
    t2 - upper time boundary
    lc_t - array corresponding to time values for the light curve
    lc_counts - array corresponding to counts for the light curve
    diag_t - array corresponding to times for the diagnostic variable
    diag_counts - array corresponding to counts for the diagnostic variable
    filetype = '.att', '.mkf', '.cl' or ['.hk',mpuno]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(diag_var) != str:
        raise TypeError("diag_var should be a string!")
    if t2<t1:
        raise ValueError("t2 should be greater than t1!")
    if filetype not in ['.att','.mkf','.cl'] and type(filetype) != list and type(filetype) != np.ndarray:
        raise ValueError("filetype should be one of '.att','.mkf','.hk','.eventcl'! Or filetype = ['.hk',mpuno]")
    obj_name = Lv2_sources.obsid_to_obj(obsid)

    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(10,8))
    if filetype == '.att' or filetype == '.mkf' or filetype == '.cl':
        fig.suptitle('Diagnostic plots for ' + obj_name + ', ObsID ' + str(obsid) + '\n Comparing binned light curve and ' + diag_var + ' from ' + filetype + '\n for time interval: ' + str(t1) +'s-' + str(t2) + 's', fontsize=12)
    elif len(filetype) == 2:
        fig.suptitle('Diagnostic plots for ' + obj_name + ', ObsID ' + str(obsid) + '\n MPU='+filetype[1] + ' - Comparing binned light curve and ' + diag_var + ' from ' + filetype[0] + '\n for time interval: ' + str(t1) +'s-' + str(t2) + 's', fontsize=12)

    ax1.plot(lc_t[:-1],lc_counts,'b')
    ax1.set_xlabel('Time (s)',fontsize=12)
    ax1.set_ylabel('Counts',fontsize=12)

    ax2.plot(diag_t,diag_counts,'rx-')
    ax2.set_xlabel('Time (s)',fontsize=12)
    ax2.set_ylabel(diag_var)

    plt.subplots_adjust(hspace=0.2)

if __name__ == "__main__":
    print('hi') #placeholder, but this is more of a methods script
