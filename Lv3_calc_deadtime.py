#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 1:42pm 2019

Calculating the deadtime using the corresponding GTIs

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs,Lv0_call_eventcl,Lv0_call_ufa
import Lv1_data_gtis

Lv0_dirs.global_par()

def deadtime(obsid,mpu_no,par_list):
    """
    Calculate the accumulated deadtime for a given observation ID

    obsid - Observation ID of the object of interest (10-digit str)
    mpu_no - Will be '7' for the combined file
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    """
    datadict = Lv0_call_ufa.get_ufa(obsid,mpu_no,par_list) #Use ufa, to include all events that could
    #contribute to the deadtime; both the good events included in the pulse profile and the non-source
    #events filtered out of the profile
    times = datadict['TIME']
    deadtimes = datadict['DEADTIME']
    gtis = Lv1_data_gtis.raw_ufa_gtis(obsid,mpu_no)
    obs_deadtime = 0
    gti_exposure = 0
    count = 0
    for i in range(len(gtis)): #for each GTI interval
        gti_lowerbound = gtis[i][0] #the lower bound of the current GTI interval
        gti_upperbound = gtis[i][1] #the upper bound of the current GTI interval
        deadtime_interval = sum(deadtimes[(times>=gti_lowerbound)&(times<=gti_upperbound)])
        count += len(times[(times>=gti_lowerbound)&(times<=gti_upperbound)])
        obs_deadtime += deadtime_interval
        gti_exposure += gti_upperbound-gti_lowerbound

    return obs_deadtime

def exposure(obsid,bary,par_list):
    """
    Get the on-source, exposure time

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    """
    datadict = Lv0_call_eventcl.get_eventcl(obsid,bary,par_list)
    times = datadict['TIME']
    gtis = Lv1_data_gtis.raw_gtis(obsid,bary)
    gti_exposure = 0
    count = 0
    for i in range(len(gtis)): #for each GTI interval
        gti_lowerbound = gtis[i][0] #the lower bound of the current GTI interval
        gti_upperbound = gtis[i][1] #the upper bound of the current GTI interval
        count += len(times[(times>=gti_lowerbound)&(times<=gti_upperbound)])
        gti_exposure += (gti_upperbound-gti_lowerbound)

    return count,gti_exposure

if __name__ == "__main__":
    print(deadtime('1050390148','7',['TIME','DEADTIME']))
    print(exposure('1050390148',True,['TIME','DEADTIME','PI','PI_FAST']))
