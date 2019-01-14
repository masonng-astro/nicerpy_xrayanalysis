#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 8 2:11pm 2019

Extracting the GTIs from the FITS files. Use the event_cl files.

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv0_call_eventcl

Lv0_dirs.global_par() #obtaining the global parameters

def get_gtis(obsid,bary,gap):
    """
    Obtaining the GTIs corresponding to the ObsID

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    gap - if time between start of next subregion and the end of the previous
    subregion is greater than the gap, consider that a new subregion!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    event = Lv0_call_eventcl.open_fits(obsid,bary)
    gtis = event[2].data

    data_starts = np.array([gtis[i][0] for i in range(len(gtis))])
    data_stops = np.array([gtis[i][1] for i in range(len(gtis))])

    starts = data_starts - data_starts[0] #so starting at t = 0
    stops = data_stops - data_starts[0]

    index_subreg = [0]
    for i in range(len(starts)-1):
        if starts[i+1]-stops[i] >= 50:
            #if the time between start of next subregion and the end of the
            #previous subregion is > 50s, consider that a new subregion!
            index_subreg.append(int(round(stops[i])))
            index_subreg.append(int(round(starts[i+1])))
    index_subreg.append(int(round(stops[-1])))

    if len(index_subreg)%2 != 0: #if the length of the array is not an even number
        raise ValueError("Note that there are an odd number (so not all are pairs) of values for the subregion definitions!")

    return index_subreg
