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

def get_gtis(obsid,bary):
    """
    Obtaining the GTIs corresponding to the ObsID

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    event = Lv0_call_eventcl.open_fits(obsid,bary)
    gtis = event[2].data

    data_starts = np.array([gtis[i][0] for i in range(len(gtis))])
    data_stops = np.array([gtis[i][1] for i in range(len(gtis))])

    return data_starts, data_stops
