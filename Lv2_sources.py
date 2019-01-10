#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 8 2:11pm 2019

Extracting the GTIs from the FITS files. Use the event_cl files.

"""
from __future__ import division, print_function

def obsid_to_obj(obsid):
    """
    'Convert' the ObsID into the corresponding object name.

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise ValueError("Make sure the ObsID is given as a string!")
    cen_x3 = ['0034070101','0034070102','0034070103','0034070104',
              '1034070101','1034070102','1034070103','1034070104',
              '1034070105','1034070106'] #for Centaurus X-3!
    crab = ['1013010125'] #for Crab Pulsar!
    grs1915p105 = ['1103010117'] #for GRS1915+105

    if obsid in cen_x3:
        return 'Cen X-3'
    elif obsid in crab:
        return 'Crab'
    elif obsid in grs1915p105:
        return 'GRS 1915+105'
    else:
        raise ValueError("Object name not found. Make sure you've uploaded the ObsID + corresponding object name to the module Lv2_sources!")
