#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 12:03pm 2019

Returns some of the descriptive information of the ObsID, such as the
date/time of observation start and end, and the peak frequency! Got to manually
put these in. Can also return RA/DEC.

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
from scipy import stats
import os

def obstime(obsid):
    """
    Returns a list showing observation start and end times for a desired ObsID

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    represented = ['0034070101','0034070102','0034070103','0034070104','1034070101',
                   '1034070102','1034070103','1034070104','1034070105','1034070106']

    if obsid not in represented:
        return 'This ObsID does not have an associated start/end time. Either the ObsID is entered incorrectly, or it should be added into the function obstime!'

    obstime_dict = {}
    obstime_dict['0034070101'] = ['2017-06-20T19:06:45','2017-06-20T19:23:53']
    obstime_dict['0034070102'] = ['2017-06-21T02:43:09','2017-06-21T18:33:37']
    obstime_dict['0034070103'] = ['2017-06-22T14:28:56','2017-06-22T18:59:23']
    obstime_dict['0034070104'] = ['2017-06-23T07:34:02','2017-06-23T18:08:23']

    obstime_dict['1034070101'] = ['2017-07-24T23:13:18','2017-07-24T23:21:18']
    obstime_dict['1034070102'] = ['2017-07-25T07:59:34','2017-07-25T17:45:05']
    obstime_dict['1034070103'] = ['2017-07-26T06:00:34','2017-07-26T19:47:25']
    obstime_dict['1034070104'] = ['2017-07-27T01:41:00','2017-07-27T22:01:50']
    obstime_dict['1034070105'] = ['2017-07-28T00:49:40','2017-07-28T10:47:20']
    obstime_dict['1034070106'] = ['2018-11-07T15:23:04','2018-11-07T15:40:03']

    return obstime_dict[obsid]

def ra_dec(obsid):
    """
    Returns the NOMINAL RA,DEC of the corresponding object for a desired ObsID

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    cenx3 = ['0034070101','0034070102','0034070103','0034070104','1034070101',
             '1034070102','1034070103','1034070104','1034070105','1034070106']

    if obsid in cenx3:
        return ['170.3158', '-60.62297']

    if obsid not in cenx3:
        return 'Either you entered the wrong ObsID, or this has not yet been added into the list.'

def peak_freq(obsid):
    """
    Returns the peak frequency (+ error) from the power spectrum over the WHOLE observation.

    Future note: If I want to do it for EACH subsection, then it might be smarter to
    somehow add that into Lv2_ps or Lv2_ps_method?

    obsid - Observation of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    represented = ['0034070101','0034070102','0034070103','0034070104','1034070101',
                   '1034070102','1034070103','1034070104','1034070105','1034070106']

    if obsid not in represented:
        return 'This ObsID does not have an associated peak frequency. Either the ObsID is entered incorrectly, or it should be added into the function peak_freq!'


    peaks = {}
    peaks['0034070101'] = ['0.20846118761251825', '2.501210254644596e-05']
    peaks['0034070102'] = ['0.2080718358508059', '2.516066798883511e-07']
    peaks['0034070103'] = ['0.20854506857953284', '8.951151703709624e-07']
    peaks['0034070104'] = ['0.20796037636355533', '7.456511321335684e-07']
    peaks['1034070101'] = ['0.20799884640430513', '4.969633513422536e-05']
    peaks['1034070102'] = ['Cannot find peak frequency yet.' , ' ']
    peaks['1034070103'] = ['Cannot find peak frequency yet.' , ' ']
    peaks['1034070104'] = ['Cannot find peak frequency yet.' , ' ']
    peaks['1034070105'] = ['Cannot find peak frequency yet.' , ' ']
    peaks['1034070106'] = ['0.20854272550784736', '0.0003410777890108904']

    return peaks[obsid]
