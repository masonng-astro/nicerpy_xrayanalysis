#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 8 10:22am 2019

Opening FITS files

"""
from __future__ import division,print_function
from astropy.io import fits
import Lv0_dirs

Lv0_dirs.global_par() #obtaining the global parameters

def open_fits(obsid,bary):
    """
    Opening the FITS file for the cleaned event file

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")

    if bary == True: #if we're using barycentered data
        bary_str = '_bary'
    else: #if not using barycentered data
        bary_str = ''

    event = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl' + bary_str + '.evt'
    event = fits.open(event)
    #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return event

def get_eventcl(obsid,bary,par_list):
    """
    Getting data from the FITS files, e.g., PI_FAST, TIME, PI, PI_RATIO, FLAGS

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = open_fits(obsid,bary)
    data_dict = {}

    for i in range(len(par_list)):
        data_dict[par_list[i]] = event[1].data[par_list[i]]

    return data_dict

################################################################################
if __name__ == "__main__":
    datadict = get_eventcl('1050390140',True,['PI','PI_FAST','TIME'])
    times = datadict['TIME']
    pi = datadict['PI']
    print(len(times)) #1170671 counts?
    print(len(times[(pi>=20)&(pi<=1200)]))
    print(len(pi))


# Variables (TTYPE) from the FITS file headers that I printed

# event_cl_bary

# TIME, RAWX, RAWY, PHA, PHA_FAST, DET_ID, DEADTIME, EVENT_FLAGS, TICK, MPU_A_TEMP,
# MPU_UNDER_COUNT, PI_FAST, PI, PI_RATIO

# gtis

# START, STOP
