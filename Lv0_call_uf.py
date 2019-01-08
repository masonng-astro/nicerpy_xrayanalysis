#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 8 11:33am 2019

Opening FITS files and obtaining data from the .uf event files

"""
from __future__ import division,print_function
from astropy.io import fits
import Lv0_dirs

Lv0_dirs.global_par() #obtaining the global parameters

def open_fits(obsid,mpu_no):
    """
    Opening the FITS file for the unfiltered event file

    obsid - Observation ID of the object of interest (10-digit str)
    mpu_no - MPU number, from 0 to 6 inclusive. For the 7 MPUs. str!!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(mpu_no) != str:
        raise TypeError("mpu_no should be a string!")

    event = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_uf/ni' + obsid + '_0mpu' + mpu_no + '_uf.evt'
    event = fits.open(event)
    #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return event

def get_uf(obsid,mpu_no,par_list):
    """
    Getting data from the FITS files, e.g., PI_FAST, times, PI, PI_RATIO, FLAGS

    obsid - Observation ID of the object of interest (10-digit str)
    mpu_no - MPU number, from 0 to 6 inclusive. For the 7 MPUs.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(mpu_no) != str:
        raise TypeError("mpu_no should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = open_fits(obsid,mpu_no)
    data_dict = {}

    for i in range(len(par_list)):
        data_dict[par_list[i]] = event[1].data[par_list[i]]

    return data_dict
