#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 8 11:32am 2019

Opening FITS files and obtaining data from the .hk files

"""
from __future__ import division,print_function
from astropy.io import fits
import Lv0_dirs
import numpy as np

Lv0_dirs.global_par() #obtaining the global parameters

def open_fits(obsid,mpu_no):
    """
    Opening the FITS file for the housekeeping file

    obsid - Observation ID of the object of interest (10-digit str)
    mpu_no - MPU number, from 0 to 6 inclusive. For the 7 MPUs. str!!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(mpu_no) != str:
        raise TypeError("mpu_no should be a string!")

    event = Lv0_dirs.NICER_DATADIR + obsid + '/xti/hk/ni' + obsid + '_0mpu' + mpu_no + '.hk'
    event = fits.open(event)
    #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return event

def get_hk(obsid,mpu_no,par_list):
    """
    Getting data from the FITS files, e.g., PI_FAST, TIME, PI, PI_RATIO, FLAGS

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

################################################################################

# Variables (TTYPE) from the FITS file headers that I printed

# TIME, TICK_LOW32, MPU_P33D_VOLT, MPU_P5D_VOLT, MPU_M5D_VOLT, MPU_P33TEC_VOLT
# MPU_P33M_VOLT, MPU_GNDD_VOLT, MPU_HV_VOLT, MPU_D_TEMP, MPU_GNDA_VOLT,
# MPU_M5A_VOLT, MPU_P5A_VOLT, MPU_P25R_VOLT, MPU_P3R_VOLT, MPU_A_TEMP,
# MPU_PWRBRDG_TEMP, MPU_BAD_CSUM, MPU_INVALID_PKT, MPU_INVALID_ID, MPU_INVALID_MPU,
# MPU_LOWMEM_FIFO, MPU_LOWMEM,SCI, MPU_LOWMEM_OTHER, MPU_ALL_COUNT, MPU_OVER_COUNT,
# MPU_UNDER_COUNT, MPU_XRAY_COUNT, MPU_FPM_TEMP, MPU_PID_TEMP, MPU_HOTSIDE_TEMP,
# MPU_TEC_I, MPU_TEC_VOLT, MPU_BIAS_VOLT, MPU_FAST_LLD, MPU_SLOW_LLD, PACKET_FORMAT
# GIT_HASH, COMPILE_DATE
