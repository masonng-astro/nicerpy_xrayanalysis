#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 8 12.02pm 2019

Opening FITS files and obtaining data from the .mkf files

"""
from __future__ import division,print_function
from astropy.io import fits
import Lv0_dirs
import numpy as np

Lv0_dirs.global_par() #obtaining the global parameters

def open_fits(obsid):
    """
    Opening the FITS file for the mkf filter file

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    event = Lv0_dirs.NICER_DATADIR + obsid + '/auxil/ni' + obsid + '.mkf'
    event = fits.open(event)
    #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return event

def get_mkf(obsid,par_list):
    """
    Getting data from the FITS files, e.g., PI_FAST, TIME, PI, PI_RATIO, FLAGS

    obsid - Observation ID of the object of interest (10-digit str)
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = open_fits(obsid)
    data_dict = {}

    for i in range(len(par_list)):
        data_dict[par_list[i]] = event[1].data[par_list[i]]

    return data_dict

################################################################################

# Variables (TTYPE) from the FITS file headers that I printed

# prefilter

# TIME, POSITION, VELOCITY, QUATERNION, PNTUNIT, POLAR, RA, DEC, ROLL, SAT_LAT,
# SAT_LON, SAT_ALT, ELV, BR_EARTH, SUNSHINE, FOV_FLAG, SUN_ANGLE, MOON_ANGLE,
# RAM_ANGLE, ANG_DIST, SAA, SAA_TIME, COR_ASCA, COR_SAX, MCILWAIN_L, SUN_RA, SUN_DEC,
# MOON_RA, MOON_DEC, EARTH_RA, EARTH_DEC, TIME_ADJ, ST_BBO, ST_VALID, ST_OBJECTS,
# ST_VIDEO_VDC, ATT_ANG_AZ, ATT_ANG_EL, RA_CMD, DEC_CMD, ATT_ERR_AZ, ATT_ERR_EL,
# ATT_STATE, ATT_MODE, ATT_SUBMODE_AZ, ATT_SUBMODE_EL, TARG_CMD, PPS_SOURCE,
# PPS_ERR_LOWPASS, GPS_INIT, GPS_CONVERGED, NICER_SAA, ST_STARS, ST_FAILCODE,
# MPU_ALL_COUNT, MPU_OVER_COUNT, MPU_UNDER_COUNT, MPU_XRAY_COUNT, TOT_ALL_COUNT,
# TOT_UNDER_COUNT, TOT_OVER_COUNT, TOT_XRAY_COUNT, FPM_ON, NUM_FPM_ON,
# FPM_RATIO_REJ_COUNT, FPM_XRAY_PI_0000_0025, FPM_XRAY_PI_0035_0200,
# FPM_XRAY_PI_0200_0800, FPM_XRAY_PI_0800_1200, FPM_XRAY_PI_1200_1500,
# FPM_XRAY_PI_1500_1700, FPM_XRAY_PI_COUNT, MPU_DEADTIME, MPU_DOUBLE_COUNT,
# MPU_FT_COUNT, MPU_NOISE25_COUNT, MPU_OVERONLY_COUNT, MPU_UNDERONLY_COUNT,
# FPM_DOUBLE_COUNT, FPM_OVERONLY_COUNT, FPM_UNDERONLY_COUNT, FPM_FT_COUNT,
# FPM_NOISE25_COUNT, MPU_FT_PI_AVG, MPU_FT_PI_ERR, MPU_FT_PI_FAST_AVG,
# MPU_FT_PI_FAST_ERR, MPU_NOISE25_PI_AVG, MPU_NOISE25_PI_ERR, ISS_ATT_STATE,
# ROBO_STATE, VEHICLE_SOYUZ_DC1, VEHICLE_SOYUZ_MRM1, VEHICLE_SOYUZ_MRM2,
# VEHICLE_SOYUZ_SM

# orig_prefilter

# TIME, POSITION, VELOCITY, QUATERNION, PNTUNIT, POLAR, RA, DEC, ROLL, SAT_LAT,
# SAT_LON, SAT_ALT, ELV, BR_EARTH, SUNSHINE, FOV_FLAG, SUN_ANGLE, MOON_ANGLE,
# RAM_ANGLE, ANG_DIST, SAA, SAA_TIME, COR_ASCA, COR_SAX, MCILWAIN_L, SUN_RA, SUN_DEC,
# MOON_RA, MOON_DEC, EARTH_RA, EARTH_DEC, TIME_ADJ, ST_BBO, ST_VALID, ST_OBJECTS,
# ST_VIDEO_VDC, ATT_ANG_AZ, ATT_ANG_EL, RA_CMD, DEC_CMD, ATT_ERR_AZ, ATT_ERR_EL,
# ATT_STATE, ATT_MODE, ATT_SUBMODE_AZ, ATT_SUBMODE_EL, TARG_CMD, PPS_SOURCE,
# PPS_ERR_LOWPASS, GPS_INIT, GPS_CONVERGED, NICER_SAA, ST_STARS, ST_FAILCODE,
# MPU_ALL_COUNT, MPU_OVER_COUNT, MPU_UNDER_COUNT, MPU_XRAY_COUNT, TOT_ALL_COUNT,
# TOT_UNDER_COUNT, TOT_OVER_COUNT, TOT_XRAY_COUNT, FPM_ON, NUM_FPM_ON,
