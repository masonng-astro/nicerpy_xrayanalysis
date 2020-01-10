#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thurs Jul 18 10:59am 2020

Generic script to open FITS files for NICER housekeeping data. It is objectively
redundant to have different functions for what is generally the same routine, but
it will help me keep track of things.

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs

Lv0_dirs.global_par()

def get_att(eventfile,par_list):
    """
    Getting data from the .att FITS file! Just provide a path to the event file.
    I can also input the NICERsoft-output event files because this script will
    just search for the ObsID and trawl through the NICER-data folder.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = fits.open(eventfile)
    obsid = event[0].header['OBS_ID']

    att_file = fits.open(Lv0_dirs.NICER_DATADIR+obsid+'/auxil/ni'+obsid+'.att')
    att_dict = {par_list[i]:att_file[1].data[par_list[i]] for i in range(len(par_list))}
    #returns a dictionary in the form {'TIME':[array of TIME values],'QPARAM':[array of QPARAM values]}

    return att_dict

def get_cat(eventfile,par_list):
    """
    Getting data from the .cat FITS file! Just provide a path to the event file.
    I can also input the NICERsoft-output event files because this script will
    just search for the ObsID and trawl through the NICER-data folder.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = fits.open(eventfile)
    obsid = event[0].header['OBS_ID']

    cat_file = fits.open(Lv0_dirs.NICER_DATADIR+obsid+'/auxil/ni'+obsid+'.cat')
    cat_dict = {par_list[i]:cat_file[1].data[par_list[i]] for i in range(len(par_list))}
    #returns a dictionary in the form {'FILENAME':[array of FILENAME values],'FORMAT':[array of FORMAT values]}

    return cat_dict

def get_mkf(eventfile,par_list):
    """
    Getting data from the .mkf FITS file! Just provide a path to the event file.
    I can also input the NICERsoft-output event files because this script will
    just search for the ObsID and trawl through the NICER-data folder.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = fits.open(eventfile)
    obsid = event[0].header['OBS_ID']

    mkf_file = fits.open(Lv0_dirs.NICER_DATADIR+obsid+'/auxil/ni'+obsid+'.mkf')
    mkf_dict = {par_list[i]:mkf_file[1].data[par_list[i]] for i in range(len(par_list))}
    #returns a dictionary in the form {'TIME':[array of TIME values],'FPM_UNDERONLY_COUNT':[array of FPM_UNDERONLY_COUNT values]}

    return mkf_dict

def get_orb(eventfile,par_list):
    """
    Getting data from the .orb FITS file! Just provide a path to the event file.
    I can also input the NICERsoft-output event files because this script will
    just search for the ObsID and trawl through the NICER-data folder.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = fits.open(eventfile)
    obsid = event[0].header['OBS_ID']

    orb_file = fits.open(Lv0_dirs.NICER_DATADIR+obsid+'/auxil/ni'+obsid+'.orb')
    orb_dict = {par_list[i]:orb_file[1].data[par_list[i]] for i in range(len(par_list))}
    #returns a dictionary in the form {'TIME':[array of TIME values],'Vx':[array of Vx values]}

    return orb_dict

def get_hk(eventfile,mpu_no,par_list):
    """
    Getting data from the .hk FITS file! Just provide a path to the event file.
    I can also input the NICERsoft-output event files because this script will
    just search for the ObsID and trawl through the NICER-data folder.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    mpu_no - MPU number, from 0 to 6 inclusive. For the 7 MPUs.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(mpu_no) != str:
        raise TypeError("mpu_no should be a string!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = fits.open(eventfile)
    obsid = event[0].header['OBS_ID']

    hk_file = fits.open(Lv0_dirs.NICER_DATADIR+obsid+'/xti/hk/ni'+obsid+'_0mpu'+mpu_no+'.hk')
    hk_dict = {par_list[i]:hk_file[1].data[par_list[i]] for i in range(len(par_list))}
    #returns a dictionary in the form {'TIME':[array of TIME values],'MPU_UNDER_COUNT':[array of MPU_UNDER_COUNT values]}

    return hk_dict

def get_uf(eventfile,mpu_no,ext,par_list):
    """
    Getting data from the .uf FITS file! Just provide a path to the event file.
    I can also input the NICERsoft-output event files because this script will
    just search for the ObsID and trawl through the NICER-data folder.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    mpu_no - MPU number, from 0 to 6 inclusive. For the 7 MPUs.
    ext - which extension number; 1 for EVENTS, 2 for GTI, 3 for PPS_TREND
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(mpu_no) != str:
        raise TypeError("mpu_no should be a string!")
    if type(ext) != int:
        raise TypeError("ext should be an integer!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = fits.open(eventfile)
    obsid = event[0].header['OBS_ID']

    uf_file = fits.open(Lv0_dirs.NICER_DATADIR+obsid+'/xti/event_uf/ni'+obsid+'_0mpu'+mpu_no+'_uf.evt')
    uf_dict = {par_list[i]:uf_file[ext].data[par_list[i]] for i in range(len(par_list))}
    #returns a dictionary in the form {'TIME':[array of TIME values],'DEADTIME':[array of DEADTIME values]}

    return uf_dict

def get_ufa(eventfile,mpu_no,ext,par_list):
    """
    Getting data from the .ufa FITS file! Just provide a path to the event file.
    I can also input the NICERsoft-output event files because this script will
    just search for the ObsID and trawl through the NICER-data folder.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    mpu_no - MPU number, from 0 to 6 inclusive. For the 7 MPUs.
             MPU number 7 corresponds to the COMBINED file!
    ext - which extension number; 1 for EVENTS, 2 for GTI
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(mpu_no) != str:
        raise TypeError("mpu_no should be a string!")
    if type(ext) != int:
        raise TypeError("ext should be an integer!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = fits.open(eventfile)
    obsid = event[0].header['OBS_ID']

    ufa_file = fits.open(Lv0_dirs.NICER_DATADIR+obsid+'/xti/event_cl/ni'+obsid+'_0mpu'+mpu_no+'_ufa.evt')
    ufa_dict = {par_list[i]:ufa_file[ext].data[par_list[i]] for i in range(len(par_list))}
    #returns a dictionary in the form {'TIME':[array of TIME values],'MPU_UNDER_COUNT':[array of MPU_UNDER_COUNT values]}

    return ufa_dict

if __name__ == "__main__":
    eventfile = '/Volumes/Samsung_T5/NICER-data/1034070101/xti/event_cl/ni1034070101_0mpu7_cl.evt'

    #print(get_att(eventfile,['TIME','QPARAM']))
    #print(get_cat(eventfile,['FILENAME','FORMAT']))
    #print(get_mkf(eventfile,['TIME','NICER_SAA','ANG_DIST']))
    #print(get_orb(eventfile,['TIME','Vx','Vy']))
    #print(get_hk(eventfile,'3',['TIME','GIT_HASH']))
    #print(get_uf(eventfile,'2',1,['TIME','RAWX','DEADTIME']))
    #print(get_ufa(eventfile,'7',1,['TIME','MPU_UNDER_COUNT','PI_RATIO']))

################################################################################
### RELEVANT PARAMETERS:

##### ATT file:
# Extension: ATTITUDE
# TIME, QPARAM, STATE, MODE, SUBMODE_AZ, SUBMODE_EL, ST_VALID, QUATSRC, FINEMEAS
# Extension: INST_ATTITUDE
# TIME, QPARAM, STATE, MODE, SUBMODE_AZ, SUBMODE_EL, ST_VALID, QUATSRC, FINEMEAS

##### CAT file:
# Extension: CATALOG
# FILENAME, FORMAT, TYPE, FILECLAS, DESCRIP, FILESIZE, ARCHSIZE, CHECKSUM, GZIP_CRC, CKSUM_B4

##### MKF file:
# Extension: ORIG_PREFILTER
# TIME, POSITION, VELOCITY, QUATERNION, PNTUNIT, POLAR, RA, DEC, ROLL, SAT_LAT,
# SAT_LON, SAT_ALT, ELV, BR_EARTH, SUNSHINE, FOV_FLAG, SUN_ANGLE, MOON_ANGLE,
# RAM_ANGLE, ANG_DIST, SAA, SAA_TIME, COR_ASCA, COR_SAX, MCILWAIN_L, SUN_RA, SUN_DEC,
# MOON_RA, MOON_DEC, EARTH_RA, EARTH_DEC, TIME_ADJ, ST_BBO, ST_VALID, ST_OBJECTS,
# ST_VIDEO_VDC, ATT_ANG_AZ, ATT_ANG_EL, RA_CMD, DEC_CMD, ATT_ERR_AZ, ATT_ERR_EL,
# ATT_STATE, ATT_MODE, ATT_SUBMODE_AZ, ATT_SUBMODE_EL, TARG_CMD, PPS_SOURCE,
# PPS_ERR_LOWPASS, GPS_INIT, GPS_CONVERGED, NICER_SAA, ST_STARS, ST_FAILCODE,
# MPU_ALL_COUNT, MPU_OVER_COUNT, MPU_UNDER_COUNT, MPU_XRAY_COUNT, TOT_ALL_COUNT,
# TOT_UNDER_COUNT, TOT_OVER_COUNT, TOT_XRAY_COUNT, FPM_ON, NUM_FPM_ON
# Extension: PREFILTER
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

##### ORB file:
# Extension: ORBIT
# TIME, X, Y, Z, Vx, Vy, Vz, GEONS_J2K_TIME_RAW0, GEONS_J2K_WEEK0, ORIG_TIME, QUALITY
# Extension: SPS_ORBIT
# TIME, X, Y, Z, Vx, Vy, Vz, TIME_VALID, QUALITY, SPS_SECS, SPS_WEEK, SDS_GDOP, ORIG_TIME

##### HK file:
# Extension: MPU_HK
#TIME, TICK_LOW32, MPU_P33D_VOLT, MPU_P5D_VOLT, MPU_M5D_VOLT, MPU_P33TEC_VOLT
# MPU_P33M_VOLT, MPU_GNDD_VOLT, MPU_HV_VOLT, MPU_D_TEMP, MPU_GNDA_VOLT,
# MPU_M5A_VOLT, MPU_P5A_VOLT, MPU_P25R_VOLT, MPU_P3R_VOLT, MPU_A_TEMP,
# MPU_PWRBRDG_TEMP, MPU_BAD_CSUM, MPU_INVALID_PKT, MPU_INVALID_ID, MPU_INVALID_MPU,
# MPU_LOWMEM_FIFO, MPU_LOWMEM,SCI, MPU_LOWMEM_OTHER, MPU_ALL_COUNT, MPU_OVER_COUNT,
# MPU_UNDER_COUNT, MPU_XRAY_COUNT, MPU_FPM_TEMP, MPU_PID_TEMP, MPU_HOTSIDE_TEMP,
# MPU_TEC_I, MPU_TEC_VOLT, MPU_BIAS_VOLT, MPU_FAST_LLD, MPU_SLOW_LLD, PACKET_FORMAT
# GIT_HASH, COMPILE_DATE

##### UF file:
# Extension: EVENTS
# TIME, RAWX, RAWY, PHA, PHA_FAST, DET_ID, DEADTIME, EVENT_FLAGS, TICK,
# Extension: GTI
# START, STOP
# Extension: PPS_TREND
# TIME, TICK, PKT_TICK

##### UFA file:
# Extension: EVENTS
# TIME, RAWX, RAWY, PHA, PHA_FAST, DET_ID, DEADTIME, EVENT_FLAGS, TICK, MPU_A_TEMP
# MPU_UNDER_COUNT, PI_FAST, PI, PI_RATIO
# Extension: GTI
# START, STOP
