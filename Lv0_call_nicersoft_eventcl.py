#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 3 9:12am 2019

Opening cleaned event FITS files processed through NICERsoft

"""
from __future__ import division,print_function
from astropy.io import fits
import Lv0_dirs,Lv0_nicersoft_evt_filename

Lv0_dirs.global_par() #obtaining the global parameters

def open_fits(obsid,name_par_list):
    """
    Opening the FITS file for the cleaned event file

    obsid - Observation ID of the object of interest (10-digit str)
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range

    name_par_list should be [GTI_true,E_true,GTIno,segment_length,PI1,PI2]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(name_par_list) != list and type(name_par_list) != np.ndarray:
        raise TypeError("name_par_list should either be a list or an array!")
    if len(name_par_list) != 6:
        raise ValueError("There seems to be fewer or more values in the list/array than there should be! You should have [GTI_true, E_true, GTIno, segment length, PI1, PI2]")

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/' + Lv0_nicersoft_evt_filename.evt_filename(obsid,name_par_list)
    print(event)
    event = fits.open(event)
    #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return event

def get_eventcl(obsid,name_par_list,par_list):
    """
    Getting data from the FITS files, e.g., PI_FAST, TIME, PI, PI_RATIO, FLAGS

    obsid - Observation ID of the object of interest (10-digit str)
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)

    name_par_list should be [GTI_true,E_true,GTIno,segment_length,PI1,PI2]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(name_par_list) != list and type(name_par_list) != np.ndarray:
        raise TypeError("name_par_list should either be a list or an array!")
    if len(name_par_list) != 6:
        raise ValueError("There seems to be fewer or more values in the list/array than there should be! You should have [GTI_true, E_true, GTIno, segment length, PI1, PI2]")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = open_fits(obsid,name_par_list)
    data_dict = {}

    for i in range(len(par_list)):
        data_dict[par_list[i]] = event[1].data[par_list[i]]

    return data_dict

################################################################################
if __name__ == "__main__":
    datadict = get_eventcl('0034070101',[True,'',0,100,'',''],['PI','PI_FAST','TIME'])
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
