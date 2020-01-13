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

def fits2dict(fits_file,ext,par_list):
    """
    'Converts' a FITS file to a Python dictionary, with a list of the original
    FITS table's columns, of the user's choosing. Can use this for ANY FITS file,
    but is meant for mkf/orb files in $OBSID_pipe folders from NICER, or event files
    (be it from NICER-data or NICERsoft_outputs or any other mission which has this format)

    fits_file - path to the FITS file
    ext - which extension number; 1 for EVENTS, 2 for GTI, 3 for PPS_TREND
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI)
    """
    if type(fits_file) != str:
        raise TypeError("fits_file should be a string!")
    if type(ext) != int:
        raise TypeError("ext should be an integer!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")

    event = fits.open(fits_file)

    fits_dict = {par_list[i]:event[ext].data[par_list[i]] for i in range(len(par_list))}
    #returns a dictionary in the form {'TIME':[array of TIME values],'VARIABLE':[array of VARIABLE values]}

    return fits_dict

if __name__ == "__main__":
    eventfile = "/Volumes/Samsung_T5/NICERsoft_outputs/1034070101_pipe/cleanfilt.evt"
    #print(fits2dict(eventfile,1,['TIME','RAWX']))
    print(fits2dict(eventfile,2,['START','STOP']))
