#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 2 5:55pm 2019

Obtaining the file name for event files in NICERsoft directories...

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
import glob

Lv0_dirs.global_par()

def evt_filename(obsid,name_par_list):
    """
    Getting the desired file name for the event file.

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

    prefix = 'ni' + obsid + '_nicersoft_bary'
    suffix = '.evt'
    if name_par_list[0] == True: #i.e., if we have time segments
        add1 = '_GTI' + str(name_par_list[2]) + '_' + str(name_par_list[3]) + 's'
    else:
        add1 = ''

    if name_par_list[1] == True: #i.e., if we have energy segments
        add2 = '_E' + str(name_par_list[4]) + '-' + str(name_par_list[5])
    else:
        add2 = ''

    return prefix + add1 + add2 + suffix

#print(evt_filename('0034070101',[True,True,1,500,200,800]))
