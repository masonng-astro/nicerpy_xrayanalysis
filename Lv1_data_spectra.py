#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 9 10:23am 2019

Obtain the spectra (photons/s) from the counts/s. This would involve folding
in the response matrix. Recall this is found in ciri/mason - placed by Jack.
"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv0_call_eventcl

Lv0_dirs.global_par() #obtaining the global parameters

def read_redist():
    """
    Opening the FITS file corresponding to the redistribution matrix

    Describes the detector response

    Has CHANNEL, E_MIN, E_MAX for first card
    Has ENERG_LO, ENERG_HI, N_GRP, F_CHAN, N_CHAN, MATRIX
    """
    redist = Lv0_dirs.BASE_DIR + 'nicer-data/nicer.rmf' #reading in redistribution matrix file (RMF)
    redist = fits.open(redist)
    #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return redist

def read_anc():
    """
    Opening the FITS file corresponding to the ancillary response

    Describes the efficiency vs energy; gives telescope area x filter efficiency
    x detector quantum efficiency vs energy

    Has ENERG_LO, ENERG_HI, SPECRESP, ENERGY, XRCAREA, QE, WINDOW, THERMALSD
    """
    anc = Lv0_dirs.BASE_DIR + 'nicer-data/nicer.arf' #reading in ancillary response file (ARF)
    anc = fits.open(anc)
    #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return anc
