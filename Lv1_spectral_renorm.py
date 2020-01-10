#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues Feb 5 1.19pm 2019

Renormalizing the NICER spectral data by the NICER observations of the Crab Nebula.
This is to mitigate residuals.

"""
from __future__ import division, print_function
import numpy as np
import time

from astropy.io import fits
import Lv0_dirs

Lv0_dirs.global_par()

def renorm_spectra(obsid,bary):
    """
    Renormalizing the NICER spectral data with the 'ratios of residuals' sent by
    Jack. See email from 4 Feb 2019.

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")

    renorm_file = Lv0_dirs.NICER_DATADIR + 'crabcor.dat' #file name for the renormalization ratios
    chan,ratio = np.loadtxt(renorm_file,usecols=(0,1),unpack=True) #loading the channel and ratio values

    orig_spectral = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_cl/' + obsid+'_spec.pi' #spectral file you want to renormalize
    event = fits.open(orig_spectral)

    orig_chan = event[1].data['CHANNEL'] #extracting the channel values from the spectral file
    orig_counts = event[1].data['COUNTS'] #extracting the count rate from the spectral file

    if len(orig_chan) != len(chan):
        return ValueError(obsid+ ": The number of channels from the renormalization ratio file is not equal to that of the number of channels in the spectral file. Make the necessary alterations!")

    renormalized = orig_counts/ratio
    event[1].data['COUNTS'] = renormalized

    event.writeto(Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_cl/' + obsid + '_renormspec.pi',overwrite=True)

if __name__ == "__main__":
    obsids = ['0034070101']#,'0034070102','0034070103','0034070104','1034070101','1034070102','1034070103','1034070104','1034070105','1034070106']
    for i in range(len(obsids)):
        renorm_spectra(obsids[i],True)
