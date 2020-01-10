#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:18am 2019

Finding the boundary energies

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import Lv0_dirs,Lv1_data_filter
from scipy import stats
from PyAstronomy.pyasl import foldAt
import matplotlib.pyplot as plt
import os

def E_bound(eventfile,par_list,E1,E2,cut_type,bound):
    """
    Gives the energy bound corresponding to either a custom cut or a median cut.
    Could add more cuts in the future!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    par_list - A list of parameters we'd like to extract from the FITS file
    (e.g., from eventcl, PI_FAST, TIME, PI,)
    >> e.g., tbin_size = 2 means bin by 2s
    >> e.g., tbin_size = 0.05 means bin by 0.05s!
    E1 - energy value for the lower boundary (in keV)
    E2 - energy value for the upper boundary (in keV)
    cut_type - 'manual' or 'median'
    bound - boundary energy for when cut_type = 'manual'
    """
    if 'PI' and 'TIME' not in par_list:
        raise ValueError("You should have BOTH 'PI' and 'TIME' in the parameter list!")
    if type(par_list) != list and type(par_list) != np.ndarray:
        raise TypeError("par_list should either be a list or an array!")
    if (E1<0) or (E2>20):
        raise ValueError("You're asking for boundaries <0 keV or > 20 keV. Do check!")
    if E2<E1:
        raise ValueError("E2 should be greater than E1!")
    if cut_type != 'manual' and cut_type != 'median':
        raise ValueError("Should be 'manual' or 'median', or time to add a new type of cut into Lv3_E_boundary!")

    if cut_type == 'manual':
        return bound

    if cut_type == 'median':
        t_cut,E_cut = Lv1_data_filter.filter_energy(eventfile,par_list,E1,E2)
        boundary_index = int(len(E_cut)/2)
        boundary_E = E_cut[boundary_index]

        return boundary_E

if __name__ == "__main__":
    obsid = '1034070101'
    eventfile = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl_bary.evt'
    print(E_bound(eventfile,['PI','TIME','PI_FAST'],0.0,20,'median',2.7))
