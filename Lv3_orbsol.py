#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 8 2:20pm 2019

Obtaining the orbital solution for some given ObsID.

"""
from __future__ import division, print_function
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os
from astropy.time import Time

### use the info from Bildsten 1997!

#orbital epoch
E_0 = 48561.65670200 # T_pi/2 , epoch of 90 deg. mean orbital longitude
E_0_error = 0.00000071

P_orb = 2.0870653300
P_orb_error = 0.0000000049

P_orb_dot = -9.93e-9 #days/day ; from Nagase et al. 1992
P_orb_dot_error = 0.02e-9 #days/day

def get_phase(obs_epoch):
    obs_epoch_obj = [obs_epoch]
    t = Time(obs_epoch_obj,format='isot',scale='utc')
    obs_epoch_mjd = t.mjd[0]

    obs_phase_day = (obs_epoch_mjd-E_0)%P_orb #modulo!
    obs_phase = obs_phase_day/P_orb

    return obs_phase

starts = ['2017-06-20T19:06:45', '2017-06-21T02:43:09', '2017-06-22T14:28:56', '2017-06-23T07:34:02', '2017-07-24T23:13:18', '2018-11-07T15:23:04']
stops = ['2017-06-20T19:23:53','2017-06-21T18:33:37', '2017-06-22T18:59:23', '2017-06-23T18:08:23', '2017-07-24T23:21:18', '2018-11-07T15:40:03']

for i in range(len(starts)):
    print(get_phase(starts[i])*360,get_phase(stops[i])*360)
