#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:56am 2019

Time to execute the programs!

"""
from __future__ import division, print_function
import numpy as np
import time

import Lv0_dirs
import Lv1_data_bin
import Lv2_lc,Lv2_ps,Lv2_color,Lv2_phase
import Lv3_E_boundary,Lv3_diagnostics

import matplotlib.pyplot as plt

Lv0_dirs.global_par() #obtaining the global parameters

### parameters used EVERYWHERE
basedir = Lv0_dirs.NICERSOFT_DATADIR
#eventfiles = [basedir + str(i) + '_pipe/ni' + str(i) + '_nicersoft_bary.evt' for i in ['0034070101']]
#eventfiles = [Lv0_dirs.NICER_DATADIR + 'rxj0209/rxj0209kgfilt_bary.evt']
eventfiles = ['/Volumes/Samsung_T5/NICER-data/1030180113/xti/event_cl/ni1030180113_0mpu7_ufa.evt']
par_list = ['PI','PI_FAST','TIME'] #parameter list from event_cl

tbin_size = 1 #how you want to bin the light curve data
Ebin_size = 0.05 #in keV
mode = 'show'
truncations = 'all' #'all', 't', 'E', or 'tE', depending on whether we want to look at entire time series (all), or truncation by time interval (t), or time truncation by energy range (E), or truncation by both (tE)

lc = False
ps = False
phase = False
color = False
###############################################################################

#### DEFINE DESIRED TIME INTERVALS AND ENERGY RANGES HERE FOR:
# Lv2_ps - partial_t, partial_E, partial_tE
# Lv2_phase - partial_t, partial_E, partial_tE
# Lv2_color - plotting_t

t1 = 39651
t2 = 39680
E1 = 0.3
E2 = 12
#0 1800 ; 5400-7360 ; 10960-12925 ; 16380-18500 ; 22050 - 22350 ; 28360 - 28600
###############################################################################

#for Lv2_ps
ps_type = 'both' # 'period' (for periodogram) or 'manual' (for FFT) or 'both'
oversampling = [False,5] # [False to NOT oversample, oversampling factor - 5 to oversample by factor of 5. (factor-1) sets of 0s are padded.]
xlims = [True,0,1000] # [False to NOT impose xlimit on plots; 2nd/3rd entries are the desired x-limits if needed.]
vlines = [True,581.127] # [False to NOT draw a vertical line on the plot; 2nd entry is the equation for the vertical line, e.g. x=2]

#for Lv2_phase
### For an unknown observation, one should run JUST Lv2_lc and Lv2_ps first to get
### the pulsation frequencies. Pulse profiles come LATER.
### If I have pulse_pars[1] and pulse_pars[2] != 0, then time binning DOES NOT MATTER, i.e., it'll be counts/s!
pulse_pars = [581.127,0,0]
shift = 0.4 # how much to shift the pulse by in the phase axis. It only affects how the pulse profile is 'displaced'.
no_phase_bins = 20 # number of phase bins desired

#for Lv2_color
E1_data = 0.3 #data is reliable between 0.3 and 12 keV
E2_data = 12 # in keV
cut_type = 'manual' # 'manual' cut for boundary energy, or 'median' - for half number of counts
bound = 2.7 # boundary energy for when cut_type = 'manual'!

for i in range(len(eventfiles)):
    E_bound = Lv3_E_boundary.E_bound(eventfiles[i],par_list,E1_data,E2_data,cut_type,bound) #use Lv3_E_boundary to get boundary energy
    ############################ FOR WHOLE OBSERVATION ############################
    if truncations == 'all':
        if lc == True:
            Lv2_lc.whole(eventfiles[i],par_list,tbin_size,mode) #light curve
            time.sleep(1)
        if ps == True:
            Lv2_ps.whole(eventfiles[i],par_list,tbin_size,mode,ps_type,oversampling,xlims,vlines) #power spectra
            time.sleep(1)
        if phase == True:
            Lv2_phase.whole(eventfiles[i],par_list,tbin_size,pulse_pars,shift,no_phase_bins,mode)
            time.sleep(1)
        if color == True:
            Lv2_color.plotting(eventfiles[i],par_list,E_bound,tbin_size,mode)

    ########################## FOR DESIRED TIME INTERVAL ##########################
    if truncations == 't':
        if lc == True:
            Lv2_lc.partial_t(eventfiles[i],par_list,tbin_size,t1,t2,mode) #light curve
            time.sleep(1)
        if ps == True:
            Lv2_ps.partial_t(eventfiles[i],par_list,tbin_size,t1,t2,mode,ps_type,oversampling,xlims,vlines) #power spectra
            time.sleep(1)
        if phase == True:
            Lv2_phase.partial_t(eventfiles[i],par_list,tbin_size,pulse_pars,shift,no_phase_bins,t1,t2,mode)
            time.sleep(1)
        if color == True:
            Lv2_color.plotting_t(eventfiles[i],par_list,E_bound,tbin_size,t1,t2,mode)

    ########################### FOR DESIRED ENERGY RANGE ##########################
    # won't anticipate that this will be used much?
    if truncations == 'E':
        if lc == True:
            Lv2_lc.partial_E(eventfiles[i],par_list,tbin_size,Ebin_size,E1,E2,mode)
            time.sleep(1)
        if ps == True:
            Lv2_ps.partial_E(eventfiles[i],par_list,tbin_size,Ebin_size,E1,E2,mode,ps_type,oversampling,xlims,vlines)
            time.sleep(1)
        if phase == True:
            Lv2_phase.partial_E(eventfiles[i],par_list,tbin_size,Ebin_size,pulse_pars,shift,no_phase_bins,E1,E2,mode)

    ################# FOR DESIRED TIME INTERVAL AND ENERGY RANGE #################
    if truncations == 'tE':
        if lc == True:
            Lv2_lc.partial_tE(eventfiles[i],par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode)
            time.sleep(1)
        if ps == True:
            Lv2_ps.partial_tE(eventfiles[i],par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode,ps_type,oversampling,xlims,vlines)
            time.sleep(1)
        if phase == True:
            Lv2_phase.partial_tE(eventfiles[i],par_list,tbin_size,Ebin_size,pulse_pars,shift,no_phase_bins,t1,t2,E1,E2,mode)
            time.sleep(1)
        if color == True:
            Lv2_color.plotting_t(eventfiles[i],par_list,E_bound,tbin_size,t1,t2,mode)

################################################################################
################################# DIAGNOSTICS ##################################
################################################################################

att_var = ['TIME','MODE','SUBMODE_AZ','SUBMODE_EL']
mkf_var = ['TIME','ELV', 'BR_EARTH', 'SUNSHINE', 'FOV_FLAG', 'SUN_ANGLE',
           'MOON_ANGLE', 'ANG_DIST', 'SAA', 'SAA_TIME', 'COR_ASCA', 'COR_SAX',
           'MCILWAIN_L', 'NICER_SAA', 'TOT_ALL_COUNT', 'TOT_UNDER_COUNT',
           'TOT_OVER_COUNT', 'TOT_XRAY_COUNT']
hk_var = ['TIME','MPU_D_TEMP','MPU_A_TEMP','MPU_PWRBRDG_TEMP']
par_list = ['TIME','PI','PI_FAST','PI_RATIO']

diag_vars = {}
diag_vars['att'] = att_var
diag_vars['mkf'] = mkf_var
diag_vars['hk'] = hk_var

# [Best to do mode='save' with the many plots]
# Getting diagnostic plots over the entire observation; looks at how variables
# like TOT_OVER_COUNT changes over time!
mode = 'save'
for i in range(len(eventfiles)):
    Lv3_diagnostics.diag_all(eventfiles[i],par_list,tbin_size,mode,diag_vars)

# [Best to do mode='save' with the many plots]
# Getting diagnostic plots over the desired time interval; looks at how variables
# like TOT_OVER_COUNT changes over time for a desired time interval!
#Lv3_diagnostics.diag_t(obsid,bary,par_list,tbin_size,t1,t2,mode,diag_vars)

"""
###############################################################################
###############################################################################
###############################################################################
#### do for March 25 2019 for now - see broad overview of pulse shapes
#oh bless, 'overlap' works exactly as I want it to!
f_pulses = [0.101464,0.10151,0.101738,0.101993,0.102029,0.102029,0.102082,0.102163,0.102103,0.102109,0.102158,0.101977]
subplot_Es = [(0.2,1),(1,2),(2,3),(3,5),(5,8),(8,12)]
if type(obsids) == list or type(obsids) == np.array:
    for i in range(len(obsids)):
        Lv2_phase.partial_subplots_E(eventfiles[i],par_list,tbin_size,Ebin_size,f_pulses[i],shift,no_phase_bins,subplot_Es,E1,E2,mode)
#        Lv2_phase.partial_E(eventfiles[i],par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,E1,E2,mode)
    plt.show()
"""
### ONE SHOULD NOT NEED TO EDIT THE SCRIPTS BEYOND THIS POINT, UNLESS YOU KNOW
### WHAT YOU'RE DOING. ARGUABLY, I DON'T EITHER, BUT I'LL TRY MY BEST.
