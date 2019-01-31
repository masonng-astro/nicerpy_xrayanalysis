#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:56am 2019

Time to execute the programs!

"""
from __future__ import division, print_function
import numpy as np
import time

import Lv0_dirs
import Lv0_call_eventcl,Lv0_call_att,Lv0_call_hk,Lv0_call_mkf,Lv0_call_orb,Lv0_call_uf,Lv0_call_ufa
import Lv1_data_bin,Lv1_data_gtis,Lv1_data_spectra
import Lv2_lc,Lv2_ps,Lv2_color,Lv2_phase
import Lv3_E_boundary,Lv3_diagnostics

import matplotlib.pyplot as plt

### parameters used EVERYWHERE
obsid = '0034070102' #observation ID.
#obsids = ['0034070101','0034070102','0034070103','0034070104','1034070101','1034070102','1034070103','1034070104','1034070105','1034070106']
bary = True #whether the data you want is barycenter-corrected or not
par_list = ['PI','PI_FAST','TIME'] #parameter list from event_cl
tbin_size = 0.1 #how you want to bin the light curve data
Ebin_size = 0.05 #in keV
mode = 'show' # 'show' the plots or 'save' the plots

truncations = 'all' #'all', 't', 'E', or 'tE', depending on whether we want to look at entire time series (all), or truncation by time interval (t), or time truncation by energy range (E), or truncation by both (tE)

###############################################################################

#### DEFINE DESIRED TIME INTERVALS AND ENERGY RANGES HERE FOR:
# Lv2_ps - partial_t, partial_E, partial_tE
# Lv2_phase - partial_t, partial_E, partial_tE
# Lv2_color - plotting_t

t1 = 0
t2 = 810
E1 = 0.3
E2 = 2.7

###############################################################################

### more 'obscure' parameters
#for Lv1_data_gtis
gap = 50

#for Lv2_ps
ps_type = 'both' # 'period' (for periodogram) or 'manual' (for FFT) or 'both'
oversampling = [False,5] # [False to NOT oversample, oversampling factor - 5 to oversample by factor of 5. (factor-1) sets of 0s are padded.]
xlims = [True,0,1] # [False to NOT impose xlimit on plots; 2nd/3rd entries are the desired x-limits if needed.]
vlines = [True,0.2084] # [False to NOT draw a vertical line on the plot; 2nd entry is the equation for the vertical line, e.g. x=2]

#for Lv2_phase
### For an unknown observation, one should run JUST Lv2_lc and Lv2_ps first to get
### the pulsation frequencies. Pulse profiles come LATER.
f_pulse = 0.2080718358508059 #frequency of the pulse
shift = 0.4 # how much to shift the pulse by in the phase axis. It only affects how the pulse profile is 'displaced'.
no_phase_bins = 51 # number of phase bins desired

#for Lv2_color
E1_data = 0.3 #data is reliable between 0.3 and 12 keV
E2_data = 12 # in keV
cut_type = 'manual' # 'manual' cut for boundary energy, or 'median' - for half number of counts
bound = 2.7 # boundary energy for when cut_type = 'manual'!
E_bound = Lv3_E_boundary.E_bound(obsid,bary,par_list,E1_data,E2_data,cut_type,bound) #use Lv3_E_boundary to get boundary energy

### first get GTIs for the observation
gti_array = Lv1_data_gtis.get_gtis(obsid,bary,gap)
#print(gti_array)
"""
for i in range(len(obsids)):
    print(obsids[i])
    print(Lv1_data_gtis.get_gtis(obsids[i],bary,gap))
"""
# is in the form: [gti_1_start,gti_1_stop,gti_2_start,gti_2_stop,...]

############################ FOR WHOLE OBSERVATION ############################
if truncations == 'all':
    Lv2_lc.whole(obsid,bary,par_list,tbin_size,mode) #light curve
    time.sleep(1)
    Lv2_ps.whole(obsid,bary,par_list,tbin_size,mode,ps_type,oversampling,xlims,vlines) #power spectra
    time.sleep(1)
    Lv2_phase.whole(obsid,bary,par_list,tbin_size,f_pulse,shift,no_phase_bins,mode)
    time.sleep(1)
    Lv2_color.plotting(obsid,bary,par_list,E_bound,tbin_size,mode)

########################## FOR DESIRED TIME INTERVAL ##########################
if truncations == 't':
    Lv2_lc.partial_t(obsid,bary,par_list,tbin_size,t1,t2,mode) #light curve
    time.sleep(1)
    Lv2_ps.partial_t(obsid,bary,par_list,tbin_size,t1,t2,mode,ps_type,oversampling,xlims,vlines) #power spectra
    time.sleep(1)
    Lv2_phase.partial_t(obsid,bary,par_list,tbin_size,f_pulse,shift,no_phase_bins,t1,t2,mode)
    time.sleep(1)
    Lv2_color.plotting_t(obsid,bary,par_list,E_bound,tbin_size,t1,t2,mode)

########################### FOR DESIRED ENERGY RANGE ##########################
# won't anticipate that this will be used much?
if truncations == 'E':
    Lv2_lc.partial_E(obsid,bary,par_list,tbin_size,Ebin_size,E1,E2,mode)
    time.sleep(1)
    Lv2_ps.partial_E(obsid,bary,par_list,tbin_size,Ebin_size,E1,E2,mode,ps_type,oversampling,xlims,vlines)
    time.sleep(1)
    Lv2_phase.partial_E(obsid,bary,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,E1,E2,mode)

################# FOR DESIRED TIME INTERVAL AND ENERGY RANGE #################
if truncations == 'tE':
    Lv2_lc.partial_tE(obsid,bary,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode)
    time.sleep(1)
    Lv2_ps.partial_tE(obsid,bary,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode,ps_type,oversampling,xlims,vlines)
    time.sleep(1)
    Lv2_phase.partial_tE(obsid,bary,par_list,tbin_size,Ebin_size,f_pulse,shift,no_phase_bins,t1,t2,E1,E2,mode)
    time.sleep(1)
    Lv2_color.plotting_t(obsid,bary,par_list,E_bound,tbin_size,t1,t2,mode)

###############################################################################
################################# DIAGNOSTICS #################################
###############################################################################

att_var = ['TIME','MODE','SUBMODE_AZ','SUBMODE_EL']
mkf_var = ['TIME','ELV', 'BR_EARTH', 'SUNSHINE', 'FOV_FLAG', 'SUN_ANGLE',
           'MOON_ANGLE', 'ANG_DIST', 'SAA', 'SAA_TIME', 'COR_ASCA', 'COR_SAX',
           'MCILWAIN_L', 'NICER_SAA', 'TOT_ALL_COUNT', 'TOT_UNDER_COUNT',
           'TOT_OVER_COUNT', 'TOT_XRAY_COUNT']
hk_var = ['TIME','MPU_D_TEMP','MPU_A_TEMP','MPU_PWRBRDG_TEMP']
eventcl_var = ['TIME','DEADTIME','MPU_A_TEMP','MPU_UNDER_COUNT','PI_RATIO']

diag_vars = {}
diag_vars['att'] = att_var
diag_vars['mkf'] = mkf_var
diag_vars['hk'] = hk_var
diag_vars['cl'] = eventcl_var

# [Best to do mode='save' with the many plots]
# Getting diagnostic plots over the entire observation; looks at how variables
# like TOT_OVER_COUNT changes over time!
#Lv3_diagnostics.diag_all(obsid,bary,par_list,tbin_size,mode,diag_vars)

# [Best to do mode='save' with the many plots]
# Getting diagnostic plots over the desired time interval; looks at how variables
# like TOT_OVER_COUNT changes over time for a desired time interval!
#Lv3_diagnostics.diag_t(obsid,bary,par_list,tbin_size,t1,t2,mode,diag_vars)


###############################################################################
###############################################################################
###############################################################################

### ONE SHOULD NOT NEED TO EDIT THE SCRIPTS BEYOND THIS POINT, UNLESS YOU KNOW
### WHAT YOU'RE DOING. ARGUABLY, I DON'T EITHER, BUT I'LL TRY MY BEST.
