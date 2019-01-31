#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 7:13pm 2019

Used to find the peak frequency

"""
from __future__ import division, print_function
import numpy as np
import time
import matplotlib.pyplot as plt
import Lv0_call_eventcl,Lv2_ps,Lv2_ps_method,Lv2_phase
from scipy import stats,signal

### step 1: plot power spectrum from periodogram and manual method; just to visually check them.
### step 2: use manual method to get the estimated peak frequency oversampled by a factor of 5!
### step 3: zoom in to the region hosting the peak frequency
### step 4: find the peak frequency and note the bin number
### step 5: normalize the power P by the local noise power... where \bar{P}_k = P_k/mean(P)
### step 6: get the 'best estimate of the signal frequency' and the uncertainty from Eq 3.11 and 3.12 in Deepto's PhD thesis

###############################################################################

### parameters used EVERYWHERE
obsid = '1034070106' #observation ID.
#obsids = ['0034070101','0034070102','0034070103','0034070104','1034070101','1034070102','1034070103','1034070104','1034070105','1034070106']
bary = True #whether the data you want is barycenter-corrected or not
par_list = ['PI','PI_FAST','TIME'] #parameter list from event_cl
tbin_size = 1 #how you want to bin the light curve data
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

#for Lv2_ps
ps_type = 'both' # 'period' (for periodogram) or 'manual' (for FFT) or 'both'
oversampling = [True,5] # [False to NOT oversample, oversampling factor - 5 to oversample by factor of 5. (factor-1) sets of 0s are padded.]
xlims = [True,0,5] # [False to NOT impose xlimit on plots; 2nd/3rd entries are the desired x-limits if needed.]
vlines = [True,0.2084] # [False to NOT draw a vertical line on the plot; 2nd entry is the equation for the vertical line, e.g. x=2]
"""
############################ FOR WHOLE OBSERVATION ############################
if truncations == 'all':
    if ps_type == 'both':
        pdgm_f,pdgm_ps,manual_f,manual_ps = Lv2_ps.whole(obsid,bary,par_list,tbin_size,mode,ps_type,oversampling,xlims,vlines) #power spectra
    else:
        f,ps = Lv2_ps.whole(obsid,bary,par_list,tbin_size,mode,ps_type,oversampling,xlims,vlines) #power spectra

########################## FOR DESIRED TIME INTERVAL ##########################
if truncations == 't':
    if ps_type == 'both':
        pdgm_f,pdgm_ps,manual_f,manual_ps = Lv2_ps.partial_t(obsid,bary,par_list,tbin_size,t1,t2,mode,ps_type,oversampling,xlims,vlines) #power spectra
    else:
        f,ps = Lv2_ps.partial_t(obsid,bary,par_list,tbin_size,t1,t2,mode,ps_type,oversampling,xlims,vlines) #power spectra

########################### FOR DESIRED ENERGY RANGE ##########################
# won't anticipate that this will be used much?
if truncations == 'E':
    if ps_type == 'both':
        pdgm_f,pdgm_ps,manual_f,manual_ps = Lv2_ps.partial_E(obsid,bary,par_list,tbin_size,Ebin_size,E1,E2,mode,ps_type,oversampling,xlims,vlines)
    else:
        f,ps = Lv2_ps.partial_E(obsid,bary,par_list,tbin_size,Ebin_size,E1,E2,mode,ps_type,oversampling,xlims,vlines)

################# FOR DESIRED TIME INTERVAL AND ENERGY RANGE #################
if truncations == 'tE':
    if ps_type == 'both':
        pdgm_f,pdgm_ps,manual_f,manual_ps = Lv2_ps.partial_tE(obsid,bary,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode,ps_type,oversampling,xlims,vlines)
    else:
        f,ps = Lv2_ps.partial_tE(obsid,bary,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode,ps_type,oversampling,xlims,vlines)
"""
################################### STEP 2+ ###################################
"""
harmonic = 2 #harmonic number used!
data_dict = Lv0_call_eventcl.get_eventcl(obsid,bary,par_list)
times = data_dict['TIME'] #extract the timestamps
shifted_t = times-times[0] #so that we start at t = 0
T = shifted_t[-1]

max_ps = max(manual_ps[(manual_f>=0.40)&(manual_f<=0.45)])
print(np.mean(manual_ps[int(0.9*len(manual_f)):]))
index = np.where(manual_ps==max_ps)[0][0]
print(manual_f[index])
print(max_ps)
print(manual_ps[index+oversampling[1]])
best_f = manual_f[index] + 3/(4*np.pi**2*5*T) * (manual_ps[index+oversampling[1]]-manual_ps[index-oversampling[1]])/manual_ps[index]
print('The best estimate for the signal frequency is: ' + str(best_f/harmonic))
error = 1/(2*np.pi*T) * np.sqrt(6/(manual_ps[index]/np.mean(manual_ps[(manual_f>=0.40)&(manual_f<=0.45)])))
print('The error is: ' + str(error/harmonic))
"""
################################################################################
### Fold the time series by the average period and manually

f_ave = 0.20826334005309413
T_ave = 1/f_ave

obsids = ['0034070101','0034070102','0034070103','0034070104','1034070101','1034070106']
f = np.array([0.20846118761251825,0.2080718358508059,0.20854506857953284,0.20796037636355533,0.20799884640430513,0.20854272550784736])
for i in range(len(obsids)):
    data_dict = Lv0_call_eventcl.get_eventcl(obsids[i],bary,par_list)
    times = data_dict['TIME']
    counts = np.ones(len(times))
    print(len(times),len(counts))
    shifted_t = times-times[0]
    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])*1/tbin_size+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

    phase_bins, pulse_profile = Lv2_phase.pulse_profile(f_ave,times,counts,0.5,40)
    print(sum(pulse_profile))
    phase_bins_f, pulse_profile_f = Lv2_phase.pulse_profile(f[i],times,counts,0.5,40)

    plt.figure(i)
    plt.title(obsids[i],fontsize=12)
    plt.xlabel('Phase',fontsize=12)
    plt.ylabel('Count/'+str(tbin_size)+'s')
    plt.plot(phase_bins[:-1],pulse_profile,'rx-',alpha=0.5)
    plt.plot(phase_bins_f[:-1],pulse_profile_f,'bx-',alpha=0.5)
    plt.legend(('Average f','Pulse f',),loc='best')

plt.show()
