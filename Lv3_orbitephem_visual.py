#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 11:14am 2020

Given an orbital ephemeris (really just T0 and orbital period), show visually/graphically
where in the orbit an observation is (e.g., for NICER's NGC 300 ULX-1, whether it's in
eclipse of X-1 or otherwise)

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import Lv0_dirs,Lv2_dj_lsp,Lv2_swift_lc,Lv2_phase
import os
from scipy import stats
from scipy.optimize import curve_fit
from tqdm import tqdm
import subprocess
import mplcursors
import pathlib
from stingray.pulse.pulsar import pulse_phase,phase_exposure,fold_events

def folding(eventfile,Porb,nbins):
    """
    Folding the events by some orbital period
    """
    times = fits.open(eventfile)[1].data['TIME'] #getting array of times
    gtis_data = fits.open(eventfile)[2].data #getting GTIs
    T = sum([ gtis_data[i]['STOP']-gtis_data[i]['START'] for i in range(len(gtis_data)) ]) #exposure time

    gtis_conform = []
    for i in range(len(gtis_data)):
        gtis_conform.append([gtis_data[i][0],gtis_data[i][1]]) #conform to the input that Stingray uses

    phase_sr,prof_sr,err_sr = fold_events(times,1/Porb,gtis=np.array(gtis_conform),ref_time=times[0],nbin=nbins)
    phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times,1/Porb,gtis=np.array(gtis_conform),ref_time=times[0],expocorr=True,nbin=nbins)

    total_phase_sr = list(phase_sr) + list(phase_sr+1)
    total_prof_sr = list(prof_sr)*2
    total_err_sr = list(err_sr)*2

    total_phase_sr_expo = list(phase_sr_expo) + list(phase_sr_expo+1)
    total_prof_sr_expo = list(prof_sr_expo)*2
    total_err_sr_expo = list(err_sr_expo)*2

    plt.figure()
    plt.errorbar(x=total_phase_sr,y=total_prof_sr/T,yerr=total_err_sr/T,color='r',drawstyle='steps-mid')
    plt.errorbar(x=total_phase_sr_expo,y=total_prof_sr_expo/T,yerr=total_err_sr_expo/T,color='b',drawstyle='steps-mid')
    plt.legend(('Folded profile','Exposure-corrected'),loc='best',fontsize=12)
    plt.title(str(pathlib.Path(eventfile).name) +', exposure-corrected (using Stingray fold_events)',fontsize=12)
    plt.xlabel('Phase',fontsize=12)
    plt.ylabel('Counts/s',fontsize=12)

    return total_phase_sr_expo,total_prof_sr_expo/T,total_err_sr_expo/T

def ephemeris(x,y,yerr,eventfile,T0,Porb,nbins):
    """
    Plotting the light curve (usually), along with the ephemeris visually
    """

    ##### Figure 1 shows the folded profile ; using stingray.pulse.par's fold_events
    phase,prof,prof_err = folding(eventfile,Porb,nbins)

    ##### Figure 2 will superpose the light curve and the ephemeris
    plt.figure()
    plt.errorbar(x,y,yerr=yerr,color='r',fmt='x')
    plt.xlabel('Time (MJD)',fontsize=12)
    plt.ylim([np.min(y),1.1*np.max(y)])
    plt.ylabel('Rate (counts/s)',fontsize=12)

    intervals = np.arange(T0,x[-1]+Porb/86400,Porb/86400) #defining EACH orbital cycle, starting with T0
    for i in range(0,len(intervals)-1):
        subintervals = np.linspace(intervals[i],intervals[i+1],nbins+1)
        plt.axvline(x=intervals[i],color='k')
        plt.axvline(x=intervals[i+1],color='k')
        plt.errorbar(subintervals,prof[:nbins+1]*(0.5*np.max(y)/np.mean(prof)),yerr=prof_err[:nbins+1]*(0.5*np.max(y)/np.mean(prof)),color='b',drawstyle='steps-mid',alpha=0.5)

        for j in range(len(subintervals)):
            plt.axvline(x=subintervals[j],alpha=0.5,lw=0.5,color='k')

    plt.show()

if __name__ == "__main__":
    eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge_niceroverlap_all.evt'
    ##### Running ephemeris
    nbins = 20
    Porb = 1/8.4712e-6

    ##### Calling X-1 light curve data from Swift
    bary_outputfolder = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/lightcurve/'
    obsids = [str(i) for i in range(49834027,49834042)] + [str(i) for i in range(49834043,49834062)] + [str(i) for i in range(49834063,49834066)] + ['88810002'] + [str(i) for i in range(49834066,49834069)] + [str(i) for i in range(49834070,49834079)] + [str(i) for i in range(49834080,49834088)]
    corr_lc_files = [bary_outputfolder + 'sw000' + obsids[i] + '_corr.lc' for i in range(len(obsids))]
    corr_bg_files = [bary_outputfolder + 'sw000' + obsids[i] + '_bg_corr.lc' for i in range(len(obsids))]
    bg_scale_x1 = (30/120)**2
    rebinned_t, rebinned_rate, rebinned_err, rebinned_fracexp = Lv2_dj_lsp.rebin_lc(corr_lc_files,corr_bg_files,bg_scale_x1,100,0.0)

    tstart_49834027 = 546830295.758713
    mjd = fits.open(eventfile)[1].header['MJDREFI'] + fits.open(eventfile)[1].header['MJDREFF'] + (tstart_49834027+rebinned_t)/86400
    T0_MJD = fits.open(eventfile)[1].header['MJDREFI'] + fits.open(eventfile)[1].header['MJDREFF'] + (tstart_49834027)/86400

    #ephemeris(mjd,rebinned_rate,rebinned_err,eventfile,T0_MJD,Porb,nbins)

    ##### Calling NICER's ULX-1 observations
    s1,s2,a,b,c,d,inband,mjd = np.genfromtxt('/Volumes/Samsung_T5/n300_ulx_2020/n300_ulx.bgsub_cl50_g2020norm.fffphot',usecols=(3,4,5,6,7,8,9,11),unpack=True)
    ephemeris(mjd,inband,np.zeros(len(inband)),eventfile,T0_MJD,Porb,nbins)
