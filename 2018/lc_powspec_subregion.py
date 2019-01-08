#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 17:21:04 2018

@author: masonng

To produce:
    
1st page – Light curve for WHOLE subregion
2nd page – Power spectrum of WHOLE subregion
3rd page – Light curve for SPECIFIC SECTION OF SUBREGION – sub-subregion
4th page – Light curve for SPECIFIC SECTION OF SUBREGION – sub-subregion
"""
from __future__ import division
import datetime
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy import signal

timestart = time.time()

print(datetime.datetime.now())

work_dir = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/'
obsid = '1034070104' 
obsname = 'Cen_X-3'
t_interval = '0p01s'
place = 100 #for t_interval = 0p1s

def open_fits(work_dir, obsid):
    ## use barycentered data
    
    event = work_dir + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl_bary.evt'
    event = fits.open(event) #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return event

def get_data(work_dir,obsid):
    event = open_fits(work_dir,obsid)
    
    #recall: energy stuff is given in PI, so don't do PI_FAST truncation
    pi_fast = event[1].data['PI_FAST']
#    print 'Done with PI_FAST'
    times = event[1].data['TIME']
#    print 'Done with TIME'
    pi_data = event[1].data['PI']
#    print 'Done with PI'
    pi_ratio = event[1].data['PI_RATIO']
#    print 'Done with PI_RATIO'
    flags = event[1].data['EVENT_FLAGS']
#    print 'Done with FLAGS'
    
    return pi_fast, times, pi_data, pi_ratio, flags

def lightcurve(work_dir,obsid,doplot,impose_xlim,xlim1,xlim2):
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    
    shifted_t = times-times[0]
    counts_data = np.ones(len(pi_data))
        
    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts_data,statistic='sum',bins=t_bins)
    
    if doplot==True:
        plt.figure(figsize=(10,8))
        plt.plot(t_bins[:-1],summed_data,'r-')
        if impose_xlim==True:
            plt.xlim([xlim1,xlim2])
            plt.ylim([min(summed_data[xlim1:xlim2]),max(summed_data[xlim1:xlim2])])
            plt.title('Light curve in 1s bins for t='+str(xlim1)+'-'+str(xlim2)+'\n ObsID : ' + str(obsid),fontsize=15)
        else:
            plt.title('Light curve in 1s bins for whole time series, t=0-'+str(t_bins[-1])+'\n ObsID: ' + str(obsid), fontsize=15)
        plt.xlabel('Time (s)',fontsize=15)
        plt.ylabel('Counts/s',fontsize=15)

    return t_bins[:-1],summed_data

def get_txt_data(work_dir,obsid,t_interval):
    """
    Returns the barycentered data for the times and corresponding counts, where
    the time interval is t_interval.
    """
    
    times, counts = np.loadtxt(work_dir+'sumcounts_'+obsid+'_'+t_interval+'_bary.txt', usecols=(0,1), unpack=True) 
    return times, counts

def indiv_pow_spec(work_dir,obsid,t_interval,cut1,cut2):
    #takes as input the light curve (times,counts), as well as the indices
    #(cut1,cut2) for the truncated light curves to obtain the power spectra
    #individually
    times,counts = get_txt_data(work_dir,obsid,t_interval)
    
    times_trunc = times[cut1*100:cut2*100]
    counts_trunc = counts[cut1*100:cut2*100]
    
    T = times_trunc[-1]-times_trunc[0]
    dt = T/len(times_trunc)
    
    mean_corrected = counts_trunc - np.mean(counts_trunc)
    freqs = np.fft.fftfreq(mean_corrected.size, dt)
    N = len(freqs)

    power_spec = 2.0/sum(counts_trunc)**2*np.abs(np.fft.fft(mean_corrected))**2
    
    return freqs, power_spec, N

def plot_powspec(freqs,power_spec,N,indicator,xlim1,xlim2):
    
    plt.semilogy(freqs[1:int(N/2)],power_spec[1:int(N/2)]/np.mean(power_spec[1:int(N/2)]),'r-')
    plt.xlabel('Hz')
    plt.ylabel('Some normalized power spectrum')
    plt.xlim([xlim1,xlim2])
    plt.axvline(x=indicator,color='k',alpha=0.5,lw=0.5)

def obtain_pdf(obsid,obsname,t_interval,subr_cut1,subr_cut2,cut1,cut2,indicator,ps_xlim1,ps_xlim2):
    
    filename = work_dir + obsid + '/' + obsid + '_' + obsname + '_' + str(subr_cut1) + '-' +str(subr_cut2) + 's.pdf'
    with PdfPages(filename) as pdf:
    #### Truncated time series
        lightcurve_t, lightcurve_data = lightcurve(work_dir,obsid,True,True,subr_cut1,subr_cut2) #but it's ok as impose_xlim is False
        pdf.savefig()
        plt.close()
        print 'Whole light curve for t=' + str(subr_cut1) + '-' + str(subr_cut2) + 's'
        
        freqs_all, powspec_all, N_all = indiv_pow_spec(work_dir,obsid,t_interval,subr_cut1,subr_cut2)
        plot_powspec(freqs_all,powspec_all,N_all,indicator,ps_xlim1,ps_xlim2)
        pdf.savefig()
        plt.close()
        print 'Whole power spectrum done for t=' + str(subr_cut1) + '-' + str(subr_cut2) + 's'

        lightcurve_t, lightcurve_data = lightcurve(work_dir,obsid,True,True,cut1,cut2)
        pdf.savefig()
        plt.close()
        print 'Subregion light curve for t=' + str(cut1) + '-' + str(cut2) + 's'
        
        freqs_sub, powspec_sub, N_sub = indiv_pow_spec(work_dir,obsid,t_interval,cut1,cut2)
        plot_powspec(freqs_sub,powspec_sub,N_sub,indicator,ps_xlim1,ps_xlim2)
        pdf.savefig()
        plt.close()
        print 'Power spectrum for t=' + str(cut1) + '-' + str(cut2) + 's'
        

## for 1034070104 - 11113, 11945 ; 11400, 11500,
## for 1034070102 - 0,620, 380, 620
obtain_pdf(obsid,obsname,t_interval,11113,11945,11113,11700,0.20865682023238702,0,1)
        
timeend = time.time()

print (timeend-timestart)