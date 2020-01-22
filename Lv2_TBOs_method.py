#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 14 1:36pm 2020

Script that operates on data so that it's ready for burst oscillation searches

"""
from __future__ import division, print_function
import numpy as np
import subprocess
from astropy.io import fits
from scipy import stats
from tqdm import tqdm
import pathlib

import Lv0_dirs,Lv0_fits2dict,Lv2_ps_method,Lv3_detection_level
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.tri as tri
from matplotlib import cm
import mplcursors

Lv0_dirs.global_par()

def burst_cands(eventfile):
    """
    Identifies 1s bins which are burst candidates by adopting the method from Galloway+ 08;
    that is, looking for bins where the number of counts exceeds 4 sigma of the overall mean

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    """
    parent_folder = str(pathlib.Path(eventfile).parent)

    event_header = fits.open(eventfile)[1].header
    obj_name = event_header['OBJECT']
    obsid = event_header['OBS_ID']

    data_dict = Lv0_fits2dict.fits2dict(eventfile,1,['TIME'])
    times = data_dict['TIME']
    counts = np.ones(len(times))

    gtis = Lv0_fits2dict.fits2dict(eventfile,2,['START','STOP'])
    starts = gtis['START']
    stops = gtis['STOP']
    T_obs = sum(np.array([ (gtis['STOP'][i]-gtis['START'][i]) for i in range(len(gtis)) ]))

    shifted_t = times-times[0]
    t_bins = np.linspace(0,np.ceil(shifted_t[-1]),np.ceil(shifted_t[-1])*1/1+1)
    binned_counts, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

    mean_counts = np.mean(binned_counts)
    std_counts = np.std(binned_counts)

    burst_cand_text = open(parent_folder+'/' + obsid + '_burst_cands.txt','w')
    burst_cand_times = []
    for i in range(len(binned_counts)):
        if binned_counts[i] >= mean_counts + 3*std_counts: #e.g., mean + 4 sigma
            burst_cand_times.append(t_bins[i])
            burst_cand_text.write(str(t_bins[i]) + ' ')
            if i%10 == 0:
                burst_cand_text.write('\n')
    burst_cand_text.close()

    if len(burst_cand_times) != 0:
        with PdfPages(parent_folder+'/'+obsid+'_burst_cands.pdf') as pdf:
            plt.figure()
            plt.plot(burst_cand_times,np.ones(len(burst_cand_times))*np.max(binned_counts),'rx')
            plt.plot(t_bins[:-1],binned_counts,'b-')
            plt.legend(('Burst candidate bins','Data'),loc='best',fontsize=12)
            plt.title('Burst candidates for ' + obj_name + ' with ObsID ' + str(obsid) + '\n Mean counts: ' + str(round(mean_counts,2)) + ' ; std: ' + str(round(std_counts,2)),fontsize=12)
            plt.xlabel('Time (s)',fontsize=12)
            plt.xlim([burst_cand_times[0]-10,burst_cand_times[-1]+10])
            plt.ylabel('Count rate in counts/s',fontsize=12)
            pdf.savefig()
            plt.close()

    return

def dynamic_ps(eventfile,search_window,T,dt,tbin_size,df,f_central,mode):
    """
    Plotting the dynamic power spectrum with both a colormap and a contour map.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    search_window - array of two values: [start time for burst searches, end time for burst searches]
    T - array of window sizes (not time interval)
    dt - array of time steps between consecutive time windows
    tbin_size - size of time bins
    df - frequency window width for the search
    f_central - central frequency of the search
    mode - "save" or "show"
    """
    if mode != 'show' and mode != 'save':
        raise ValueError("Mode should either be 'show' or 'save'!")
    if len(search_window) != 2:
        raise ValueError("search_window should have two values only - start and end times for burst searches")

    parent_folder = str(pathlib.Path(eventfile).parent)

    ev_header = fits.open(eventfile)[0].header
    MJDREFI = ev_header['MJDREFI']
    MJDREFF = ev_header['MJDREFF']
    source_name = ev_header['OBJECT']
    obsid = ev_header['OBS_ID']

    ### get the time series and zero-ise it
    #define an array of start times? So do in steps of dt from search_start to search_end
    times = fits.open(eventfile)[1].data['TIME']
    T_zeroized = times-times[0]
    counts = np.ones(len(T_zeroized))

    T_bins = np.linspace(0,np.ceil(T_zeroized[-1]),np.ceil(T_zeroized[-1])*1/tbin_size+1)
    binned_counts, bin_edges, binnumber = stats.binned_statistic(T_zeroized,counts,statistic='sum',bins=T_bins) #binning the photons

    for i in tqdm(range(len(T))): #for every window size:
        output_file = open(parent_folder + '/' + obsid + '_TBO_search_' + str(T[i]) + 's.txt','w')
        output_file.write('Source name: ' + source_name + ' ; ObsID: ' + obsid + '\n')
        output_file.write('Window size: T = ' + str(T[i]) + 's, stepping size = ' + str(dt[i]) + ' ; dt = ' + str(tbin_size) + '\n')

        T_start = np.arange(search_window[0],search_window[1],dt[i]) #start time of each sliding window
        T_end = T_start + T[i] #end time of each sliding window
        N = T[i]/tbin_size #number of trials for each window
        sig3 = Lv3_detection_level.power_for_sigma(3,N,1,1)
        sig4 = Lv3_detection_level.power_for_sigma(4,N,1,1)
        sig5 = Lv3_detection_level.power_for_sigma(5,N,1,1)

        output_file.write('Power needed for: 3 sigma - ' + str(sig3) + ' ; 4 sigma - ' + str(sig4) + ' ; 5 sigma - ' + str(sig5) + '\n')
        output_file.write('Starting/Ending MJD of TBO search scheme: ' + str(MJDREFI+MJDREFF+(times[0]+search_window[0])/86400) + '/' + str(MJDREFI+MJDREFF+(times[0]+search_window[1])/86400) + '\n')

        fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True) #dynamic power spectrum #define a 2x1 subplot or something
        fig.subplots_adjust(hspace=0)

        f_max = [] #corresponding frequencies to the maximum power
        ps_max = [] #to store the maximum power from each power spectrum of each sliding time series
        for j in tqdm(range(len(T_start))): #for every sliding window
            T_search = T_bins[:-1][(T_bins[:-1]>=T_start[j])&(T_bins[:-1]<=T_end[j])] #time series to search for burst oscillations
            binned_search = binned_counts[(T_bins[:-1]>=T_start[j])&(T_bins[:-1]<=T_end[j])]

            f,ps = Lv2_ps_method.manual(T_search,binned_search,[False,400,500],[False,400],False,[False,5]) #calculating Leahy-normalized power spectra
            f_window = f[(f>=f_central-df)&(f<=f_central+df)] #get frequency values 'df' Hz about f_central
            ps_window = ps[(f>=f_central-df)&(f<=f_central+df)] #get powers 'df' Hz about f_central

            scatt = ax1.scatter(x=np.ones(len(f_window))*T_start[j],y=f_window,s=12,c=ps_window,marker='o',cmap=cm.gist_heat,vmin=1,vmax=50)

            f_max.append(f_window[ps_window==np.max(ps_window)][0])
            ps_max.append(ps_window[ps_window==np.max(ps_window)][0])

            output_file.write('Start time for this window: zeroized - ' + str(T_start[j]) + ' ; MJD - ' + str(MJDREFI+MJDREFF+(times[0]+T_start[j])/86400) + '\n')
            for k in range(len(f_window)):
                output_file.write(str(f_window[k]) + ' ' + str(ps_window[k]) + ' ' + str(Lv3_detection_level.signal_significance(N,1,1,ps_window[k])) + '\n')

        output_file.close()

        if mode == "show":
            mplcursors.cursor(hover=True)
        ax1.set_title('Window size: ' + str(T[i]) + 's, dt='+str(dt[i])+'s \n' + 'Central freq. = '+str(f_central) + 'Hz, df = ' + str(df) + 'Hz \n Power required for 3 sigma: ' + str(sig3),fontsize=12)
        ax1.set_ylabel('Frequency (Hz)',fontsize=12)
        ax1.set_ylim([f_central-df,f_central+df])

        ax2.set_ylabel('Frequency (Hz)',fontsize=12)

        scat = ax2.scatter(x=T_start,y=f_max,s=12,c=ps_max,marker='o',cmap=cm.gist_heat,vmin=np.min(ps_max),vmax=np.max(ps_max),edgecolors='k')
        if mode == "show":
            mplcursors.cursor(hover=True)

        #fig.colorbar(scat,ax=ax1)
        #fig.colorbar(scat,ax=ax2)
        ax2.set_ylim([f_central-df,f_central+df])

        ps_contour = ax3.tricontour(T_start,f_max,ps_max,levels=30,linewidths=0.5,colors='k')
        ax3.clabel(ps_contour,fontsize=8)
        ax3.set_xlabel('Time (s)')
        ax3.set_ylim([f_central-df,f_central+df])

        if mode == "show":
            mplcursors.cursor(hover=True)
            plt.show()

        if mode == "save":
            filename = obsid + "_TBO_plots_" + str(T[i]) + 's.pdf'
            plt.savefig(parent_folder + '/' + filename,dpi=900)
            plt.close()

if __name__ == "__main__":
    eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/2584010501_pipe/ni2584010501_nicersoft_bary.evt'
    #burst_cands(eventfile)

    search_window = [5560,5660] #start/end times of burst oscillation search
    T = np.array([10]) #window sizes
    dt = T/10  #i.e., do overlapping windows of T seconds in steps of dt seconds
    tbin_size = 0.001 #size of time bins in seconds
    df = 10 #tolerance of 10 Hz, say
    f_central = 401 #central frequency
    mode = "show"
    dynamic_ps(eventfile,search_window,T,dt,tbin_size,df,f_central,mode)
