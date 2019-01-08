#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 17:30:42 2018

----------------------- !!! SPECTRAL ANALYSIS !!! -----------------------------

Columns of data: TIME, RAWX, RAWY, PHA, PHA_FAST, DET_ID, DEADTIME, EVENT_FLAGS,
TICK, MPU_A_TEMP, MPU_UNDER_COUNT, PI_FAST, PI, PI_RATIO

Recall that PI=energy/10 eV, pi=110 is 1.10keV

Spectral Analysis ONLY for Cen X-3! This is because of all the truncations
in the time series that I will be making.

003407010[1-4] - Cen X-3 [Commissioning phase] [20-23 June 2017]
103407010[1-6] - Cen X-3 [Science phase] [24-28 July 2017 and 7 Nov 2018]

Just to check my code for creating the hardness ratios:
https://arxiv.org/pdf/1806.02342.pdf?fbclid=IwAR0H78FWyi_sIyIWk0Z0oO09nyzs329xS9sJHiVBK1OHfNX2MFcA-POxkOg
"""
from __future__ import division
import datetime
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from scipy import stats
from scipy import signal

timestart = time.time()

print(datetime.datetime.now())

work_dir = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/'
obsid = '1034070104'  #Cen X-3
obsname = 'Cen_X-3'
#obsid = '1103010117' ## GRS 1915+105, from that paper by Neilsen et al. 2018
#obsname = 'GRS_1915+105'
t_interval = '0p01s'
extra = ''
place = 100 #for t_interval = 0p1s

def open_fits(work_dir, obsid):
    ## use barycentered data
    
    event = work_dir + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl_bary.evt'
    event = fits.open(event) #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return event

#def power_spectra(work_dir,obsid,t_interval) 
#    times, counts = np.loadtxt(work_dir + 'sumcounts_'+obsid+'_'+t_interval+'_bary.txt', usecols=(0,1), unpack=True) 

def get_data(work_dir,obsid):
    event = open_fits(work_dir,obsid)
    
    #PI gives energy info from 0.2keV-12keV, PI_FAST is supplementary info,
    #gives energy data about >1 keV
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


def get_cutoff(work_dir,obsid,doplot,manual,manual_cut):
    #given a SORTED array of counts vs energy (in keV), find the corresponding
    #E value such that you get 50% of counts on either side
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0]

    if manual == False: #if we let DATA choose the energy cut-off
        pi_data = np.array(pi_data,dtype=float)
    #    print pi_data[0:10]
    #    print type(pi_data[0:10])
    #    print type(pi_data[0])
        #sort the energies, then find the value which corresponds to half its length    
        sort_ind = np.argsort(pi_data) 
        sorted_E = pi_data[sort_ind]*10.0/1000.0 #into keV
        truncated_E = sorted_E[(sorted_E>=0.3)&(sorted_E<=12)]
        truncated_counts = np.ones(len(truncated_E))
        
        cutoff_index = int(len(truncated_E)/2)
        cutoff_E = truncated_E[cutoff_index]
        
        E_bins = np.linspace(0.3,12,1170+1)
        summed_counts, bin_edges, binnumber = stats.binned_statistic(truncated_E,truncated_counts,statistic='sum',bins=E_bins)
    
#    print sum(summed_counts[(E_bins[:-1]>=0.3)&(E_bins[:-1]<=4.3)])
#    print sum(truncated_counts[(truncated_E>=0.3)&(truncated_E<=4.3)])
#    print 
#    print sum(summed_counts[(E_bins[:-1]>4.3)&(E_bins[:-1]<=12)])
#    print sum(truncated_counts[(truncated_E>4.3)&(truncated_E<=12)])
    
    elif manual == True:
        cutoff_E = manual_cut
        cutoff_index = 0
        
    ### probably pointless to plot?
    if doplot == True:
        plt.figure(figsize=(10,8))
        plt.plot(E_bins[:-1],summed_counts,'rx')
        plt.xlabel('keV')
        plt.ylabel('Counts')   
        plt.axvline(x=cutoff_E,alpha=0.5,lw=0.5,color='k')
  
    return cutoff_E, cutoff_index ##CHANGE TO BOUNDARY ENERGY

#print get_cutoff(work_dir,obsid,True)


######################### GETTING THE ENERGIES ################################
def soft_counts(work_dir,obsid,cutoff,counts_data,truncated,truncated_pi):
    cutoff_PI = cutoff*1000/10 #get the cut-off in terms of PI; cutoff is in keV!
    #we do assume that there will be NO bins in which you get ZERO counts
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    if truncated == True:
        pi_data = truncated_pi
    np.place(counts_data,pi_data>=cutoff_PI,0) #get the counts for soft photons
    return counts_data

def hard_counts(work_dir,obsid,cutoff,counts_data,truncated,truncated_pi):
    cutoff_PI = cutoff*1000/10 #get the cut-off in terms of PI; cutoff is in keV!
    #we do assume that there will be NO bins in which you get ZERO counts
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    if truncated == True:
        pi_data = truncated_pi
    np.place(counts_data,pi_data<cutoff_PI,0) #get the counts for soft photons
    return counts_data

def soft_spec_lc(work_dir,obsid,manual,manual_cutoff,impose_xlim,xlim1,xlim2,impose_ylim,ylim1,ylim2):
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0]
    counts_data_all = np.ones(len(pi_data))
    cutoff_E = get_cutoff(work_dir,obsid,False,manual,manual_cutoff)[0]
    
    if impose_xlim == True:
        pi_data = pi_data[(shifted_t>=xlim1)&(shifted_t<=xlim2)]
        counts_data = counts_data_all[(shifted_t>=xlim1)&(shifted_t<=xlim2)]
        shifted_t = shifted_t[(shifted_t>=xlim1)&(shifted_t<=xlim2)]
    
    soft_c = soft_counts(work_dir,obsid,cutoff_E,counts_data,True,pi_data)
    pi_data = pi_data*10.0/1000.0 #to put into keV
#    sort_ind = np.argsort(pi_data)
#    sort_pi = pi_data[sort_ind]*10.0/1000.0 #to put into keV
#    soft_c = soft_c[sort_ind]
    
    t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1])+2)
    E_bins = np.linspace(0.3,12,117*2+1)
    
    summed_counts_lc, bin_edges, binnumber = stats.binned_statistic(shifted_t,soft_c,statistic='sum',bins=t_bins)
    summed_counts_spec, bin_edges, binnumber = stats.binned_statistic(pi_data,soft_c,statistic='sum',bins=E_bins)
    
    print max(summed_counts_lc[xlim1:xlim2])
    print t_bins[:-1][summed_counts_lc==max(summed_counts_lc[xlim1:xlim2])]
    
    checker_lc,lol1,lol2 = stats.binned_statistic(shifted_t,np.ones(len(soft_c)),statistic='sum',bins=t_bins)
    checker_spec,lol1,lol2 = stats.binned_statistic(pi_data,np.ones(len(soft_c)),statistic='sum',bins=E_bins)
    
    plt.figure(figsize=(10,8))
    plt.plot(t_bins[:-1],checker_lc,'b-')
    plt.xlim([11100,12000])
    plt.figure(figsize=(10,8))
    plt.plot(E_bins[:-1],checker_spec,'b-')
    
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(10,8)) 
    fig.suptitle('[Obs ID: ' + obsid+ ']; Object: ' + obsname + '\n Light curve + Spectra for soft photons \n for t='+str(xlim1)+'-'+str(xlim2)+'s; Boundary E = ' + str(cutoff_E) + ' keV', fontsize=15)
    if impose_xlim==False:
        fig.suptitle('[Obs ID: ' + obsid+ ']; Object: ' + obsname + '\n Light curve + Spectra for soft photons \n for t='+str(shifted_t[0])+'-'+str(shifted_t[-1])+'s; Boundary E = ' + str(cutoff_E) + ' keV', fontsize=15)
        
    ax1.plot(t_bins[:-1],summed_counts_lc,'r-')
    if impose_xlim==True:
        ax1.set_xlim([xlim1,xlim2])
        ax1.set_ylim([min(summed_counts_lc[xlim1:xlim2]),max(summed_counts_lc[xlim1:xlim2])])
    ax1.set_ylabel('Counts/s',fontsize=15)
    ax1.set_xlabel('Time',fontsize=15)
    
    ax2.plot(E_bins[:-1],summed_counts_spec,'r-')
#    if impose_xlim==True:
    ax2.set_xlim([0.3,cutoff_E])
    if impose_ylim==True:
        ax2.set_ylim([ylim1,ylim2])
    ax2.set_ylabel('Soft flux',fontsize=15)
    ax2.set_xlabel('Energy (keV)',fontsize=15)

    plt.subplots_adjust(hspace=0.2)
    
soft_spec_lc(work_dir,'1034070102',True,2.7,True,0,625,False,0,1)
    
def hard_spec_lc(work_dir,obsid,manual,manual_cutoff,impose_xlim,xlim1,xlim2,impose_ylim,ylim1,ylim2):
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0]
    counts_data_all = np.ones(len(pi_data))
    cutoff_E = get_cutoff(work_dir,obsid,False,manual,manual_cutoff)[0]
    
    if impose_xlim == True:
        pi_data = pi_data[(shifted_t>=xlim1)&(shifted_t<=xlim2)]
        counts_data = counts_data_all[(shifted_t>=xlim1)&(shifted_t<=xlim2)]
        shifted_t = shifted_t[(shifted_t>=xlim1)&(shifted_t<=xlim2)]
    
    hard_c = hard_counts(work_dir,obsid,cutoff_E,counts_data,True,pi_data)
    pi_data = pi_data*10.0/1000.0 #to put into keV
#    sort_ind = np.argsort(pi_data)
#    sort_pi = pi_data[sort_ind]*10.0/1000.0 #to put into keV
#    hard_c = hard_c[sort_ind]
    
    t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1])+2)   
    E_bins = np.linspace(0.3,12,117*2+1)
    summed_counts_lc, bin_edges, binnumber = stats.binned_statistic(shifted_t,hard_c,statistic='sum',bins=t_bins)
    summed_counts_spec, bin_edges, binnumber = stats.binned_statistic(pi_data,hard_c,statistic='sum',bins=E_bins)
    
#    print max(summed_counts_lc[xlim1:xlim2])
#    print t_bins[:-1][summed_counts_lc==max(summed_counts_lc[xlim1:xlim2])]
    
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(10,8)) 
    fig.suptitle('[Obs ID: ' + obsid+ ']; Object: ' + obsname + '\n Light curve + Spectra for hard photons \n for t='+str(xlim1)+'-'+str(xlim2)+'s; Boundary E = ' + str(cutoff_E) + ' keV', fontsize=15)
    if impose_xlim==False:
        fig.suptitle('[Obs ID: ' + obsid+ ']; Object: ' + obsname + '\n Light curve + Spectra for hard photons \n for t='+str(shifted_t[0])+'-'+str(shifted_t[-1])+'s; Boundary E = ' + str(cutoff_E) + ' keV', fontsize=15)
        
    ax1.plot(t_bins[:-1],summed_counts_lc,'r-')
    if impose_xlim==True:
        ax1.set_xlim([xlim1,xlim2])
        ax1.set_ylim([min(summed_counts_lc[xlim1:xlim2]),max(summed_counts_lc[xlim1:xlim2])])
    ax1.set_ylabel('Counts/s',fontsize=15)
    ax1.set_xlabel('Time (s)',fontsize=15)
    
    ax2.plot(E_bins[:-1],summed_counts_spec,'r-')
#    if impose_xlim==True:
    ax2.set_xlim([cutoff_E,12])
    if impose_ylim==True:
        ax2.set_ylim([ylim1,ylim2])
    ax2.set_ylabel('Hard flux',fontsize=15)
    ax2.set_xlabel('keV',fontsize=15)

    plt.subplots_adjust(hspace=0.2)
    
hard_spec_lc(work_dir,'1034070102',True,2.7,True,0,625,False,0,1)

def color(work_dir,obsid,manual,manual_cutoff):
    # cutoff is the energy cutoff for the hardness ratio/color
    # shifted_t is where the time series starts at t=0
    # pi_data is the corresponding energy tag for each event

    #for example, for the soft photons below some cutoff energy, replace each
    #element in the array corresponding to COUNTS by 0, IF the corresponding
    #energy value is ABOVE the cutoff (i.e., filter out higher E photons)
    
    ## getting the cutoff energy in keV
    cutoff,i = get_cutoff(work_dir,obsid,False,manual,manual_cutoff)
    
    ## getting the data
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0] 
    
    counts_data = np.ones(len(pi_data))
    counts_soft = soft_counts(work_dir,obsid,cutoff,counts_data,False,pi_data)
    ## re-obtain the data
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0]
    
    counts_data = np.ones(len(pi_data))
    counts_hard = hard_counts(work_dir,obsid,cutoff,counts_data,False,pi_data)

    t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1])+2) #get 1-second bins
    sum_soft, bin_edges_soft, binnumber_soft = stats.binned_statistic(shifted_t,counts_soft,statistic='sum',bins=t_bins)
    sum_hard, bin_edges_hard, binnumber_hard = stats.binned_statistic(shifted_t,counts_hard,statistic='sum',bins=t_bins)
    np.place(sum_soft,sum_soft==0,1) #so that you get 0/1 instead

    color = sum_hard/sum_soft
    color_diff = (sum_hard-sum_soft)/(sum_soft+sum_hard)
    
    return t_bins[:-1], color, color_diff

#color(work_dir,obsid,True,2.7)
    
def plotting(work_dir,obsid,manual,manual_cutoff,impose_xlim,impose_ylim,xlim1,xlim2,ylim1,ylim2,ylim1_diff,ylim2_diff):
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0] 
    cutoff_E = get_cutoff(work_dir,obsid,False,manual,manual_cutoff)[0]
    
    lightcurve_t, lightcurve_data = lightcurve(work_dir,obsid,False,False,0,0)
    color_t, color_data, color_diff = color(work_dir,obsid,manual,manual_cutoff)
 
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(10,8)) 
    fig.suptitle('[Obs ID: ' + obsid+ ']; Object: ' + obsname + '\n Light curve + Color (H/S) vs time + Color (H-S)/(H+S) vs time \n for t='+str(xlim1)+'-'+str(xlim2)+'s; Boundary E = ' + str(cutoff_E) + ' keV', fontsize=15)
    if impose_xlim==False:
        fig.suptitle('[Obs ID: ' + obsid+ ']; Object: ' + obsname + '\n Light curve + Color (H/S) vs time + Color (H-S)/(H+S) vs time \n for t='+str(shifted_t[0])+'-'+str(shifted_t[-1])+'s; Boundary E = ' + str(cutoff_E) + ' keV', fontsize=15)
        
    ax1.plot(lightcurve_t,lightcurve_data,'r-')
    if impose_xlim==True:
        ax1.set_xlim([xlim1,xlim2])
        ax1.set_ylim([min(lightcurve_data[xlim1:xlim2]),max(lightcurve_data[xlim1:xlim2])])
    ax1.set_ylabel('Counts/s',fontsize=15)
    
    ax2.plot(color_t,color_data,'r-')
    if impose_xlim==True:
        ax2.set_xlim([xlim1,xlim2])
    if impose_ylim==True:
        ax2.set_ylim([ylim1,ylim2])
    ax2.set_ylabel('Color \n (H/S)',fontsize=15)

    ax3.plot(color_t,color_diff,'r-')
    if impose_xlim==True:
        ax3.set_xlim([xlim1,xlim2])
    if impose_ylim==True:
        ax3.set_ylim([ylim1_diff,ylim2_diff])
    ax3.set_ylabel('Color \n (H-S)/(H+S)',fontsize=15)
    ax3.set_xlabel('Time (s)',fontsize=15)
    
    plt.subplots_adjust(hspace=0.1)
#    fig.tight_layout()
    
    return
    
#cutoff,i = get_cutoff(work_dir,obsid)
#soft_spec_lc(work_dir,obsid,True,2.7,True,11100,11950,False,0,1)
#hard_spec_lc(work_dir,obsid,True,2.7,True,11100,11950,False,0,1)
#plotting(work_dir,obsid,True,2.7,True,True,11100,11950,0,2,-0.8,0.6)

######################### GET INDICES OF INTERVALS ###########################
######## to get the indices of the boundaries of the different regions ########
#
#def get_subregion_indices(work_dir,obsid):#,count_cutoff):
#    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
#    shifted_t = times-times[0] 
#    
#    ## rebin the counts into 1-sec intervals
#    t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1]+2)) #get 1-second bins
#    pi_fast_clean = pi_fast
#    np.place(pi_fast_clean,pi_fast<0,0) #replace negative values by 0
#    intensity_1s, bin_edges, binnumber = stats.binned_statistic(shifted_t,pi_fast_clean,statistic='sum',bins=t_bins)
#
#    ### get indices where the counts are > 5:
#    index_subreg = [0]
#    threshold = np.where(intensity_1s>=5)
#    for i in range(len(threshold[0])-1):
#        if threshold[0][i+1]-threshold[0][i] > 50:
#            index_subreg.append(threshold[0][i])
#            index_subreg.append(threshold[0][i+1])
#    index_subreg.append(threshold[0][-1])
#    
#    return index_subreg

def get_subregion_indices(work_dir,obsid):
    ### FOR NOW, CONSIDER 1s BINS! 
    gtis = open_fits(work_dir,obsid)[2].data
    data_starts = np.array([gtis[i][0] for i in range(len(gtis))])
    data_stops = np.array([gtis[i][1] for i in range(len(gtis))])
    
    starts = data_starts-data_starts[0]
    stops = data_stops-data_starts[0]
    
    ### Use the list of GTIs to get subregion indices!
    index_subreg = [0]
    for i in range(len(starts)-1):
        if starts[i+1]-stops[i] >= 50:
            #if the time between start of next subregion and the end of the
            #previous subregion is > 50s, consider that a new subregion!
            index_subreg.append(int(round(stops[i])))
            index_subreg.append(int(round(starts[i+1])))
    index_subreg.append(int(round(stops[-1])))
    
    return index_subreg
#print get_subregion_indices(work_dir,obsid,590000) ##dependent on the observation!!
#indices = get_subregion_indices(work_dir,obsid)

############################### GETTING PDFs ##################################
#

def get_pdf(work_dir,obsid,manual,manual_cutoff,obsname,doplot,impose_xlim,impose_ylim,ylim1,ylim2,ylim1_diff,ylim2_diff):
    """ 
    Given the working directory, OBSID and object name, create a PDF which has 
    ...
    """
    
    from matplotlib.backends.backend_pdf import PdfPages
    
    filename = work_dir + obsid + '/' + obsid + '_' + obsname + '_v2.pdf'
    with PdfPages(filename) as pdf:
    #### WHOLE time series
        doplot = True
        impose_xlim = False
        lightcurve_t, lightcurve_data = lightcurve(work_dir,obsid,True,False,0,1) #but it's ok as impose_xlim is False
        pdf.savefig()
        plt.close()
        #    
        subr_i = get_subregion_indices(work_dir,obsid)
    #    #
        for i in range(0,len(subr_i),2): #can do a dictionary! 
    #    
            xlim1 = subr_i[i]
            xlim2 = subr_i[i+1]
            doplot = True
            impose_xlim = True
        
            lightcurve_t, lightcurve_data = lightcurve(work_dir,obsid,doplot,impose_xlim,xlim1,xlim2)
            pdf.savefig()
            plt.close()
            plotting(work_dir,obsid,manual,manual_cutoff,impose_xlim,False,xlim1,xlim2,ylim1,ylim2,ylim1_diff,ylim2_diff)
            pdf.savefig()
            plt.close()
            
            
            plotting(work_dir,obsid,manual,manual_cutoff,impose_xlim,True,xlim1,xlim2,ylim1,ylim2,ylim1_diff,ylim2_diff)
            pdf.savefig()
            plt.close()
             
ylim1 = 0
ylim2 = 2
ylim1_diff = -0.8
ylim2_diff = 0.8
doplot = True
impose_xlim = True
impose_ylim = True

#cenx3_obsid_list = ['0034070101','0034070102','0034070103','0034070104','1034070101','1034070102','1034070103','1034070104','1034070105','1034070106']
#ylims_list = [(0.6,1.2,-0.25,0.5),(0.5,1.5,-0.5,0.25),(0.6,1.5,-0.25,0.25),(0.5,1.5,-0.3,0.2),(0.6,1.5,-0.25,0.25),(0,2,-0.8,0.5),(0.5,1.5,-0.5,0.5),(0,100,0.9,1),(0,4,-1,0.5),(0,2.5,-0.5,0.3)]
#cenx3_obsid_list = ['1034070104']
#ylims_list = [(0,10,-0.5,0.9)]
####get_pdf(work_dir,'1034090111','GX 349+2',True,True,False,ylim1,ylim2,ylim1_diff,ylim2_diff)
###
#for i in range(len(cenx3_obsid_list)):
#    get_pdf(work_dir,cenx3_obsid_list[i],True,2.7,obsname,doplot,impose_xlim,impose_ylim,ylims_list[i][0],ylims_list[i][1],ylims_list[i][2],ylims_list[i][3])
#    print 'Done with ' + str(cenx3_obsid_list[i])    
        
timeend = time.time()

print (timeend-timestart)

## Old draft for get_subregion_indices. Can be improved!!

#        if subreg == 1:
#            if zeros_counter == 500: #i.e., if there are no counts for 500s. CAN CHANGE!!
#                subreg += 1 #got the indices for the 1st subregion
#                zeros_counter = 0
#                index_subreg.append(index_counts)
#            if intensity_1s[i] > 5.0:
#                index_counts = i
#                zeros_counter = 0
#            elif intensity_1s[i] <= 5.0:
#                zeros_counter += 1
#        else:
#            ### get starting position of next subregion
#            if intensity_1s[i] >= 5.0:
#                index_counts = i
#                zeros_counter = 0
#                index_subreg.append(index_counts)
#            
#            elif intensity_1s[i] <= 5.0:
#                zeros_counter += 1
#            if zeros_counter == 500: #i.e., if there are no counts for 500s. CAN CHANGE!!
#                index_subreg.append(index_counts)
#                if intensity_1s[i] > 5.0: #so after a 'long' time with no counts, if there are finally counts again, RECORD the index!
#                    index_subreg.append(i)

