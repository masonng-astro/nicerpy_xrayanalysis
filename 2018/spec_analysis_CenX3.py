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
obsid = '0034070101'  #Cen X-3
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
    pi_fast_clean = pi_fast
    np.place(pi_fast_clean,pi_fast<0,0) #replace negative values by 0
    
    shifted_t = times-times[0]
        
    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,pi_fast_clean,statistic='sum',bins=t_bins)
    
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


def get_cutoff(work_dir,obsid,doplot):
    #given a SORTED array of counts vs energy (in keV), find the corresponding
    #E value such that you get 50% of counts on either side
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    pi_fast_clean = pi_fast
    np.place(pi_fast_clean,pi_fast<0,0) #replace negative values by 0

    sort_ind = np.argsort(pi_data) 
    sorted_E = pi_data[sort_ind]*10/1000 #into keV
    sorted_counts = pi_fast_clean[sort_ind]
    
    half_counts = int(0.5*sum(sorted_counts)) 
    tally = 0
    for i in tqdm(range(len(sorted_counts))):
        if tally < half_counts:
            tally += sorted_counts[i]
            energy = sorted_E[i]
        else:
            break
  
    ### probably pointless to plot?
    if doplot == True:
        plt.figure(figsize=(10,8))
        plt.plot(sorted_E,sorted_counts,'rx')
        plt.xlabel('keV')
        plt.ylabel('Counts')   
        plt.axvline(x=cutoff,alpha=0.5,lw=0.5,color='k')
  
    return energy, i

######################### GETTING THE ENERGIES ################################
def soft_counts(cutoff,shifted_t,pi_data,pi_fast_clean_soft):
    cutoff_PI = cutoff*1000/10 #get the cut-off in terms of PI; cutoff is in keV!

    #we do assume that there will be NO bins in which you get ZERO counts
    np.place(pi_fast_clean_soft,pi_data>=cutoff_PI,0) #get the counts for soft photons

    return pi_fast_clean_soft

def hard_counts(cutoff,shifted_t,pi_data,pi_fast_clean_hard):
    cutoff_PI = cutoff*1000/10 #get the cut-off in terms of PI; cutoff is in keV!

    #we do assume that there will be NO bins in which you get ZERO counts
    np.place(pi_fast_clean_hard,pi_data<cutoff_PI,0) #get the counts for soft photons

    return pi_fast_clean_hard

def color(work_dir,obsid):
    # cutoff is the energy cutoff for the hardness ratio/color
    # shifted_t is where the time series starts at t=0
    # pi_data is the corresponding energy tag for each event
    # pi_fast_noneg is the corresponding counts 

    #for example, for the soft photons below some cutoff energy, replace each
    #element in the array corresponding to COUNTS by 0, IF the corresponding
    #energy value is ABOVE the cutoff (i.e., filter out higher E photons)
    
    ## getting the cutoff energy in keV
    cutoff,i = get_cutoff(work_dir,obsid,False)
    
    ## getting the data
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0] 
    
    pi_fast_clean = pi_fast
    np.place(pi_fast_clean,pi_fast<0,0) #replace negative values by 0
    pi_fast_soft = soft_counts(cutoff,shifted_t,pi_data,pi_fast_clean)

    ## re-obtain the data
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0]
    
    pi_fast_clean = pi_fast
    np.place(pi_fast_clean,pi_fast<0,0) #replace negative values by 0
    pi_fast_hard = hard_counts(cutoff,shifted_t,pi_data,pi_fast_clean)

    t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1]+2)) #get 1-second bins
    sum_soft, bin_edges_soft, binnumber_soft = stats.binned_statistic(shifted_t,pi_fast_soft,statistic='sum',bins=t_bins)
    sum_hard, bin_edges_hard, binnumber_hard = stats.binned_statistic(shifted_t,pi_fast_hard,statistic='sum',bins=t_bins)
    np.place(sum_soft,sum_soft==0,1) #so that you get 0/1 instead

    color_diff = (sum_hard-sum_soft)/(sum_soft+sum_hard)
    color = sum_hard/sum_soft
    
    return t_bins[:-1], color, color_diff
    
def plotting(work_dir,obsid,impose_xlim,impose_ylim,xlim1,xlim2,ylim1,ylim2,ylim1_diff,ylim2_diff):
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0] 
    cutoff_E = get_cutoff(work_dir,obsid,False)[0]
    
    lightcurve_t, lightcurve_data = lightcurve(work_dir,obsid,False,False,0,0)
    color_t, color_data, color_diff = color(work_dir,obsid)
 
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(10,8)) 
    fig.suptitle('[Obs ID: ' + obsid+ ']; Object: ' + obsname + '\n Light curve + Color (H/S) vs time + Color (H-S)/(H+S) vs time \n for t='+str(xlim1)+'-'+str(xlim2)+'s; Cut-off E = ' + str(cutoff_E) + ' keV', fontsize=15)
    if impose_xlim==False:
        fig.suptitle('[Obs ID: ' + obsid+ ']; Object: ' + obsname + '\n Light curve + Color (H/S) vs time + Color (H-S)/(H+S) vs time \n for t='+str(shifted_t[0])+'-'+str(shifted_t[-1])+'s; Cut-off E = ' + str(cutoff_E) + ' keV', fontsize=15)
        
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
    
    plt.subplots_adjust(hspace=0.05)
#    fig.tight_layout()
    
    return
    
#cutoff,i = get_cutoff(work_dir,obsid)
#plotting(work_dir,obsid,True,True,5000,5830,0.6,1.3,-0.3,0.2)

######################### GET INDICES OF INTERVALS ###########################
######## to get the indices of the boundaries of the different regions ########

def get_subregion_indices(work_dir,obsid):#,count_cutoff):
    pi_fast, times, pi_data, pi_ratio, flags = get_data(work_dir,obsid)
    shifted_t = times-times[0] 
    
    ## rebin the counts into 1-sec intervals
    t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1]+2)) #get 1-second bins
    pi_fast_clean = pi_fast
    np.place(pi_fast_clean,pi_fast<0,0) #replace negative values by 0
    intensity_1s, bin_edges, binnumber = stats.binned_statistic(shifted_t,pi_fast_clean,statistic='sum',bins=t_bins)

    ### get indices where the counts are > 5:
    index_subreg = [0]
    threshold = np.where(intensity_1s>=5)
    for i in range(len(threshold[0])-1):
        if threshold[0][i+1]-threshold[0][i] > 50:
            index_subreg.append(threshold[0][i])
            index_subreg.append(threshold[0][i+1])
    index_subreg.append(threshold[0][-1])
    
    return index_subreg
    
#print get_subregion_indices(work_dir,obsid,590000) ##dependent on the observation!!
indices = get_subregion_indices(work_dir,obsid)

############################### GETTING PDFs ##################################
#

def get_pdf(work_dir,obsid,obsname,doplot,impose_xlim,impose_ylim,ylim1,ylim2,ylim1_diff,ylim2_diff):
    """ 
    Given the working directory, OBSID and object name, create a PDF which has 
    ...
    """
    
    from matplotlib.backends.backend_pdf import PdfPages
    
    filename = work_dir + obsid + '/' + obsid + '_' + obsname + '.pdf'
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
            plotting(work_dir,obsid,impose_xlim,False,xlim1,xlim2,ylim1,ylim2,ylim1_diff,ylim2_diff)
            pdf.savefig()
            plt.close()
            
            
            plotting(work_dir,obsid,impose_xlim,True,xlim1,xlim2,ylim1,ylim2,ylim1_diff,ylim2_diff)
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
#ylims_list = [(0.6,1.2,-0.25,0.5),(0.5,1.5,-0.5,0.25),(0.6,1.5,-0.25,0.25),(0.5,1.5,-0.3,0.2),(0.6,1.5,-0.25,0.25),(0,2,-0.8,0.5),(0.5,1.5,-0.5,0.5),(0,3,-1,0.5),(0,4,-1,0.5),(0,2.5,-0.5,0.3)]

#get_pdf(work_dir,'1034090111','GX 349+2',True,True,False,ylim1,ylim2,ylim1_diff,ylim2_diff)

cenx3_obsid_list = ['1034070104']
ylims_list = [(0,3,-1,0.5)]
for i in range(len(cenx3_obsid_list)):
    get_pdf(work_dir,cenx3_obsid_list[i],obsname,doplot,impose_xlim,impose_ylim,ylims_list[i][0],ylims_list[i][1],ylims_list[i][2],ylims_list[i][3])
    print 'Done with ' + str(cenx3_obsid_list[i])    
        
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

