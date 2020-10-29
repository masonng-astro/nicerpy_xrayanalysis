#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 2:59pm 2020

Simulating a time series similar to NGC 300 X-1, but with a larger S/N...

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from PyAstronomy.pyasl import foldAt
import Lv0_dirs,Lv2_dj_lsp,Lv2_swift_lc,Lv2_phase
from matplotlib.backends.backend_pdf import PdfPages
import os
from scipy import stats
from tqdm import tqdm
import subprocess
import pathlib
from stingray.pulse.pulsar import pulse_phase,phase_exposure,fold_events

from astropy.utils import iers
iers.conf.auto_download = False

def rebin_lc(corr_lc_files,corr_bg_files,bg_scale,tbin):
    rebinned_time = []
    rebinned_rate = []
    rebinned_errs = []
    rebinned_fracexp = []
    completeness = []

    times,rates,errors,fracexp = Lv2_swift_lc.get_bgsub(corr_lc_files,corr_bg_files,bg_scale)
    trunc_times = times-times[0]

    time_bins = np.arange(0,trunc_times[-1]+tbin,tbin)
    print('Rebinning...')
    for i in tqdm(range(len(time_bins)-1)):
        time_interval = trunc_times[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
        rate_interval = rates[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
        error_interval = errors[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
        fracexp_interval = fracexp[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]

        comp = len(time_interval)/(tbin/10)

        if len(time_interval) != 0:# and comp >= 0.99:
            mean_time = np.mean(time_interval)
            mean_rate = np.mean(rate_interval)
            mean_error = np.sqrt(np.sum(error_interval**2))/np.size(error_interval)
            sum_fracexp = sum(fracexp_interval)

            rebinned_time.append(mean_time)
            rebinned_rate.append(mean_rate)
            rebinned_errs.append(mean_error)
            rebinned_fracexp.append(sum_fracexp)
            completeness.append(len(time_interval)/(tbin/10))

    return np.array(rebinned_time),np.array(rebinned_rate),np.array(rebinned_errs),np.array(rebinned_fracexp)

def rebin_txt(times,rates,errors,fracexp,tbin,cmpltness):
    """
    Very similar to rebin_lc, but takes in the 10s-binned light curve from running
    rebin_lc with tbin=10s. This saves time!
    """

    rebinned_time = []
    rebinned_rate = []
    rebinned_errs = []
    rebinned_fracexp = []
    completeness = []

    trunc_times = times-times[0]

    time_bins = np.arange(0,trunc_times[-1]+tbin,tbin)
    print('Rebinning...')
    for i in tqdm(range(len(time_bins)-1)):
        time_interval = trunc_times[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
        rate_interval = rates[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
        error_interval = errors[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]
        fracexp_interval = fracexp[(trunc_times>=time_bins[i])&(trunc_times<time_bins[i+1])]

        comp = len(time_interval)/(tbin/10)

        if len(time_interval) != 0 and comp >= cmpltness:
            mean_time = np.mean(time_interval)
            mean_rate = np.mean(rate_interval)
            mean_error = np.sqrt(np.sum(error_interval**2))/np.size(error_interval)
            sum_fracexp = sum(fracexp_interval)

            rebinned_time.append(mean_time)
            rebinned_rate.append(mean_rate)
            rebinned_errs.append(mean_error)
            rebinned_fracexp.append(sum_fracexp)
            completeness.append(len(time_interval)/(tbin/10))

    return np.array(rebinned_time),np.array(rebinned_rate),np.array(rebinned_errs),np.array(rebinned_fracexp)

def phase_folding(t,y,T,T0,f,nbins):
    """
    Calculating the folded profile
    Goes from 0 to 2.

    x - array of time values
    y - flux array
    T - sum of all the GTIs
    T0 - reference epoch in MJD
    f - folding frequency
    nbins - number of phase bins desired
    """
    MJDREFI = 51910
    MJDREFF = 7.428703700000000E-04
    TIMEZERO = 0

    t_MJDs =  MJDREFI + MJDREFF + (TIMEZERO+t)/86400

    tau = (t_MJDs-T0)*86400
    #phase = (f*tau + fdot/2 *tau**2 + fdotdot/6*tau**3)%1
    phase = (f*tau)%1

    phase_bins = np.linspace(0,1,nbins+1)
    summed_profile,bin_edges,binnumber = stats.binned_statistic(phase,y,statistic='mean',bins=phase_bins)
    error = np.sqrt(summed_profile*100)/100

    phase_bins_total = np.array(list(phase_bins[:-1]) + list(phase_bins+1))
    summed_profile_total = np.array(list(summed_profile)*2)
    error_total = np.array(list(error)*2) #the 100 should change!!!

    return phase_bins_total, summed_profile_total, error_total

##### Part 1 of the simulation

eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge_niceroverlap_all.evt'
times = fits.open(eventfile)[1].data['TIME']
gtis_data = fits.open(eventfile)[2].data

gtis_conform = []
for i in range(len(gtis_data)):
    gtis_conform.append([gtis_data[i][0],gtis_data[i][1]]) #conform to the input that Stingray uses

T = sum([ gtis_data[i]['STOP']-gtis_data[i]['START'] for i in range(len(gtis_data)) ]) #exposure time
x = times-times[0]
pb = 120e3 #120ks orbital period
longx = np.linspace(x[0],x[-1],10001)

y = 5*np.sin(2*np.pi/pb * x) + 5 #+ np.random.normal(0,10,size=len(x))
longy = 5*np.sin(2*np.pi/pb * longx) + 5# + np.random.normal(0,1,size=len(longx))

"""
omega,psd,prob3,prob4,prob5 = Lv2_dj_lsp.lsp(longx,longy)
freq = omega/(2*np.pi)

print(psd[psd>=0.98*np.max(psd)],freq[psd>=0.98*np.max(psd)])

plt.figure()
plt.plot(freq,psd,'rx-')
plt.xlabel('Frequency (Hz)',fontsize=12)
plt.ylabel('Normalized Power',fontsize=12)
plt.axhline(y=prob3,lw=0.5,alpha=0.5)
plt.axhline(y=prob4,lw=0.5,alpha=0.5)
plt.axhline(y=prob5,lw=0.5,alpha=0.5)
plt.show()
"""

nphase = 20

##### using foldAt
phases = foldAt(x,pb,T0=0)

##### using phase mod 1
phase_mod = (1/pb * x)%1 #or shift by -0.7*pb

##### using stingray.pulse.pulsar
phase_stingray = pulse_phase(x,[1/pb])#,ph0=1-0.7)

expocorr = Lv2_phase.phase_exposure(times[0]-times[0],times[-1]-times[0],period=pb,nbin=nphase,gtis=gtis_conform)

################################################################################

##### Testing the 3 different routines for calculating phase

phase_bins = np.linspace(0,1,nphase+1)
profile,bin_edges,binnumber = stats.binned_statistic(phases,y,statistic='mean',bins=nphase)
profile_mod,bin_edges,binnumber = stats.binned_statistic(phase_mod,y,statistic='mean',bins=nphase)
profile_sr,bin_edges,binnumber = stats.binned_statistic(phase_stingray,y,statistic='mean',bins=nphase)

phase_to_2 = np.array(list(phase_bins[:-1]) + list(phase_bins+1))
profile_to_2 = np.array(list(profile)*2)
profile_mod_to_2 = np.array(list(profile_mod)*2)
profile_sr_to_2 = np.array(list(profile_sr)*2)

plt.step(phase_to_2[:-1],profile_to_2,'b-')
plt.step(phase_to_2[:-1],profile_mod_to_2,'r-')
plt.step(phase_to_2[:-1],profile_sr_to_2,'k-')
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Flux',fontsize=12)
plt.legend(('foldAt','f*t mod 1','stingray'),fontsize=12,loc='best')
#plt.show()

##### Doing phase shifts...

offset = 0.7*nphase

##### Shifting pulse profiles through a shifted FT (see Deepto's 7/20/2020 email)
if nphase % 2 == 0:
    fft_x = np.array(list(np.arange(int(nphase/2)+1)) + list(np.arange(int(nphase/2)-1) - (int(nphase/2)-1)))
else:
    fft_x = np.array(list(np.arange(int(nphase/2)+1)) + list(np.arange(int(nphase/2)) - int(nphase/2)))

shift = np.exp(-2j*np.pi*fft_x*offset/nphase)
shifted_prof = np.real(np.fft.ifft(np.fft.fft(profile_mod)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile
#shifted_err_sr = np.real(np.fft.ifft(np.fft.fft(err_sr)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile

shifted_prof_to_2 = np.array(list(shifted_prof)*2)

plt.figure()
plt.step(phase_to_2[:-1],profile_mod_to_2,'r-')
plt.errorbar(phase_to_2[:-1],shifted_prof_to_2,color='c',drawstyle='steps-mid')
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)

plt.show()

##### Part 2 of the simulation
bary_outputfolder = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/lightcurve/'
obsids = [str(i) for i in range(49834027,49834042)] + [str(i) for i in range(49834043,49834062)] + [str(i) for i in range(49834063,49834066)] + ['88810002'] + [str(i) for i in range(49834066,49834069)] + [str(i) for i in range(49834070,49834079)] + [str(i) for i in range(49834080,49834088)]
corr_lc_files = [bary_outputfolder + 'sw000' + obsids[i] + '_corr.lc' for i in range(len(obsids))]
corr_bg_files = [bary_outputfolder + 'sw000' + obsids[i] + '_bg_corr.lc' for i in range(len(obsids))]
bg_scale_x1 = (30/120)**2

"""
times,rates,errors,fracexp = Lv2_swift_lc.get_bgsub(corr_lc_files,corr_bg_files,bg_scale_x1)

simfile = open(bary_outputfolder + 'simulate_10s.txt','w')
print('Writing into the text file...')
for i in tqdm(range(len(times))):
    simfile.write(str(times[i]-times[0]) + ' ' + str(rates[i]) + ' ' + str(errors[i]) + ' ' + str(fracexp[i]) +  '\n')
simfile.close()
"""

"""
##### Calling the text file (faster than always running "rebin_lc")
txt_t,txt_rate,txt_err,txt_fracexp = np.genfromtxt(bary_outputfolder + 'simulate_10s.txt',usecols=(0,1,2,3),unpack=True)
rebinned_t,rebinned_rate,rebinned_err,rebinned_fracexp = rebin_txt(txt_t,txt_rate,txt_err,txt_fracexp,100,0)
#new_rebinned_rate = np.array([ 5*np.sin(2*np.pi/pb * rebinned_t[i]) + 5 + np.random.normal(0,5) if rebinned_rate[i] > 0 else 0 for i in range(len(rebinned_t)) ])

#for i in range(50):
#    print(rebinned_t[i],rebinned_rate[i],new_rebinned_rate[i])

print(len(rebinned_t),len(rebinned_rate))
##### CHECKING THE PERIODOGRAM - IT IS THE SAME!
omega,psd,prob3,prob4,prob5 = Lv2_dj_lsp.lsp(rebinned_t,rebinned_rate)
#omega_new,psd_new,prob3_new,prob4_new,prob5_new = Lv2_dj_lsp.lsp(rebinned_t,new_rebinned_rate)

freq = omega/(2*np.pi)
#freq_new = omega_new/(2*np.pi)

plt.figure()
plt.plot(freq,psd,'rx-')
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel('Frequency (Hz)',fontsize=12)
plt.ylabel('Normalized Power',fontsize=12)
plt.axhline(y=prob3,lw=0.5,alpha=0.5)
plt.axhline(y=prob4,lw=0.5,alpha=0.5)
plt.axhline(y=prob5,lw=0.5,alpha=0.5)
#print(prob3,prob4,prob5)
print(freq[psd==np.max(psd)][0],psd[psd==np.max(psd)][0])

#plt.figure()
#plt.plot(freq_new,psd_new,'rx-')
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel('Frequency (Hz)',fontsize=12)
#plt.ylabel('Normalized Power',fontsize=12)
#plt.axhline(y=prob3_new,lw=0.5,alpha=0.5)
#plt.axhline(y=prob4_new,lw=0.5,alpha=0.5)
#plt.axhline(y=prob5_new,lw=0.5,alpha=0.5)
#print(prob3,prob4,prob5)
#print(freq_new[psd_new==np.max(psd_new)][0],psd_new[psd_new==np.max(psd_new)][0])

plt.show()
"""


##### Testing the folding...
tstart_49834027 = 546830295.758713
tstart_49834027_MJD = fits.open(eventfile)[1].header['MJDREFI'] + fits.open(eventfile)[1].header['MJDREFF'] + tstart_49834027/86400

"""
phases,profile,error = phase_folding(rebinned_t+tstart_49834027,rebinned_rate,T,tstart_49834027_MJD,8.46465218853785e-06,nphase)
expocorr = Lv2_phase.phase_exposure(times[0]-times[0],times[-1]-times[0],period=1/8.46465218853785e-06,nbin=nphase,gtis=gtis_conform-times[0])
plt.figure()
plt.errorbar(x=phases[:-1],y=profile,yerr=error,color='r',drawstyle='steps-mid')
plt.errorbar(x=phases[:-1],y=profile/np.array(list(expocorr)*2),yerr=error/np.array(list(expocorr)*2),color='b',drawstyle='steps-mid')
plt.legend(('Folded','Expo-corr'),fontsize=12)

print(expocorr)

phases_new,profile_new,error_new = phase_folding(rebinned_t+tstart_49834027,new_rebinned_rate,T,tstart_49834027_MJD,8.334549348081744e-06,nphase)
expocorr_new = Lv2_phase.phase_exposure(times[0]-times[0],times[-1]-times[0],period=1/8.334549348081744e-06,nbin=nphase,gtis=gtis_conform-times[0])
plt.figure()
plt.errorbar(x=phases_new[:-1],y=profile_new,yerr=error_new,color='r',drawstyle='steps-mid')
plt.errorbar(x=phases_new[:-1],y=profile_new/np.array(list(expocorr_new)*2),yerr=error_new/np.array(list(expocorr_new)*2),color='b',drawstyle='steps-mid')
plt.legend(('Folded','Expo-corr'),fontsize=12)

print(expocorr_new)

plt.show()
"""

"""
##### Doing the chi^2 exploration

chi2 = []
freqs = np.arange(8.3e-6,8.4e-6,0.0001e-6)
for i in tqdm(range(len(freqs))):
    phases,profile,error = phase_folding(rebinned_t+tstart_49834027,new_rebinned_rate,T,tstart_49834027_MJD,freqs[i],nphase)
    expocorr = Lv2_phase.phase_exposure(times[0]-times[0],times[-1]-times[0],period=1/freqs[i],nbin=nphase,gtis=gtis_conform-times[0])
    chi2.append( Lv2_phase.get_chi2(profile/np.array(list(expocorr)*2),error/np.array(list(expocorr)*2) ) )

plt.figure()
plt.plot(freqs,chi2,'rx-')
plt.yscale('log')
plt.xlabel('Frequency (Hz)',fontsize=12)
plt.ylabel('chi^2 [ sum( (profile-mean)^2/error^2) ]',fontsize=12)
plt.show()
"""

"""
print('Looking at completeness...')
completeness = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
txt_t,txt_rate,txt_err,txt_fracexp = np.genfromtxt(bary_outputfolder + 'simulate_10s.txt',usecols=(0,1,2,3),unpack=True)
noise_rate = np.array([ 0.005*np.sin(2*np.pi/pb * txt_t[i]) + 0.005 + np.random.normal(0,0.005) if txt_rate[i] > 0 else 0 for i in range(len(txt_t)) ])
noise_err = np.array([ 0 if noise_rate[i] == 0 else np.sqrt(np.abs(noise_rate[i])*10)/10 for i in range(len(noise_rate)) ])

for i in range(len(completeness)):
    rebinned_t,rebinned_rate,rebinned_err,rebinned_fracexp = rebin_txt(txt_t,noise_rate,noise_err,txt_fracexp,100,completeness[i])
    omega_new,psd_new,prob3_new,prob4_new,prob5_new = Lv2_dj_lsp.lsp(rebinned_t,rebinned_rate)
    freq_new = omega_new/(2*np.pi)

    freqs_list, psd_list = Lv2_dj_lsp.psd_error(rebinned_t,rebinned_rate,rebinned_err)
    print(str(completeness[i]*100) + '%')
    print(freq_new[psd_new==np.max(psd_new)][0],psd_new[psd_new==np.max(psd_new)][0])
    print(len(rebinned_t))
    print('Median frequency: ' + str(np.median(freqs_list)))
    print('Error in frequency: ' + str(np.std(freqs_list)))
"""
