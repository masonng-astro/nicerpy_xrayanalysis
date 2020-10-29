#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 7:07pm 2020

Meant to interface with Lv2_dj_lsp and functions from stingray.pulse.pulsar to
analyze Swift data pertaining to NGC 300 X-1 in one place, instead of having
the analysis spread between Lv2_dj_lsp.py and test.py

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
from matplotlib import cm
from PyAstronomy.pyasl import foldAt
from mpl_toolkits.mplot3d import Axes3D
import mplcursors
import pathlib
from stingray.pulse.pulsar import pulse_phase,phase_exposure,fold_events

#####
## Noting here first that all the barycentering, time-ordering, extracting events
## (with XSELECT), doing exposure corrections (xrtlccorr), and subsequently the
## background subtraction, are all done in Lv2_swift_lc. There's no need to do so here.
#####

##### Parameters
eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge_niceroverlap_all.evt' #1 year of data; overlaps with NICER
#eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_swift_dec16_may19.evt'
#eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge.evt' #all 14 years
eventfile_xmm = '/Volumes/Samsung_T5/NGC300_XMMdata/ngc300x1_pn.evt'

times = fits.open(eventfile)[1].data['TIME'] #getting array of times
times_xmm = fits.open(eventfile_xmm)[1].data['TIME']

gtis_data = fits.open(eventfile)[2].data #getting GTIs
gtis_data_xmm = fits.open(eventfile_xmm)[59].data #59 for pn, 15 for mos1, 19 for mos2

T = sum([ gtis_data[i]['STOP']-gtis_data[i]['START'] for i in range(len(gtis_data)) ]) #exposure time
T_xmm = sum([ gtis_data_xmm[i]['STOP']-gtis_data_xmm[i]['START'] for i in range(len(gtis_data_xmm)) ]) #exposure time
print(T_xmm)
T0_MJD = fits.open(eventfile)[1].header['MJDREFI'] + fits.open(eventfile)[1].header['MJDREFF'] + fits.open(eventfile)[1].header['TSTART']/86400 #SWIFT
T0_MJD_eclipse = 58239.3498 #mid-eclipse!
T0_MJD_xmm = fits.open(eventfile_xmm)[1].header['MJDREF'] + fits.open(eventfile_xmm)[1].header['TSTART']/86400 #XMM-NEWTON

MJDREFI = fits.open(eventfile)[1].header['MJDREFI'] #Swift
MJDREFF = fits.open(eventfile)[1].header['MJDREFF'] #Swift
MJDREF = fits.open(eventfile_xmm)[1].header['MJDREF'] #XMM-Newton
diff_swiftxmm = (MJDREFI+MJDREFF-MJDREF)*86400

##### Get the phase offset between Swift eclipse time and XMM's first event time:
Porb_days = (1/8.4712e-6)/86400
xmm_first = MJDREF + times_xmm[0]/86400
no_cycles = (T0_MJD_eclipse - T0_MJD_xmm)/Porb_days
xmm_ecl = T0_MJD_eclipse - int(no_cycles)*Porb_days #time of the mid-eclipse BEFORE the first XMM event
if xmm_ecl > xmm_first:
    xmm_ecl -= Porb_days
phaseoff = (xmm_first-xmm_ecl)/Porb_days
print('Phase offset is ' + str(phaseoff))

##### Be careful here, as Swift and XMM have different MJDREFs!!!
gtis_conform = []
for i in range(len(gtis_data)):
    gtis_conform.append([gtis_data[i][0],gtis_data[i][1]]) #conform to the input that Stingray uses

gtis_conform_xmm = []
for i in range(len(gtis_data_xmm)):
    gtis_conform_xmm.append([gtis_data_xmm[i][0],gtis_data_xmm[i][1]]) #conform to the input that Stingray uses

#bary_outputfolder = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/lightcurve/'
#obsids = [str(i) for i in range(49834027,49834042)] + [str(i) for i in range(49834043,49834062)] + [str(i) for i in range(49834063,49834066)] + ['88810002'] + [str(i) for i in range(49834066,49834069)] + [str(i) for i in range(49834070,49834079)] + [str(i) for i in range(49834080,49834088)]
#corr_lc_files = [bary_outputfolder + 'sw000' + obsids[i] + '_corr.lc' for i in range(len(obsids))]
#corr_ulx1_files = [bary_outputfolder + 'sw000' + obsids[i] + '_ulx1_corr.lc' for i in range(len(obsids))]
#corr_bg_files = [bary_outputfolder + 'sw000' + obsids[i] + '_bg_corr.lc' for i in range(len(obsids))]
#bg_scale_x1 = (30/120)**2
#bg_scale_ulx1 = (35/120)**2
#completeness = np.array([0,10,20,30,40,50,60,70,80,90,100])/100
#rebinned_t, rebinned_rate, rebinned_err, rebinned_fracexp = Lv2_dj_lsp.rebin_lc(corr_lc_files,corr_bg_files,bg_scale_x1,100,0.5)
#rebinned_t_ulx1, rebinned_rate_ulx1, rebinned_err_ulx1, rebinned_fracexp_ulx1 = rebin_lc(corr_ulx1_files,corr_bg_files,bg_scale_ulx1,3600,0)

#tstart_49834027 = 546830295.758713

"""
### Writing the data from the light curves of X-1 and ULX-1 into text files; also plotting the light curve, This is mainly for 3600s bins

x1_text = open(bary_outputfolder + 'ngc300x1_bg_exp_corr_lc_3600s.txt','w')
ulx1_text = open(bary_outputfolder + 'ngc300ulx1_bg_exp_corr_lc_3600s.txt','w')

for i in range(len(rebinned_t)):
    x1_text.write(str(51910 + 7.428703700000000E-04+(rebinned_t[i]+tstart_49834027)/86400) + ' ' + str(rebinned_rate[i]) + ' ' + str(rebinned_err[i]) + '\n')
x1_text.close()

for i in range(len(rebinned_t_ulx1)):
    ulx1_text.write(str(51910 + 7.428703700000000E-04 + (rebinned_t_ulx1[i]+tstart_49834027)/86400) + ' ' + str(rebinned_rate_ulx1[i]) + ' ' + str(rebinned_err_ulx1[i]) + '\n')
ulx1_text.close()

mjd = 51910 + 7.428703700000000E-04 + (tstart_49834027+rebinned_t)/86400
mjd_ulx1 = 51910 + 7.428703700000000E-04 + (tstart_49834027+rebinned_t_ulx1)/86400
plt.errorbar(x=mjd[rebinned_err<=0.06],y=rebinned_rate[rebinned_err<=0.06],yerr=rebinned_err[rebinned_err<=0.06],fmt='rx')
plt.errorbar(x=mjd_ulx1[rebinned_err_ulx1<=0.06],y=rebinned_rate_ulx1[rebinned_err_ulx1<=0.06],yerr=rebinned_err_ulx1[rebinned_err_ulx1<=0.06],fmt='bx')
plt.legend(('X-1','ULX-1'),fontsize=12)
plt.xlabel('Time (MJD)',fontsize=12)
plt.ylabel('[Exposure-corrected] Count rate (c/s)',fontsize=12)
plt.axhline(y=0,color='k',lw=0.5,alpha=0.5)
plt.show()
"""

### Running Lv2_dj_lsp.lsp
"""
for i in range(len(completeness)):
    rebinned_t, rebinned_rate, rebinned_err, rebinned_fracexp = Lv2_dj_lsp.rebin_lc(corr_lc_files,corr_bg_files,bg_scale_x1,100,completeness[i])

    omega,psd,prob3,prob4,prob5 = Lv2_dj_lsp.lsp(rebinned_t,rebinned_rate)
    nu_reg = omega/(2.0*np.pi)
    freq = omega/(2*np.pi)

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

    print(np.max(psd),freq[psd==np.max(psd)][0])
    #plt.show()
"""

### Doing FR/RSS

#for i in range(len(completeness)):
#    rebinned_t, rebinned_rate, rebinned_err, rebinned_fracexp = Lv2_dj_lsp.rebin_lc(corr_lc_files,corr_bg_files,bg_scale_x1,100,completeness[i])
#    freqs_list, psd_list = Lv2_dj_lsp.psd_error(rebinned_t,rebinned_rate,rebinned_err)
#    print(str(completeness[i]) + '%')
#    print('Median frequency: ' + str(np.median(freqs_list)))
#    print('Error in frequency: ' + str(np.std(freqs_list)))
#print('Powers: ' + str(psd_list))

################################################################################
################################### FOLDING ####################################
################################################################################

"""
##### Folding using my routine; confirmed that the folding of the raw data agrees with Stingray's and foldAt
nbins = 20
freq = 8.4712e-6
offset = -0.215*nbins
#freq = 8.6088e-6
freqdot = 0
freqdotdot = 0
phase_frac = (T0_MJD_eclipse-T0_MJD)/((1/freq)/86400)

#print('MID ECLIPSE TIME:')
#print(  fits.open(eventfile)[1].header['MJDREFI'] + fits.open(eventfile)[1].header['MJDREFF'] + (times[0] + 0.21569724*1/freq)/86400)
#T0_MJD = fits.open(eventfile)[1].header['MJDREF'] + times[0]/86400

##### Using Lv2_phase
plt.figure()
phase,profile,profile_error = Lv2_phase.pulse_folding(times,T,T0_MJD,freq,freqdot,freqdotdot,nbins,"SWIFT")
plt.errorbar(x=phase[:-1],y=profile,yerr=profile_error,color='r',drawstyle='steps-mid')

expos = Lv2_phase.phase_exposure(times[0]-times[0],times[-1]-times[0],1/freq,nbin=nbins,gtis=np.array(gtis_conform)-times[0])
total_expos = np.array(list(expos) + list(expos))
plt.errorbar(x=phase[:-1],y=profile/total_expos,yerr=profile_error/total_expos,color='b',drawstyle='steps-mid')
plt.title(str(pathlib.Path(eventfile).name) +', exposure-corrected (using Lv2_phase)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.legend(('Folded profile','Exposure-corrected profile'),loc='best',fontsize=12)
print('Original expos:')
print(expos)
##### Using stingray.pulse.pulsar's fold_events
phase_sr,prof_sr,err_sr = fold_events(times,freq,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times[0],nbin=nbins)
phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times,freq,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times[0],expocorr=True,nbin=nbins)

total_phase_sr = list(phase_sr) + list(phase_sr+1)
total_prof_sr = list(prof_sr)*2
total_err_sr = list(err_sr)*2

total_phase_sr_expo = list(phase_sr_expo) + list(phase_sr_expo+1)
total_prof_sr_expo = list(prof_sr_expo)*2
total_err_sr_expo = list(err_sr_expo)*2

if nbins % 2 == 0:
    fft_x = np.array(list(np.arange(int(nbins/2)+1)) + list(np.arange(int(nbins/2)-1) - (int(nbins/2)-1)))
else:
    fft_x = np.array(list(np.arange(int(nbins/2)+1)) + list(np.arange(int(nbins/2)) - int(nbins/2)))

shift = np.exp(-2j*np.pi*fft_x*offset/nbins)
shifted_prof_sr = np.real(np.fft.ifft(np.fft.fft(prof_sr_expo)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile
shifted_err_sr = np.real(np.fft.ifft(np.fft.fft(err_sr_expo)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile
a = np.array(list(shifted_prof_sr)*2)/T
b = np.array(list(shifted_err_sr)*2)/T

swift_lc = open(Lv0_dirs.NGC300_2020 + 'swift_shifted_folded_curve.txt','w')
for i in range(len(total_expos)):
    swift_lc.write(str(total_phase_sr[i]) + ' ' + str(a[i]) + ' ' + str(b[i]) + '\n')
swift_lc.close()

plt.figure()
plt.errorbar(x=total_phase_sr,y=total_prof_sr/T,yerr=total_err_sr/T,color='r',drawstyle='steps-mid')
plt.errorbar(x=total_phase_sr_expo,y=total_prof_sr_expo/T,yerr=total_err_sr_expo/T,color='b',drawstyle='steps-mid')
plt.legend(('Folded profile','Exposure-corrected'),loc='best',fontsize=12)
plt.title(str(pathlib.Path(eventfile).name) +', exposure-corrected (using Stingray fold_events)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.show()
"""

"""
##### Using foldAt by PyAstronomy
plt.figure()
phase_bins = np.linspace(0,1,21)
phases = foldAt(times,1/freq,T0=times[0]-(1-phase_frac)*1/freq)

expos = Lv2_phase.phase_exposure(times[0]-times[0],times[-1]-times[0],1/freq,nbin=nbins,gtis=np.array(gtis_conform)-times[0])
total_expos = np.array(list(expos) + list(expos))
expos_index = int(phase_frac/(phase_bins[1]-phase_bins[0])) #starting point for exposures
altered_expos = np.array(list(total_expos[expos_index:]) + list(total_expos[:expos_index]))

#print('Altered expos:')
#print(altered_expos)

profile,bin_edges,binnumber = stats.binned_statistic(phases,np.ones(len(phases)),statistic='sum',bins=phase_bins)
error = np.sqrt(profile)
phase_to_2 = np.array(list(phase_bins[:-1]) + list(phase_bins+1))
profile_to_2 = np.array(list(profile)*2)
error_to_2 = np.array(list(error)*2)
plt.errorbar(phase_to_2[:-1],profile_to_2/(T*altered_expos),yerr=error_to_2/(T*altered_expos),color='b',drawstyle='steps-mid')
plt.legend(('Folded profile','Exposure-corrected'),loc='best',fontsize=12)
plt.title(str(pathlib.Path(eventfile).name) +', exposure-corrected (using PyAstronomy foldAt)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)


##### Shifting pulse profiles through a shifted FT (see Deepto's 7/20/2020 email)
if nbins % 2 == 0:
    fft_x = np.array(list(np.arange(int(nbins/2)+1)) + list(np.arange(int(nbins/2)-1) - (int(nbins/2)-1)))
else:
    fft_x = np.array(list(np.arange(int(nbins/2)+1)) + list(np.arange(int(nbins/2)) - int(nbins/2)))

shift = np.exp(-2j*np.pi*fft_x*offset/nbins)
shifted_prof_sr = np.real(np.fft.ifft(np.fft.fft(prof_sr_expo)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile
shifted_err_sr = np.real(np.fft.ifft(np.fft.fft(err_sr_expo)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile

plt.figure()
plt.errorbar(x=total_phase_sr_expo,y=total_prof_sr_expo/T,yerr=total_err_sr_expo/T,color='b',drawstyle='steps-mid')
plt.errorbar(total_phase_sr,np.array(list(shifted_prof_sr)*2)/T,yerr=np.array(list(shifted_err_sr)*2)/T,color='r',drawstyle='steps-mid')
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.title('Exposure-corrected, folded profiles for NGC 300 X-1 from Swift over May 2018 to May 2019')
plt.legend(('Folded with T0 = time of first event','Folded with T0 = inferred eclipse time/phase'),fontsize=12)
"""

"""
nbins_t = len(times)
offset = (1-0.215)*1/freq
##### Shifting pulse profiles through a shifted FT (see Deepto's 7/20/2020 email)
if nbins_t % 2 == 0:
    fft_x = np.array(list(np.arange(int(nbins_t/2)+1)) + list(np.arange(int(nbins_t/2)-1) - (int(nbins_t/2)-1)))
else:
    fft_x = np.array(list(np.arange(int(nbins_t/2)+1)) + list(np.arange(int(nbins_t/2)) - int(nbins_t/2)))

shift = np.exp(-2j*np.pi*fft_x*offset/nbins_t)
shifted_t = np.real(np.fft.ifft(np.fft.fft(times)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile
for i in range(20):
    print(times[i],shifted_t[i])

phase_sr,prof_sr,err_sr = fold_events(shifted_t,freq,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times[0],nbin=nbins)
phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(shifted_t,freq,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times[0],expocorr=True,nbin=nbins)

plt.figure()
plt.errorbar(phase_sr,prof_sr/T,color='b',drawstyle='steps-mid')
plt.errorbar(phase_sr,prof_sr_expo/T,color='r',drawstyle='steps-mid')
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
"""

#plt.show()

"""
##### Fitting 6-model step-and-ramp parameters to the folded profile
plt.figure()
plt.errorbar(x=phase[:-1],y=profile,yerr=profile_error,color='r',drawstyle='steps-mid')
plt.errorbar(x=phase[:-1],y=profile/total_expos,yerr=profile_error/total_expos,color='b',drawstyle='steps-mid')
plt.title(str(pathlib.Path(eventfile).name) +', exposure-corrected (using Lv2_phase)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.legend(('Folded profile','Exposure-corrected profile'),loc='best',fontsize=12)

start_phase = 0.45
end_phase = 1.95
phase_model = np.linspace(start_phase,end_phase,1001)
x = phase[:-1][(phase[:-1]>=start_phase)&(phase[:-1]<=end_phase)]
y = profile[(phase[:-1]>=start_phase)&(phase[:-1]<=end_phase)]/total_expos[(phase[:-1]>=start_phase)&(phase[:-1]<=end_phase)]
y_err = profile_error[(phase[:-1]>=start_phase)&(phase[:-1]<=end_phase)]/total_expos[(phase[:-1]>=start_phase)&(phase[:-1]<=end_phase)]

def piecewise_linear(x,b1,b2,b3,b4,top,bottom):
    return np.piecewise(x, [(x>=start_phase)&(x<=b1), (x>b1)&(x<=b2), (x>b2)&(x<=b3), (x>b3)&(x<=b4), (x>b4)&(x<=end_phase)], [lambda x:top, lambda x:((bottom-top)/(b2-b1)*x+bottom-(bottom-top)/(b2-b1)*b2), lambda x:bottom, lambda x:((top-bottom)/(b4-b3)*x+top-(top-bottom)/(b4-b3)*b4), lambda x:top])

pguess = np.array([1.05,1.15,1.30,1.45,0.0011,0.0003])
popt,pcov = curve_fit(piecewise_linear,x,y,p0=pguess)#,sigma=y_err)
print(popt)
print(np.diag(np.sqrt(pcov))/popt*100)
plt.plot(phase_model,piecewise_linear(phase_model,*popt),'k-')
"""

#plt.show()

########################### DOING CHI^2 EXPLORATION ############################

def lorentzian(f, f0, a, gam,const):
    x = (f-f0)/(gam/2)
    return a * 1/(1+x**2) + const

def gaussian(f,f0,a,sig,const):
    return a * np.exp( -(f-f0)**2/(2*sig**2) ) + const

def sum(f,f0,a,gam,b,sig,const):
    x = (f-f0)/(gam/2)
    return a * 1/(1+x**2) + b * np.exp( -(f-f0)**2/(2*sig**2)) + const

"""
nbins=20
chi2 = []
freqs = np.arange(8.25e-6,8.7e-6,0.01e-6)
#freqs = np.arange(-9e-17,-1e-18,1e-20)
for i in tqdm(range(len(freqs))):
    phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times,freqs[i],gtis=np.array(gtis_conform),ref_time=times[0],expocorr=True,nbins=nbins)
    chi2_freq = Lv2_phase.get_chi2(prof_sr_expo,err_sr_expo)
    chi2.append( chi2_freq )
"""

"""
freqs_filter = freqs[(freqs>=8.4693e-6)&(freqs<=8.472e-6)] #8.47 to 8.47275 for 1-year data
chi2_filter = np.array(chi2)[(freqs>=8.4693e-6)&(freqs<=8.472e-6)]

freq_model = np.linspace(8.4693e-6,8.472e-6,1001)

pguess_l = np.array([8.4706e-6,650,0.002e-6])
popt_l,pcov_l = curve_fit(lorentzian,freqs_filter,chi2_filter,p0=pguess_l)
print(popt_l)
print(np.sqrt(np.diag(pcov_l)))

pguess_g = np.array([8.4706e-6,650,0.002e-6])
popt_g,pcov_g = curve_fit(gaussian,freqs_filter,chi2_filter,p0=pguess_g)
print(popt_g)
print(np.sqrt(np.diag(pcov_g)))

pguess_s = np.array([8.4706e-6,650,0.002e-6,600,0.002e-6])
popt_s,pcov_s = curve_fit(sum,freqs_filter,chi2_filter,p0=pguess_s)
print(popt_s)
print(np.sqrt(np.diag(pcov_s)))
"""

"""
fig,ax = plt.subplots()

def pdot_to_fdot(pdot):
    return -pdot/(1/8.4712e-6)**2

def fdot_to_pdot(fdot):
    return (-fdot/(8.4712e-6)**2)/1e-7
chi2 = np.array(chi2)
#secax = ax.secondary_xaxis('top',functions=(fdot_to_pdot,pdot_to_fdot))
#secax.set_xlabel('Period Derivative (1E-7 s/s)',fontsize=12)
print(np.max(chi2),freqs[chi2==np.max(chi2)])
ax.plot(freqs,chi2,'rx-')
#ax.axvline(x=-5.60e-17,lw=0.5,alpha=0.5,color='k')
#ax.axvline(x=-2.80e-17,lw=0.5,alpha=0.5,color='k')
ax.axhline(y=869.357,lw=0.5,alpha=0.5,color='b')
#plt.plot(freq_model,lorentzian(freq_model,popt_l[0],popt_l[1],popt_l[2]),'b-')
#plt.plot(freq_model,gaussian(freq_model,popt_g[0],popt_g[1],popt_g[2]),'k-')
#plt.plot(freq_model,sum(freq_model,popt_s[0],popt_s[1],popt_s[2],popt_s[3],popt_s[4]),'m-')
ax.set_xlabel('Frequency Derivative (Hz/s)',fontsize=12)
ax.set_ylabel('chi^2 [ sum( (profile-mean)^2/error^2) ]',fontsize=12)
#plt.legend(('manual chi^2','Lorentzian fit','Gaussian fit','L+G'),fontsize=12)
plt.show()
"""

def sinecurve(x,a,T,phi,c):
    return a*np.sin(2*np.pi/T*x+phi) + c

##### Exploring reduced data from XMM-Newton
##### Doing sine curve fitting with the RATE data
"""
xmm_lc1 = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_lccorr.lc'
rebinned_t_xmm1 = fits.open(xmm_lc1)[1].data['TIME']
rebinned_rate_xmm1 = fits.open(xmm_lc1)[1].data['RATE']
rebinned_err_xmm1 = fits.open(xmm_lc1)[1].data['ERROR']
xmm_lc2 = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_lccorr.lc'
rebinned_t_xmm2 = fits.open(xmm_lc2)[1].data['TIME']
rebinned_rate_xmm2 = fits.open(xmm_lc2)[1].data['RATE']
rebinned_err_xmm2 = fits.open(xmm_lc2)[1].data['ERROR']

mjd_x1_xmm = fits.open(xmm_lc1)[1].header['MJDREF'] + np.array(list(rebinned_t_xmm1) + list(rebinned_t_xmm2))/86400
rebinned_t_xmm = np.array(list(rebinned_t_xmm1) + list(rebinned_t_xmm2))
rebinned_rate_xmm = np.array(list(rebinned_rate_xmm1) + list(rebinned_rate_xmm2))
rebinned_err_xmm = np.array(list(rebinned_err_xmm1) + list(rebinned_err_xmm2))

pguess = np.array([0.2,120e3,-0.5,0.2])
popt,pcov = curve_fit(sinecurve,rebinned_t_xmm,rebinned_rate_xmm,sigma=rebinned_err_xmm,absolute_sigma=True,p0=pguess)
print('amplitude: ' + str(popt[0]))
print('period: ' + str(popt[1]))
print('freq: ' + str(1/popt[1]))
print('phase shift: ' + str(popt[2]))
print('offset: ' + str(popt[3]))
print(np.sqrt(np.diag(pcov)))

plt.plot(rebinned_t_xmm,rebinned_rate_xmm,'r-')
plt.plot(rebinned_t_xmm,sinecurve(rebinned_t_xmm,*popt),'b-')
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Rate (counts/s)',fontsize=12)

print('subset1')
subset_t = rebinned_t_xmm[(rebinned_t_xmm>=5.9845e8)&(rebinned_t_xmm<=5.98475e8)]
subset_rate = sinecurve(rebinned_t_xmm,*popt)[(rebinned_t_xmm>=5.9845e8)&(rebinned_t_xmm<=5.98475e8)]
print(np.min(subset_rate))
print(subset_t[subset_rate==np.min(subset_rate)][0])
print(50814 + subset_t[subset_rate==np.min(subset_rate)][0]/86400)

print('subset2')
subset_t = rebinned_t_xmm[rebinned_t_xmm>=5.9855e8]
subset_rate = sinecurve(rebinned_t_xmm,*popt)[rebinned_t_xmm>=5.9855e8]
print(np.min(subset_rate))
print(subset_t[subset_rate==np.min(subset_rate)][0])
print(50814 + subset_t[subset_rate==np.min(subset_rate)][0]/86400)

plt.show()
"""

"""
tbins = np.arange(times[0],times[-1]+100,100)
summed_data, bin_edges, binnumber = stats.binned_statistic(times,np.ones(len(times)),statistic='sum',bins=tbins)

t_used = tbins[:-1][summed_data>0]
counts_used = summed_data[summed_data>0]

pguess = np.array([10,120e3,5,15])
#popt,pcov = curve_fit(sinecurve,tbins[:-1],summed_data,sigma=np.sqrt(summed_data),absolute_sigma=True,p0=pguess)
popt,pcov = curve_fit(sinecurve,t_used,counts_used,sigma=np.sqrt(counts_used),absolute_sigma=True,p0=pguess,maxfev=10000)
print('amplitude: ' + str(popt[0]))
print('period: ' + str(popt[1]))
print('freq: ' + str(1/popt[1]))
print('phase shift: ' + str(popt[2]))
print('offset: ' + str(popt[3]))
print(np.sqrt(np.diag(pcov)))

plt.plot(t_used,counts_used,'r-')
plt.plot(t_used,sinecurve(t_used,*popt),'b-')
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Counts',fontsize=12)


print('subset1')
subset_t = t_used[(t_used>=5.9845e8)&(t_used<=5.98475e8)]
subset_rate = sinecurve(t_used,*popt)[(t_used>=5.9845e8)&(t_used<=5.98475e8)]
print(np.min(subset_rate))
print(subset_t[subset_rate==np.min(subset_rate)][0])
print(50814 + subset_t[subset_rate==np.min(subset_rate)][0]/86400)

print('subset2')
subset_t = tbins[:-1][tbins[:-1]>=5.9855e8]
subset_rate = sinecurve(tbins[:-1],*popt)[tbins[:-1]>=5.9855e8]
print(np.min(subset_rate))
print(subset_t[subset_rate==np.min(subset_rate)][0])
print(50814 + subset_t[subset_rate==np.min(subset_rate)][0]/86400)

plt.show()
"""

###############################################################################
######################### Folding the XMM-Newton data #########################
###############################################################################

pb = 1/8.4712e-6
freqdot = 0
freqdotdot = 0
nbins = 20

gtis_conform = []
for i in range(len(gtis_data_xmm)):
    gtis_conform.append([gtis_data_xmm[i][0],gtis_data_xmm[i][1]])

"""
nbins=20
chi2 = []
freqs = np.arange(8e-6,9e-6,0.001e-6)
#freqs = np.arange(-9e-17,-1e-18,1e-20)
for i in tqdm(range(len(freqs))):
    phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times_xmm,freqs[i],gtis=np.array(gtis_conform),ref_time=times_xmm[0],expocorr=True,nbins=nbins)
    chi2_freq = Lv2_phase.get_chi2(prof_sr_expo,err_sr_expo)
    chi2.append( chi2_freq )

fig,ax = plt.subplots()

chi2 = np.array(chi2)
#secax = ax.secondary_xaxis('top',functions=(fdot_to_pdot,pdot_to_fdot))
#secax.set_xlabel('Period Derivative (1E-7 s/s)',fontsize=12)
#print(np.max(chi2),freqs[chi2==np.max(chi2)])
ax.plot(freqs,chi2,'rx-')
ax.set_xlabel('Frequency (Hz)',fontsize=12)
ax.set_ylabel('chi^2 [ sum( (profile-mean)^2/error^2) ]',fontsize=12)

plt.show()
"""

##### Using Lv2_phase
plt.figure()
phase,profile,profile_error = Lv2_phase.pulse_folding(times_xmm,T_xmm,T0_MJD_xmm,1/pb,freqdot,freqdotdot,nbins,"XMM")
plt.errorbar(x=phase[:-1],y=profile,yerr=profile_error,color='r',drawstyle='steps-mid')

expos = Lv2_phase.phase_exposure(times_xmm[0]-times_xmm[0],times_xmm[-1]-times_xmm[0],pb,nbin=nbins,gtis=np.array(gtis_conform)-times_xmm[0])
total_expos = np.array(list(expos) + list(expos))
plt.errorbar(x=phase[:-1],y=profile/total_expos,yerr=profile_error/total_expos,color='b',drawstyle='steps-mid')
plt.title(str(pathlib.Path(eventfile_xmm).name) +', exposure-corrected (using Lv2_phase)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.legend(('Folded profile','Exposure-corrected profile'),loc='best',fontsize=12)

##### Using stingray.pulse.pulsar's fold_events
phase_sr,prof_sr,err_sr = fold_events(times_xmm,1/pb,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times_xmm[0]-phaseoff*pb,nbin=nbins)
phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times_xmm,1/pb,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times_xmm[0]-phaseoff*pb,expocorr=True,nbin=nbins)

total_phase_sr = np.array(list(phase_sr) + list(phase_sr+1))
total_prof_sr = np.array(list(prof_sr)*2)
total_err_sr = np.array(list(err_sr)*2)

total_phase_sr_expo = np.array(list(phase_sr_expo) + list(phase_sr_expo+1))
total_prof_sr_expo = np.array(list(prof_sr_expo)*2)
total_err_sr_expo = np.array(list(err_sr_expo)*2)

plt.figure()
plt.errorbar(x=total_phase_sr,y=total_prof_sr/T_xmm,yerr=total_err_sr/T_xmm,color='r',drawstyle='steps-mid')
plt.errorbar(x=total_phase_sr_expo,y=total_prof_sr_expo/T_xmm,yerr=total_err_sr_expo/T_xmm,color='b',drawstyle='steps-mid')
plt.legend(('Folded profile','Exposure-corrected'),loc='best',fontsize=12)
plt.title(str(pathlib.Path(eventfile_xmm).name) +', exposure-corrected (using Stingray fold_events)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)

for i in range(len(total_phase_sr_expo)):
    print(total_phase_sr_expo[i],total_prof_sr_expo[i]/T_xmm,total_err_sr_expo[i]/T_xmm)

def step_n_ramp(phase,prof,prof_err,start_phase,end_phase,pguess):
    """
    Fitting 6-model step-and-ramp parameters to the folded profile
    """
    phase_model = np.linspace(start_phase,end_phase,101)
    x = phase[(phase>=start_phase)&(phase<=end_phase)]
    y = prof[(phase>=start_phase)&(phase<=end_phase)]
    y_err = prof_err[(phase>=start_phase)&(phase<=end_phase)]

    def piecewise_linear(x,b1,b2,b3,b4,top,bottom):
        return np.piecewise(x, [(x>=start_phase)&(x<=b1), (x>b1)&(x<=b2), (x>b2)&(x<=b3), (x>b3)&(x<=b4), (x>b4)&(x<=end_phase)], [lambda x:top, lambda x:((bottom-top)/(b2-b1)*x+bottom-(bottom-top)/(b2-b1)*b2), lambda x:bottom, lambda x:((top-bottom)/(b4-b3)*x+top-(top-bottom)/(b4-b3)*b4), lambda x:top])

    plt.figure()
    popt,pcov = curve_fit(piecewise_linear,x,y,p0=pguess,sigma=y_err,absolute_sigma=True)
    pars = popt
    pars_err = np.diag(np.sqrt(pcov))
    print('Top: ' + str(pars[4]) + ' +- ' + str(pars_err[4]))
    print('Bottom: ' + str(pars[5]) + ' +- ' + str(pars_err[5]))
    print('Vertex 1: ' + str(pars[0]) + ' +- ' + str(pars_err[0]))
    print('Vertex 2: ' + str(pars[1]) + ' +- ' + str(pars_err[1]))
    print('Vertex 3: ' + str(pars[2]) + ' +- ' + str(pars_err[2]))
    print('Vertex 4: ' + str(pars[3]) + ' +- ' + str(pars_err[3]))
    plt.plot(phase_model,piecewise_linear(phase_model,*popt),'k-')

    ##### Plotting the folded profiles themselves
    plt.errorbar(x=phase,y=prof,yerr=prof_err,color='r',drawstyle='steps-mid')
    plt.title('Exposure-corrected (profiles from Stingray)',fontsize=12)
    plt.xlabel('Phase',fontsize=12)
    plt.ylabel('Counts/s',fontsize=12)
    plt.legend(('Piecewise fit','Exposure-corrected profile'),loc='best',fontsize=12)

#step_n_ramp(total_phase_sr_expo,total_prof_sr_expo/T_xmm,total_err_sr_expo/T_xmm,0.225,1.775,np.array([0.65,0.75,1,1.25,0.016,0.0035]))

plt.show()



###############################################################################
######################## Combining Swift and XMM-Newton #######################
###############################################################################
"""
#pb = 117403.24413
#pb = 1/8.47145464e-6
pb = 1/8.4712e-6
freqdot = 0
freqdotdot = 0
nbins = 20

MJDREFI = fits.open(eventfile)[1].header['MJDREFI'] #Swift
MJDREFF = fits.open(eventfile)[1].header['MJDREFF'] #Swift
MJDREF = fits.open(eventfile_xmm)[1].header['MJDREF'] #XMM-Newton
diff_swiftxmm = (MJDREFI+MJDREFF-MJDREF)*86400

gtis_conform = []
for i in range(len(gtis_data_xmm)): #for each GTI in the XMM data
    gtis_conform.append([gtis_data_xmm[i][0],gtis_data_xmm[i][1]])
for i in range(len(gtis_data)): #for each GTI in the Swift data
    gtis_conform.append([gtis_data[i][0]+diff_swiftxmm,gtis_data[i][1]+diff_swiftxmm])

times_all = np.array(list(times_xmm) + list(diff_swiftxmm + times))
T_all = T + T_xmm
T0_MJD_all = T0_MJD_xmm
"""

"""
##### chi^2 exploration
chi2 = []
freqs = np.arange(8.4e-6,8.500000e-6,0.01e-6)
for i in tqdm(range(len(freqs))):
    phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times_all,freqs[i],gtis=np.array(gtis_conform),ref_time=times_all[0],expocorr=True,nbins=nbins)
    chi2_freq = Lv2_phase.get_chi2(prof_sr_expo,err_sr_expo)
    chi2.append( chi2_freq )

plt.figure()
plt.plot(freqs/1e-6,chi2,'rx-')
plt.axvline(x=8.4712,lw=0.5,alpha=0.5)
plt.xlabel('Frequency (microHz)',fontsize=12)
plt.ylabel('chi^2 [ sum( (profile-mean)^2/error^2) ]',fontsize=12)
plt.legend(('chi^2 exploration','8.4712E-6 Hz, freq. from Swift'),loc='best')
plt.show()
"""

"""
##### Using Lv2_phase
plt.figure()
phase,profile,profile_error = Lv2_phase.pulse_folding(times_all,T_all,T0_MJD_all,1/pb,freqdot,freqdotdot,nbins,"XMM")
plt.errorbar(x=phase[:-1],y=profile,yerr=profile_error,color='r',drawstyle='steps-mid')

expos = Lv2_phase.phase_exposure(times_all[0]-times_all[0],times_all[-1]-times_all[0],pb,nbin=nbins,gtis=np.array(gtis_conform)-times_all[0])
total_expos = np.array(list(expos) + list(expos))
plt.errorbar(x=phase[:-1],y=profile/total_expos,yerr=profile_error/total_expos,color='b',drawstyle='steps-mid')
plt.title('XMM + Swift, exposure-corrected (using Lv2_phase)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.legend(('Folded profile','Exposure-corrected profile'),loc='best',fontsize=12)

##### Using stingray.pulse.pulsar's fold_events
phase_sr,prof_sr,err_sr = fold_events(times_all,1/pb,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times_all[0],nbin=nbins)
phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times_all,1/pb,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times_all[0],expocorr=True,nbin=nbins)

total_phase_sr = list(phase_sr) + list(phase_sr+1)
total_prof_sr = list(prof_sr)*2
total_err_sr = list(err_sr)*2

total_phase_sr_expo = list(phase_sr_expo) + list(phase_sr_expo+1)
total_prof_sr_expo = list(prof_sr_expo)*2
total_err_sr_expo = list(err_sr_expo)*2

plt.figure()
plt.errorbar(x=total_phase_sr,y=total_prof_sr/T_all,yerr=total_err_sr/T_all,color='r',drawstyle='steps-mid')
plt.errorbar(x=total_phase_sr_expo,y=total_prof_sr_expo/T_all,yerr=total_err_sr_expo/T_all,color='b',drawstyle='steps-mid')
plt.legend(('Folded profile','Exposure-corrected'),loc='best',fontsize=12)
plt.title('XMM + Swift, exposure-corrected (using Stingray fold_events)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)

plt.show()
"""


fig,ax = plt.subplots()

def pdot_to_fdot(pdot):
    return -pdot/(1/8.4712e-6)**2

def fdot_to_pdot(fdot):
    return (-fdot/(8.4712e-6)**2)/1e-7

##### Independent sums of chi^2 (Deepto's suggestion)
nbins = 20
"""
chi2 = []
chi2_swift_all = []
chi2_xmm_all = []
#freqs = np.arange(1.0/(40.0*3600.0),1.0/(20.0*3600.0),1e-11)
freqs = np.arange(8.45e-6,8.50e-6,1e-10)
freqdots = np.arange(1e-18,9e-17,1e-20)
chi2_all_write = open('/Volumes/Samsung_T5/NGC300_XMMdata/placeholder.txt','w')
chi2_swift_write = open('/Volumes/Samsung_T5/NGC300_XMMdata/placeholder2','w')
chi2_xmm_write = open('/Volumes/Samsung_T5/NGC300_XMMdata/placeholder3','w')

for i in tqdm(range(len(freqdots))):
    for j in tqdm(range(len(freqdots))):
        ## Swift
        phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times,freqs[i],freqdots[j],gtis=np.array(gtis_conform),ref_time=times[0],expocorr=True,nbins=nbins)
        chi2_swift = Lv2_phase.get_chi2(prof_sr_expo,err_sr_expo)
        chi2_swift_all.append(chi2_swift)

        chi2_swift_write.write(str(freqs[i]) + ' ' + str(freqdots[j]) + ' ' + str(chi2_swift) + '\n')

        ## XMM-Newton
        phase_sr_expo_xmm,prof_sr_expo_xmm,err_sr_expo_xmm = fold_events(times_xmm,freqs[i],freqdots[j],gtis=np.array(gtis_conform_xmm),ref_time=times_xmm[0],expocorr=True,nbins=nbins)
        chi2_xmm = Lv2_phase.get_chi2(prof_sr_expo_xmm,err_sr_expo_xmm)
        chi2_xmm_all.append(chi2_xmm)
        chi2_xmm_write.write(str(freqs[i]) + ' ' + str(freqdots[j]) + ' ' + str(chi2_xmm) + '\n')

        chi2.append( chi2_swift + chi2_xmm )
        chi2_all_write.write(str(freqs[i]) + ' ' + str(freqdots[j]) + ' ' + str(chi2_swift + chi2_xmm) + '\n')

chi2_all_write.close()
chi2_swift_write.close()
chi2_xmm_write.close()

secax = ax.secondary_xaxis('top',functions=(fdot_to_pdot,pdot_to_fdot))
secax.set_xlabel('Period Derivative (1E-7 s/s)',fontsize=12)

ax.plot(freqdots,chi2,'kx-')
ax.plot(freqdots,chi2_swift_all,'rx-',lw=0.5,alpha=0.5)
ax.plot(freqdots,chi2_xmm_all,'bx-',lw=0.5,alpha=0.5)
#plt.yscale('log')
ax.legend(('Swift+XMM','Swift','XMM'),fontsize=12)
ax.set_xlabel('Frequency Derivative (Hz/s)',fontsize=12)
ax.set_ylabel('chi^2 [ sum( (profile-mean)^2/error^2) ]',fontsize=12)
mplcursors.cursor(hover=True)

ax.axvline(x=5.60e-17,lw=0.5,alpha=0.5,color='k')
ax.axvline(x=2.80e-17,lw=0.5,alpha=0.5,color='k')
#ax.axhline(y=869.357,lw=0.5,alpha=0.5,color='b') for 1e-11 Hz spacing
#ax.axhline(y=11830.79495183693,lw=0.5,alpha=0.5,color='b') for 1e-11 spacing
ax.axhline(y=734.51,lw=0.5,alpha=0.5,color='b')
ax.axhline(y=11689.2,lw=0.5,alpha=0.5,color='b')

plt.show()
"""

"""
##### Plotting results from the chi^2 exploration
## secondary axis reference: https://matplotlib.org/3.1.0/gallery/subplots_axes_and_figures/secondary_axis.html

freqs,chi2_all = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_all_dec16-may19.txt',usecols=(0,1),unpack=True)
freqs,chi2_swift = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_swift_dec16-may19.txt',usecols=(0,1),unpack=True)
freqs,chi2_xmm = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_xmm_dec16-may19.txt',usecols=(0,1),unpack=True)

fig,ax = plt.subplots()

ax.plot(freqs/1e-6,chi2_all,'kx-')
ax.plot(freqs/1e-6,chi2_swift,'rx-',lw=0.5,alpha=0.5)
ax.plot(freqs/1e-6,chi2_xmm,'bx-',lw=0.5,alpha=0.5)

def time_to_freq(x):
    return (1/x)

def freq_to_time(x):
    return (1/x)/3600*1e6

secax = ax.secondary_xaxis('top',functions=(freq_to_time,time_to_freq))
secax.set_xlabel('Time (h)',fontsize=12)

freqs_fit = freqs[(freqs>=8.46e-6)&(freqs<=8.48e-6)]
chi_fit = np.array(chi2_all)[(freqs>=8.46e-6)&(freqs<=8.48e-6)]

#pguess = np.array([8.4712,400,0.02,200]) #for Swift
pguess = np.array([8.4712,500,0.02,11100]) #for all

#popt,pcov = curve_fit(lorentzian,freqs_fit/1e-6,chi_fit,p0=pguess)
#print('Lorentzian')
#print(popt)
#print(np.sqrt(np.diag(pcov)))
#plt.plot(freqs_fit/1e-6,lorentzian(freqs_fit/1e-6,*popt),'r-')

#popt,pcov = curve_fit(gaussian,freqs_fit/1e-6,chi_fit,p0=pguess)
#print('Gaussian')
#print(popt)
#print(np.sqrt(np.diag(pcov)))
#plt.plot(freqs_fit/1e-6,gaussian(freqs_fit/1e-6,*popt),'b-')

#plt.legend(('All data','Swift data','XMM data','Lorentzian fit','Gaussian fit'),loc='best')
ax.set_xlabel('Frequency (microHz)',fontsize=12)
ax.set_ylabel('chi^2 [ sum( (profile-mean)^2/error^2) ]',fontsize=12)

plt.show()
"""


##### Doing contour plots for the 2D P-Pdot exploration
### do a PDOT version too??

def summarize_2d(space,chi2_type,posneg):
    """
    Summarizing the information from the 2D chi^2 exploration involving frequency and
    the frequency derivative

    space - whether in frequency space or period space
    chi2_type - 'XMM', 'Swift', or 'all'
    posneg - positive fdot or negative fdot
    """
    plt.figure()
    if chi2_type == 'XMM':
        plt.title('XMM')
        if posneg == 'pos':
            freq,freqdot,chi2 = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_xmm_fine_ffdot.txt',usecols=(0,1,2),unpack=True)
            if space == 'frequency':
                plt.axhline(y=5.6e-17,color='w',lw=1,alpha=0.5)
                plt.axhline(y=2.8e-17,color='w',lw=1,alpha=0.5)
            elif space == 'period':
                plt.axhline(y=-3.9,color='w',lw=1,alpha=0.5)
                plt.axhline(y=-7.8,color='w',lw=1,alpha=0.5)
        elif posneg == 'neg':
            freq,freqdot,chi2 = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_xmm_fine_ffdotneg.txt',usecols=(0,1,2),unpack=True)
            if space == 'frequency':
                plt.axhline(y=-5.6e-17,color='w',lw=1,alpha=0.5)
                plt.axhline(y=-2.8e-17,color='w',lw=1,alpha=0.5)
            elif space == 'period':
                plt.axhline(y=3.9,color='w',lw=1,alpha=0.5)
                plt.axhline(y=7.8,color='w',lw=1,alpha=0.5)
    elif chi2_type == 'Swift':
        plt.title('Swift')
        if posneg == 'pos':
            freq,freqdot,chi2 = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_swift_fine_ffdot.txt',usecols=(0,1,2),unpack=True)
            if space == 'frequency':
                plt.axhline(y=5.6e-17,color='w',lw=1,alpha=0.5)
                plt.axhline(y=2.8e-17,color='w',lw=1,alpha=0.5)
            elif space == 'period':
                plt.axhline(y=-3.9,color='w',lw=1,alpha=0.5)
                plt.axhline(y=-7.8,color='w',lw=1,alpha=0.5)
        elif posneg == 'neg':
            freq,freqdot,chi2 = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_swift_fine_ffdotneg.txt',usecols=(0,1,2),unpack=True)
            if space == 'frequency':
                plt.axhline(y=-5.6e-17,color='w',lw=1,alpha=0.5)
                plt.axhline(y=-2.8e-17,color='w',lw=1,alpha=0.5)
            elif space == 'period':
                plt.axhline(y=3.9,color='w',lw=1,alpha=0.5)
                plt.axhline(y=7.8,color='w',lw=1,alpha=0.5)
    elif chi2_type == 'all':
        plt.title('all')
        if posneg == 'pos':
            freq,freqdot,chi2 = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_all_fine_ffdot.txt',usecols=(0,1,2),unpack=True)
            if space == 'frequency':
                plt.axhline(y=5.6e-17,color='w',lw=1,alpha=0.5)
                plt.axhline(y=2.8e-17,color='w',lw=1,alpha=0.5)
            elif space == 'period':
                plt.axhline(y=-3.9,color='w',lw=1,alpha=0.5)
                plt.axhline(y=-7.8,color='w',lw=1,alpha=0.5)
        elif posneg == 'neg':
            freq,freqdot,chi2 = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_all_fine_ffdotneg.txt',usecols=(0,1,2),unpack=True)
            if space == 'frequency':
                plt.axhline(y=-5.6e-17,color='w',lw=1,alpha=0.5)
                plt.axhline(y=-2.8e-17,color='w',lw=1,alpha=0.5)
            elif space == 'period':
                plt.axhline(y=3.9,color='w',lw=1,alpha=0.5)
                plt.axhline(y=7.8,color='w',lw=1,alpha=0.5)
    else:
        raise ValueError("Make sure that chi2_type is one of 'XMM', 'Swift', or 'all'! Also that posneg is either 'pos' or 'neg'.")

    freq_1d = np.arange(8.45e-6,8.5e-6,1e-10)
    if posneg == 'pos':
        freqdot_1d = np.arange(1e-18,9e-17,1e-20)
    elif posneg == 'neg':
        freqdot_1d = np.arange(-9e-17,-1e-18,1e-20)

    N_freq = len(freq_1d)
    N_freqdot = len(freqdot_1d)

    freq_grid,freqdot_grid = np.meshgrid(freq_1d,freqdot_1d)
    chi2_reshape = np.reshape(chi2,(N_freqdot,N_freq),order='F')

    if space == 'frequency':
        print('Maximum chi^2 is ' + str(round(np.max(chi2),2)) + ', at freq ' + str(round(freq[chi2==np.max(chi2)][0]/1e-6,4)) + ' microHz, at freqdot ' + str(round(freqdot[chi2==np.max(chi2)][0]/1e-17,4)) + 'E-17 Hz/s')

        plt.pcolormesh(freq_grid/1e-6,freqdot_grid,chi2_reshape,cmap='gist_heat',vmin=np.min(chi2),vmax=np.max(chi2))
        cbar = plt.colorbar()
        cbar.set_label('chi^2 [ sum( (profile-mean)^2/error^2) ]')
        mplcursors.cursor(hover=True)
        plt.xlabel('Frequency (microHz)',fontsize=12)
        plt.ylabel('Frequency Derivative (Hz/s)',fontsize=12)

        plt.show()

    elif space == 'period':
        period = 1/freq
        pdot = -freqdot/freq**2

        period_1d = np.linspace(np.min(period),np.max(period),len(freq_1d)+1)
        pdot_1d = np.linspace(np.min(pdot),np.max(pdot),len(freqdot_1d)+1)

        p_grid,pdot_grid = np.meshgrid(period_1d,pdot_1d)

        print('Maximum chi^2 is ' + str(round(np.max(chi2),2)) + ', at period ' + str(round(period[chi2==np.max(chi2)][0],4)) + ' seconds, at pdot ' + str(round(pdot[chi2==np.max(chi2)][0]/1e-7,4)) + 'E-7 s/s')

        plt.pcolormesh(p_grid,pdot_grid/1e-7,chi2_reshape,cmap='gist_heat',vmin=np.min(chi2),vmax=np.max(chi2))
        cbar = plt.colorbar()
        cbar.set_label('chi^2 [ sum( (profile-mean)^2/error^2) ]')
        mplcursors.cursor(hover=True)
        plt.xlabel('Period (ks)',fontsize=12)
        plt.ylabel('Period Derivative (1E-7 s/s)',fontsize=12)

        plt.show()

    else:
        raise ValueError("Make sure that the variable space is either 'frequency' or 'period'!")

#summarize_2d('period','Swift','pos')
#summarize_2d('period','XMM','pos')
#summarize_2d('period','all','pos')
#summarize_2d('period','Swift','neg')
#summarize_2d('period','XMM','neg')
#summarize_2d('period','all','neg')

def add_orbphase(eventfile,pb):
    """
    Add the orbital phase column into the event file

    eventfile - path to the event file
    pb - orbital period in seconds
    """
