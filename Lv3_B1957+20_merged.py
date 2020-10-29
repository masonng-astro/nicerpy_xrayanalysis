#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 7:17pm 2019

Generating averaged power spectra for B1957+20/J1959+2048! Mainly trying to
find the pulsation.

"""
from __future__ import division, print_function
import numpy as np
import time
from tqdm import tqdm
from astropy.io import fits
import glob
import pathlib
import subprocess

import matplotlib.pyplot as plt
from stingray.pulse.pulsar import pulse_phase,phase_exposure,fold_events

eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000021/merged000021_nicersoft_bary_phase.evt'
times = fits.open(eventfile)[1].data['TIME']
pulse_phase = fits.open(eventfile)[1].data['PULSE_PHASE']
gtis = fits.open(eventfile)[2].data
T = sum([ gtis[i]['STOP']-gtis[i]['START'] for i in range(len(gtis)) ]) #exposure time

gtis_conform = []
for i in range(len(gtis)):
    gtis_conform.append([gtis[i][0],gtis[i][1]]) #conform to the input that Stingray uses

freq = 622.12202499673376738
freqdot = -6.5114520166830696246e-15
freqdotdot = 0
nbins = 20

phase_sr,prof_sr,err_sr = fold_events(times,freq,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times[0],nbin=nbins)
phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times,freq,freqdot,freqdotdot,gtis=np.array(gtis_conform),ref_time=times[0],expocorr=True,nbin=nbins)

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
plt.show()


"""
DEPRECATED FOR NOW - 8/17

par_file = Lv0_dirs.NICERSOFT_DATADIR + 'B1957+20.par' #parameter file for demodulation
T_asc_paper = 51260.200512649951936
PB = 0.3819666099158615735
grid = np.arange(0,1,0.01)
#grid = np.arange(0.184,0.198,0.002)
#grid = np.arange(0.202,0.22,0.002)

print('Exploring the orbit now...')
for i in tqdm(range(len(grid))):

    print('Trying out phase = ' + str(grid[i]))
    T_asc = T_asc_paper + grid[i]*PB
    par_contents = open(par_file,'r')
    contents = par_contents.read()
    contents = contents.split('\n')
    par_contents.close()

    newstring = 'TASC ' + str(T_asc) + '    #  TASC     Epoch of ascending node passage (MJD)'
    par_contents = open(par_file,'w')
    for j in range(len(contents)):
        if j != 13:
            par_contents.write(contents[j]+'\n')
        else:
            par_contents.write(newstring + '\n')
    par_contents.close()

    pyfile = 'Lv3_average_ps_main.py'
    pyfile_contents = open(pyfile,'r')
    contents = pyfile_contents.read().split('\n')
    pyfile_contents.close()

    newstring = "    pngname = '" + Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/merged000005/' + str(grid[i]) + "PB.png'"
    pyfile_contents = open('Lv3_average_ps_main.py','w')
    for j in range(len(contents)):
        if j != 123:
            pyfile_contents.write(contents[j] + '\n')
        else:
            pyfile_contents.write(newstring + '\n')
    pyfile_contents.close()

    execfile("Lv3_average_ps_main.py")
"""

"""
obsids = ['10301801'+str(i) for i in range(49,58)]

demod = True
merged_id = '000005' #need to be very careful that I know what the next one is!
segment_length = 5000 #segment length
PI1 = 50 #lower bound for PI
PI2 = 450 #upper bound for PI
tbin = 0.00025 #bin size in s
N = Lv3_detection_level.N_trials(tbin,2) #2 because looking in between 621 and 623 Hz
threshold = 1 #threshold for counts in each segment
W = 1 #number of consecutive Fourier bins to average over
starting_freq = 10 #for noise_hist
segment_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/merged' + merged_id + '/accelsearch_' + str(segment_length).zfill(5) + 's/'

#### ONLY NEED TO RUN THESE 4 LINES ONCE!!!
#Lv2_average_merge_ps_methods.merging(obsids)
#Lv2_average_merge_ps_methods.merging_GTIs(obsids,merged_id)
#Lv2_average_merge_ps_methods.get_gti_file(merged_id,segment_length)
#Lv2_average_merge_ps_methods.niextract_gti_time_energy(merged_id,segment_length,PI1,PI2)

par_file = Lv0_dirs.NICERSOFT_DATADIR + 'B1957+20.par' #parameter file for demodulation
T_asc_paper = 48196.0635242
PB = 0.3819666069
grid = np.arange(0,1,0.02)

print('Exploring the orbit now...')
for i in tqdm(range(len(grid))):

    print('Trying out phase = ' + str(grid[i]))
    T_asc = T_asc_paper + grid[i]*PB
    par_contents = open(par_file,'r')
    contents = par_contents.read()
    contents = contents.split('\n')
    par_contents.close()

    newstring = 'TASC ' + str(T_asc) + '    #  TASC     Epoch of ascending node passage (MJD)'
    par_contents = open(par_file,'w')
    for j in range(len(contents)):
        if j != 13:
            par_contents.write(contents[j]+'\n')
        else:
            par_contents.write(newstring + '\n')
    par_contents.close()

    Lv2_average_merge_ps_methods.do_demodulate(merged_id,segment_length,par_file,PI2)

    ##### do NICERfits2presto
    demod_files = glob.glob(segment_dir + '*demod.evt')
    for i in tqdm(range(len(demod_files))):
        try:
            subprocess.check_call(['nicerfits2presto.py','--dt='+str(tbin),demod_files[i]])
        except (ValueError,subprocess.CalledProcessError):
            pass

    Lv2_average_merge_ps_methods.edit_inf(merged_id,tbin,segment_length)
    Lv2_average_merge_ps_methods.edit_binary(merged_id,tbin,segment_length)
    Lv2_average_merge_ps_methods.realfft(merged_id,segment_length)

    f,ps,ps_bins,N_greaterthanP,M = Lv2_average_merge_ps_methods.average_ps(merged_id,segment_length,demod,PI1,PI2,tbin,threshold,starting_freq,W)

    power_required_3 = Lv3_detection_level.power_for_sigma(3,N,M,W) #power required for significance
    power_required_4 = Lv3_detection_level.power_for_sigma(4,N,M,W) #power required for significance

    #plt.figure(1)
    plt.plot(f,ps,'rx-')
    plt.axhline(y=power_required_3,lw=0.8,alpha=0.5,color='b')
    plt.axhline(y=power_required_4,lw=0.8,alpha=0.5,color='k')
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.xlim([621,623])
    plt.ylim([0,10])
    plt.ylabel('Leahy-normalized power',fontsize=12)
    plt.title('Energy range: PI = ' + str(PI1) + ' - ' + str(PI2) + ', W = '+ str(W) + ', Threshold = '+str(threshold) + '%' + '\n' + 'Segment Length: ' + str(segment_length) + 's, No. Segments = ' + str(M) + '\n' + 'Demodulated: ' + str(demod),fontsize=12)
    plt.legend(('Power Spectrum','3 sigma','4 sigma'),loc='best')

#plt.figure(2)
#plt.semilogy(ps_bins,N_greaterthanP,'rx')
#plt.xlabel('Leahy-normalized power',fontsize=12)
#plt.ylabel('log[N(>P)]',fontsize=12)
#plt.title('Energy range: ' + str(PI1) + ' - ' + str(PI2) + ', W = ' + str(W),fontsize=12)

    pdfname = Lv0_dirs.NICERSOFT_DATADIR + '/merged_events/merged' + merged_id + '/' + str(grid[i]) + 'PB.pdf'
    plt.savefig(pdfname,dpi=900)
"""
