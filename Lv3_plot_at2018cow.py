#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 4:41pm 2019

Plotting f vs t and \dot{f} vs t and \dot{f} vs f for AT2018cow.
Results from 6/20 - Thurs - Evernote!

"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
import Lv0_dirs

Lv0_dirs.global_par()
obsids = ['12002501' + str(i+1).zfill(2) for i in range(26)]
gti_start_end = []
for i in range(len(obsids)):
    event = Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/ni' + obsids[i] + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data
    gti_start_end.append((gtis[0][0],gtis[-1][1]))

# time values to use for when an energy range is specified, this is because some of the
# observations span about a day, so to get an approximate time value to use for plotting,
# just use what's in between...
# list will be length 26, where each entry = time to use!
times_for_Erange = [(gti_start_end[i][0]+gti_start_end[i][1])/2 for i in range(len(gti_start_end))]

pi_30_250_t = [times_for_Erange[0],times_for_Erange[13],times_for_Erange[19],times_for_Erange[23]]
pi_30_250_f = [221.85695,237.2005,206.6032,263.1900]
pi_30_250_dotf = [3.7e-7,-2.0e-5,3.24e-5,6.0e-5]

pi_250_500_t = [times_for_Erange[23],times_for_Erange[25],times_for_Erange[25]]
pi_250_500_f = [209.5328,233.1766,205.1463]
pi_250_500_dotf = [1.99e-5,-4.8e-5,-5.5e-5]

pi_500_800_t = [times_for_Erange[6],times_for_Erange[14],times_for_Erange[15]]
pi_500_800_f = [241.7417,249.5672,245.9188]
pi_500_800_dotf = [-1.74e-5,9e-7,4.6e-5]

pi_800_1200_t = [times_for_Erange[2],times_for_Erange[7],times_for_Erange[16],times_for_Erange[23]]
pi_800_1200_f = [248.091165,246.099788,241.105972,244.0893]
pi_800_1200_dotf = [-1.3e-8,-8.3e-9,-4.5e-9,-1.06e-5]

all_t = np.array([times_for_Erange[1]+2*200,times_for_Erange[1]+500,times_for_Erange[3]+28*200,times_for_Erange[6]+2*200,
        times_for_Erange[6]+2*200,times_for_Erange[6]+500,times_for_Erange[7]+5*200,times_for_Erange[7]+28*200,
        times_for_Erange[7]+30*200,times_for_Erange[7]+85*200,times_for_Erange[7]+11*500,times_for_Erange[7]+12*500,
        times_for_Erange[7]+12*500,times_for_Erange[7]+0*1000,times_for_Erange[8]+7*200,times_for_Erange[8]+120*200,
        times_for_Erange[8]+47*500,times_for_Erange[9]+56*200,times_for_Erange[16]+0*200,times_for_Erange[16]+0*200,
        times_for_Erange[18]+28*200,times_for_Erange[22]+3*200,times_for_Erange[22]+57*200,times_for_Erange[22]+23*500,
        times_for_Erange[24]+1*200,times_for_Erange[24]+169*200,times_for_Erange[24]+280*200])
all_f = np.array([255.869,238.9078,246.9331,245.003,207.345,240.5733,232.2013,230.534,211.2431,202.464,224.3295,
                206.5125,201.0860,205.8224,222.1344,224.0169,214.0273,210.1838,217.8388,215.670,237.5006,205.7663,
                248.2513,240.0588,248.5750,249.2125,228.3013])
all_dotf = np.array([-0.00050,-0.000100,0.00048,-0.00225,0.00045,8.2e-5,-0.00047,-0.00062,0.00021,-0.00080,-7.6e-5,
                    -3.2e-5,-2.4e-5,2.3e-5,0.00033,-0.00039,-3.2e-5,0.00030,6e-5,-0.00055,0.00045,-0.00010,0.00024,
                    8.8e-5,-0.00024,-0.00059,0.00054])

index_sort = np.argsort(all_t)
all_t_sorted = all_t[index_sort]
all_f_sorted = all_f[index_sort]
all_fdot_sorted = all_dotf[index_sort]

with PdfPages('fdot_f_sequential.pdf') as pdf:
    plt.figure(figsize=(16,9))
    plt.plot(pi_30_250_f,pi_30_250_dotf,'bx-')
    plt.plot(pi_250_500_f,pi_250_500_dotf,'gx-')
    plt.plot(pi_500_800_f,pi_500_800_dotf,'mx-')
    plt.plot(pi_800_1200_f,pi_800_1200_dotf,'rx-')
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Frequency Derivative (Hz/s)',fontsize=12)
    plt.ylim([-0.0025,0.0006])
    plt.xlim([200,260])
    for i in range(len(all_f_sorted)):
        plt.plot(all_f_sorted[i],all_fdot_sorted[i],'rx-')
        pdf.savefig()

    plt.legend(('0.3-2.5 keV','2.5-5 keV','5-8 keV','8-12 keV','0.3-12 keV'),fontsize=12)
    plt.close()

plt.figure(1)
#plt.plot(pi_30_250_t,pi_30_250_f,'bx-')
#plt.plot(pi_250_500_t,pi_250_500_f,'gx-')
#plt.plot(pi_500_800_t,pi_500_800_f,'mx-')
#plt.plot(pi_800_1200_t,pi_800_1200_f,'rx-')
plt.plot(all_t_sorted,all_f_sorted,'kx-')
plt.plot(all_t_sorted[3],all_f_sorted[3],'rx-')
plt.plot(all_t_sorted[4],all_f_sorted[4],'rx-')
plt.plot(all_t_sorted[11],all_f_sorted[11],'rx-')
plt.plot(all_t_sorted[12],all_f_sorted[12],'rx-')
plt.plot(all_t_sorted[18],all_f_sorted[18],'rx-')
plt.plot(all_t_sorted[19],all_f_sorted[19],'rx-')
plt.ylim([200,260])
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Frequency (Hz)',fontsize=12)
plt.title('0.3-12 keV data', fontsize=12)
#plt.legend(('0.3-2.5 keV','2.5-5 keV','5-8 keV','8-12 keV','0.3-12 keV'),fontsize=12)

plt.figure(2)
plt.plot(pi_30_250_t,pi_30_250_dotf,'bx-')
plt.plot(pi_250_500_t,pi_250_500_dotf,'gx-')
plt.plot(pi_500_800_t,pi_500_800_dotf,'mx-')
plt.plot(pi_800_1200_t,pi_800_1200_dotf,'rx-')
plt.plot(all_t_sorted,all_fdot_sorted,'kx-')
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Frequency Derivative (Hz/s)',fontsize=12)
plt.legend(('0.3-2.5 keV','2.5-5 keV','5-8 keV','8-12 keV','0.3-12 keV'),fontsize=12)

plt.figure(3)
plt.plot(pi_30_250_f,pi_30_250_dotf,'bx-')
plt.plot(pi_250_500_f,pi_250_500_dotf,'gx-')
plt.plot(pi_500_800_f,pi_500_800_dotf,'mx-')
plt.plot(pi_800_1200_f,pi_800_1200_dotf,'rx-')
plt.plot(all_f_sorted,all_fdot_sorted,'kx-')
plt.xlabel('Frequency (Hz)',fontsize=12)
plt.ylabel('Frequency Derivative (Hz/s)',fontsize=12)
plt.legend(('0.3-2.5 keV','2.5-5 keV','5-8 keV','8-12 keV','0.3-12 keV'),fontsize=12)

plt.show()
