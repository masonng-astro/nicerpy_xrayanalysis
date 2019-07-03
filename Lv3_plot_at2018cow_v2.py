#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 6:53pm 2019

Plotting f vs t and \dot{f} vs t and \dot{f} vs f for AT2018cow (v2).
Results from 6/29 - Sat - Evernote!

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

################################################################################
all_t_256 = np.array([gti_start_end[1][0]+(2*256+1/2*256),
                  gti_start_end[7][0]+(24*256+1/2*256),
                  gti_start_end[7][0]+(26*256+1/2*256),
                  gti_start_end[9][0]+(217*256+1/2*256),
                  gti_start_end[12][0]+(0*256+1/2*256),
                  gti_start_end[20][0]+(0*256+1/2*256)])
all_f_256 = np.array([245.021,243.652,236.298,241.826,235.711,209.668])
all_f_256_error = np.array([0.004,0.004,0.004,0.004,0.004,0.004])

all_dotf_256 = np.array([6.9e-4,2.6e-4,-2.1e-4,8e-5,-1.8e-4,-1.7e-4])
all_dotf_256_error = np.array([0.2e-4,0.2e-4,0.2e-4,2e-5,0.2e-4,0.2e-4])
################################################################################
all_t_128 = np.array([gti_start_end[1][0]+(126*128+1/2*128),
                    gti_start_end[7][0]+(44*128+1/2*128),
                    gti_start_end[7][0]+(50*128+1/2*128),
                    gti_start_end[7][0]+(94*128+1/2*128),
                    gti_start_end[7][0]+(135*128+1/2*128),
                    gti_start_end[7][0]+(137*128+1/2*128),
                    gti_start_end[7][0]+(141*128+1/2*128),
                    gti_start_end[8][0]+(187*128+1/2*128),
                    gti_start_end[12][0]+(2*128+1/2*128),
                    gti_start_end[13][0]+(8*128+1/2*128),
                    gti_start_end[20][0]+(1*128+1/2*128),
                    gti_start_end[20][0]+(45*128+1/2*128),
                    gti_start_end[21][0]+(44*128+1/2*128),
                    gti_start_end[22][0]+(90*128+1/2*128)])
all_f_128 = np.array([210.789,230.537,241.512,234.387,207.354,230.820,
                    247.117,215.590,231.939,216.830,209.656,215.895,
                    248.277,235.725])
all_f_128_error = np.array([0.008]*len(all_f_128))

all_dotf_128 = np.array([-7.3E-4,-5.5E-4,-2.99E-3,-1.34E-3,-2.75E-3,
                        -1.65E-3,-1.34E-3,-5.0E-3,1.8E-4,7.3E-4,
                        -1.8E-4,1.71E-3,1.8E-4,-1.77E-3])
all_dotf_128_error = np.array([0.6E-4,0.6E-4,0.06E-3,0.06E-3,0.06E-3,
                            0.06E-3,0.06E-3,0.1E-3,0.6E-4,0.6E-4,
                            0.6E-4,0.06E-3,0.6E-4,0.06E-3])

acc_128 = np.array([-1042,-714,-3712,-1717,-3971,-2140,-1629,-6960,237,
                    1013,-262,2373,221,-2251])
acc_128_error = np.array([87,79,76,78,88,79,74,170,79,84,87,85,74,78])
################################################################################
all_t_512 = np.array([gti_start_end[7][0]+(22*512+1/2*512),
                    gti_start_end[8][0]+(46*512+1/2*512),
                    gti_start_end[8][0]+(46*512+1/2*512)])
all_f_512 = np.array([210.804,212.007,214.025])
all_f_512_error = np.array([0.002]*len(all_f_512))

all_dotf_512 = np.array([5.0E-5,-2.06E-4,-3.1E-5])
all_dotf_512_error = np.array([0.4E-5,0.08E-4,0.8E-5])
################################################################################
all_t_64 = np.array([gti_start_end[1][0]+(251*64+1/2*64),
                    gti_start_end[4][0]+(8*64+1/2*64),
                    gti_start_end[7][0]+(96*64+1/2*64),
                    gti_start_end[7][0]+(110*64+1/2*64),
                    gti_start_end[7][0]+(274*64+1/2*64),
                    gti_start_end[8][0]+(1*64+1/2*64),
                    gti_start_end[8][0]+(26*64+1/2*64),
                    gti_start_end[8][0]+(363*64+1/2*64),
                    gti_start_end[9][0]+(873*64+1/2*64),
                    gti_start_end[9][0]+(877*64+1/2*64),
                    gti_start_end[11][0]+(262*64+1/2*64),
                    gti_start_end[15][0]+(2*64+1/2*64),
                    gti_start_end[18][0]+(3*64+1/2*64),
                    gti_start_end[20][0]+(89*64+1/2*64),
                    gti_start_end[22][0]+(4*64+1/2*64),
                    gti_start_end[22][0]+(10*64+1/2*64),
                    gti_start_end[23][0]+(9*64+1/2*64),
                    gti_start_end[24][0]+(261*64+1/2*64),
                    gti_start_end[24][0]+(266*64+1/2*64),
                    gti_start_end[24][0]+(266*64+1/2*64),
                    gti_start_end[24][0]+(529*64+1/2*64),
                    gti_start_end[24][0]+(616*64+1/2*64),
                    gti_start_end[25][0]+(4*64+1/2*64)])
all_f_64 = np.array([237.273,217.094,203.457,245.063,230.875,213.297,240.617,
                    235.344,209.359,200.301,209.328,210.070,230.414,245.621,
                    227.008,242.680,215.211,227.773,247.867,218.852,205.723,
                    241.063,207.672])
all_f_64_error = np.array([0.016]*len(all_f_64))

all_dotf_64 = np.array([-2.25E-2,-1.00E-2,5.1E-3,1.61E-2,-1.5E-3,1.17E-2,
                        -1.32E-2,-9.0E-3,-2.2E-3,-1.2E-3,-1.12E-2,-3.4E-3,
                        5.4E-3,8.8E-3,5.9E-3,-1.17E-2,5E-4,8.8E-3,
                        1.22E-2,-6.8E-3,9.3E-3,6.3E-3,-1.05E-2])
all_dotf_64_error = np.array([0.05E-2,0.02E-2,0.2E-3,0.05E-2,0.2E-3,0.02E-2,
                            0.05E-2,0.2E-3,0.2E-3,0.2E-3,0.05E-2,0.2E-3,0.2E-3,
                            0.2E-3,0.2E-3,0.02E-2,2E-4,0.2E-3,0.02E-2,0.2E-3,
                            0.2E-3,0.2E-3,0.02E-2])
acc_64 = np.array([-28380,-13820,7550,19710,-1900,16470,-16430,-11510,
                -3150,-1830,-16080,-4880,6990,10730,7740,-14480,680,11570,
                14760,-9360,13520,7890,-15150])
acc_64_error = np.array([620,340,360,600,320,340,610,310,350,370,700,350,320,
                        300,320,300,340,320,300,330,360,300,350])

print('METs for 64s segments: ')
for i in range(len(all_t_64)):
    print(all_t_64[i])
print('METs for 128s segments: ')
for i in range(len(all_t_128)):
    print(all_t_128[i])
print('METs for 256s segments: ')
for i in range(len(all_t_256)):
    print(all_t_256[i])
print('METs for 512s segments: ')
for i in range(len(all_t_512)):
    print(all_t_512[i])

################################################################################
index_sort_64 = np.argsort(all_t_64)
all_t_64_sorted = all_t_64[index_sort_64]
all_f_64_sorted = all_f_64[index_sort_64]
all_fdot_64_sorted = all_dotf_64[index_sort_64]
acc_64_sorted = acc_64[index_sort_64]
acc_64_error_sorted = acc_64_error[index_sort_64]
################################################################################
index_sort_128 = np.argsort(all_t_128)
all_t_128_sorted = all_t_128[index_sort_128]
all_f_128_sorted = all_f_128[index_sort_128]
all_fdot_128_sorted = all_dotf_128[index_sort_128]
acc_128_sorted = acc_128[index_sort_128]
acc_128_error_sorted = acc_128_error[index_sort_128]
################################################################################
index_sort_256 = np.argsort(all_t_256)
all_t_256_sorted = all_t_256[index_sort_256]
all_f_256_sorted = all_f_256[index_sort_256]
all_fdot_256_sorted = all_dotf_256[index_sort_256]
################################################################################
index_sort_512 = np.argsort(all_t_512)
all_t_512_sorted = all_t_512[index_sort_512]
all_f_512_sorted = all_f_512[index_sort_512]
all_fdot_512_sorted = all_dotf_512[index_sort_512]

"""
with PdfPages('fdot_f_sequential.pdf') as pdf:
    plt.figure(figsize=(16,9))
    plt.xlabel('Frequency (Hz)',fontsize=12)
    plt.ylabel('Frequency Derivative (Hz/s)',fontsize=12)
    plt.ylim([-0.00025,0.0007])
    plt.xlim([200,260])
    for i in range(len(all_f_sorted)):
        plt.errorbar(all_f_sorted[i],all_fdot_sorted[i],all_dotf_error[i],all_f_error[i],fmt='rx-')
        pdf.savefig()
    plt.close()
"""

plt.figure(1)
plt.axvline(x=gti_start_end[24][0]+(266*64+1/2*64),lw=0.2,alpha=0.5)
plt.errorbar(all_t_64_sorted,all_f_64_sorted,all_f_64_error,fmt='mx-')
plt.errorbar(all_t_128_sorted,all_f_128_sorted,all_f_128_error,fmt='kx-')
plt.errorbar(all_t_256_sorted,all_f_256_sorted,all_f_256_error,fmt='rx-')
plt.errorbar(all_t_512_sorted,all_f_512_sorted,all_f_512_error,fmt='bx-')
plt.ylim([200,260])
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Frequency (Hz)',fontsize=12)
plt.title('0.3-12 keV data', fontsize=12)
plt.legend(('TwoCand','64s','128s','256s','512s'),loc='best')

plt.figure(2)
plt.axvline(x=gti_start_end[24][0]+(266*64+1/2*64),lw=0.2,alpha=0.5)
plt.errorbar(all_t_64_sorted,all_fdot_64_sorted,all_dotf_64_error,fmt='mx-')
plt.errorbar(all_t_128_sorted,all_fdot_128_sorted,all_dotf_128_error,fmt='kx-')
plt.errorbar(all_t_256_sorted,all_fdot_256_sorted,all_dotf_256_error,fmt='rx-')
plt.errorbar(all_t_512_sorted,all_fdot_512_sorted,all_dotf_512_error,fmt='bx-')
#plt.ylim([-5E-3,2E-3])
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Frequency Derivative (Hz/s)',fontsize=12)
plt.legend(('TwoCand','64s','128s','256s','512s'),loc='best')

plt.figure(3)
plt.errorbar(all_f_64_sorted,all_fdot_64_sorted,all_dotf_64_error,all_f_64_error,fmt='mx-')
plt.errorbar(all_f_128_sorted,all_fdot_128_sorted,all_dotf_128_error,all_f_128_error,fmt='kx-')
plt.errorbar(all_f_256_sorted,all_fdot_256_sorted,all_dotf_256_error,all_f_256_error,fmt='rx-')
plt.errorbar(all_f_512_sorted,all_fdot_512_sorted,all_dotf_512_error,all_f_512_error,fmt='bx-')
plt.xlabel('Frequency (Hz)',fontsize=12)
plt.ylabel('Frequency Derivative (Hz/s)',fontsize=12)
plt.legend(('64s','128s','256s','512s'),loc='best')

plt.figure(4)
plt.axvline(x=gti_start_end[24][0]+(266*64+1/2*64),lw=0.2,alpha=0.5)
plt.errorbar(x=all_t_64_sorted,y=acc_64_sorted,yerr=acc_64_error_sorted,fmt='mx-')
plt.errorbar(x=all_t_128_sorted,y=acc_128_sorted,yerr=acc_128_error_sorted,fmt='kx-')
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Acceleration (m/s^2)',fontsize=12)
plt.legend(('TwoCand','64s','128s'),loc='best')

plt.show()
