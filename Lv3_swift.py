#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 15 5:33pm 2019

Plotting up the daily light curve from Swift's Cen X-3 data.

Recall that Swift's Burst Alert Telescope is 15-150 keV.
https://swift.gsfc.nasa.gov/about_swift/bat_desc.html

Data array is in the form:

TIME (days)
RATE (count/cm^2/s)
ERROR (count/cm^2/s)
YEAR (yr)
DAY (d) (corresponding to YEAR)
STAT_ERR (count/cm^2/s)
SYS_ERR (count/cm^2/s)
DATA_FLAG
TIMEDEL_EXPO
TIMEDEL_CODED
TIMEDEL_DITHERED

Recall that [Converter: http://www.csgnetwork.com/julianmodifdateconv.html]

0034070101: 06/20/2017 - MJD: 57924 - flux: 0.0308 +- 0.0017
0034070102: 06/21/2017 - MJD: 57925 - flux: 0.0408 +- 0.0019
0034070103: 06/22/2017 - MJD: 57926 - flux: 0.0235 +- 0.0021
0034070104: 06/23/2017 - MJD: 57927 - flux: 0.0420 +- 0.0021
1034070101: 07/24/2017 - MJD: 57958 - flux: 0.0249 +- 0.0014
1034070102: 07/25/2017 - MJD: 57959 - flux: 0.0080 +- 0.0012
1034070103: 07/26/2017 - MJD: 57960 - flux: 0.0165 +- 0.0012
1034070104: 07/27/2017 - MJD: 57961 - flux: 0.0056 +- 0.0012
1034070105: 07/28/2017 - MJD: 57962 - flux: 0.0087 +- 0.0011
1034070106: 11/07/2018 - MJD: 58429 - flux: 0.0059 +- 0.0111
"""
from __future__ import division, print_function
import numpy as np
import time
import matplotlib.pyplot as plt
import Lv0_dirs,Lv2_phase
from scipy import signal

Lv0_dirs.global_par()

comm_phase_t = np.array([57924,57925,57926,57927])
comm_phase_flux = np.array([0.0308,0.0408,0.0235,0.0420])

sci_phase_t = np.array([57958,57959,57960,57961,57962])
sci_phase_flux = np.array([0.0249,0.0080,0.0165,0.0056,0.0087])

nov_phase_t = np.array([58429])
nov_phase_flux = np.array([0.0059])

swift_file = Lv0_dirs.BASE_DIR + 'CenX-3.lc.txt'
time,rate,error = np.loadtxt(swift_file,dtype='float',skiprows=5,usecols=(0,1,2),unpack=True)
#print(np.mean(rate[(time>=57754)&(time<=58118)])) #mean over 2017
#print(np.mean(rate[(time>=58119)&(time<=58483)])) #mean over 2018
"""
times = time-time[0]
T = times[-1]

plt.figure(1)
plt.plot(times,rate,'rx-')

plt.figure(2)

dt = T/len(times)
freq = 1/dt
f,pxx = signal.periodogram(rate,fs=freq)

plt.semilogx(f[1:],pxx[1:])
plt.xlabel('Hz',fontsize=12)
plt.ylabel('Normalized power spectrum',fontsize=12)

#plt.axvline(x=vlines[1],color='k',alpha=0.5,lw=0.5)
plt.figure(3)
phase_bins, pulse_profile = Lv2_phase.pulse_profile(1/(2.1*86400),time,rate,0.5,40)
plt.plot(phase_bins[:-1],pulse_profile,'rx-')

plt.show()
"""

slice = 1

plt.figure(1)
plt.errorbar(time[::slice],rate[::slice],yerr=error[::slice],fmt='b-',alpha=0.4)

plt.plot(comm_phase_t,comm_phase_flux,'rx')
plt.plot(sci_phase_t,sci_phase_flux,'rx')
plt.plot(nov_phase_t,nov_phase_flux,'rx')

plt.xlabel('Time(d)')
plt.ylabel('Count (count/cm^2/s)')
plt.axhline(y=0,lw=0.5,alpha=0.5)

plt.figure(2)
plt.errorbar(time[::slice],rate[::slice],yerr=error[::slice],fmt='bx',alpha=0.4)

plt.plot(comm_phase_t,comm_phase_flux,'rx')
plt.plot(sci_phase_t,sci_phase_flux,'rx')
plt.plot(nov_phase_t,nov_phase_flux,'rx')

plt.xlabel('Time(d)')
plt.ylabel('Count (count/cm^2/s)')
plt.axhline(y=0,lw=0.5,alpha=0.5)
plt.show()


def plot(xlim1,xlim2):
    """
    Plots the daily light curve from the Swift data for Cen X-3. Can generalize to other objects in the future.

    xlim1 - lower time boundary (in MJD)
    xlim2 - upper time boundary (in MJD)
    """

    swift_file = Lv0_dirs.BASE_DIR + 'CenX-3.lc.txt'

    time,rate,error = np.loadtxt(swift_file,dtype='float',skiprows=5,usecols=(0,1,2),unpack=True)

    plt.figure(1)
    plt.errorbar(time,rate,yerr=error,fmt='bx',alpha=0.4)
    plt.xlabel('Time(d)')
    plt.xlim([xlim1,xlim2])
    plt.ylabel('Count (count/cm^2/s)')
    plt.axhline(y=0,lw=0.5,alpha=0.5)

    plt.plot(sci_phase_t,sci_phase_flux,'rx')
    plt.show()

def plot_compare2(xlims1,xlims2):
    """
    Plots the daily light curve from the Swift data for Cen X-3. Can generalize to other objects in the future.
    This particular function is to compare TWO 'epochs' of observations.

    xlims1 - list which provides the lower time boundary (in MJD) and upper time boundary for the 1st epoch
    xlims2 - list which provides the lower time boundary (in MJD) and upper time boundary for the 2nd epoch
    """

    swift_file = Lv0_dirs.BASE_DIR + 'CenX-3.lc.txt'

    time,rate,error = np.loadtxt(swift_file,dtype='float',skiprows=5,usecols=(0,1,2),unpack=True)

    plt.figure(1)
    plt.errorbar(time,rate,yerr=error,fmt='x')
    plt.xlabel('Time(d)')
    plt.xlim([xlims1[0],xlims1[1]])
    plt.ylabel('Count (count/cm^2/s)')

    plt.figure(2)
    plt.errorbar(time,rate,yerr=error,fmt='x')
    plt.xlabel('Time(d)')
    plt.xlim([xlims2[0],xlims2[1]])
    plt.ylabel('Count (count/cm^2/s)')

#plot(57900,58050)
#plot_compare2([57920,57930],[57955,57965])

### BLIND PLOTTING:
"""
obs = np.array([1,2,3,4,5,6,7,8,9,10])
fluxes = np.array([0.0308,0.0408,0.0235,0.0420,0.0249,0.0080,0.0165,0.0056,0.0087,0.0059])
errors = np.array([0.0017,0.0019,0.0021,0.0021,0.0014,0.0012,0.0012,0.0012,0.0011,0.0111])

plt.figure(figsize=(10,8))
plt.ylabel('Count (count/cm^2/s)')
plt.title('Cen X-3 observations from Swift BAT. Annotations indicate date, MJD, and ObsID. ')
plt.errorbar(obs,fluxes,errors,fmt='x')
plt.annotate("06/20/2017 \n MJD 57924 \n 0034070101",(0.5,0.05))
plt.annotate("06/21/2017 \n MJD 57925 \n 0034070102",(1.5,0.05))
plt.annotate("06/22/2017 \n MJD 57926 \n 0034070103",(2.5,0.05))
plt.annotate("06/23/2017 \n MJD 57927 \n 0034070104",(3.5,0.05))
plt.annotate("07/24/2017 \n MJD 57958 \n 1034070101",(4.5,0.05))
plt.annotate("07/25/2017 \n MJD 57959 \n 1034070102",(5.5,0.05))
plt.annotate("07/26/2017 \n MJD 57960 \n 1034070103",(6.5,0.05))
plt.annotate("07/27/2017 \n MJD 57961 \n 1034070104",(7.5,0.05))
plt.annotate("07/28/2017 \n MJD 57962 \n 1034070105",(8.5,0.05))
plt.annotate("11/07/2018 \n MJD 58429 \n 1034070106",(9.5,0.05))
plt.xlim([0,11])
plt.ylim([-0.01,0.06])
plt.axvline(x=4.5,lw=0.5,alpha=0.5)
plt.axvline(x=9.5,lw=0.5,alpha=0.5)
plt.show()
"""
