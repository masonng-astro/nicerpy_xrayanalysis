#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 1:44pm 2019

Parseval's theorem test. Based on the meeting with Deepto on 5/24

"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import Lv2_phase

### Constructing the time series
freq1 = 0.01 #0.01 Hz
times1 = np.linspace(0,32768,32770)
sin1 = np.sin(2*np.pi*freq1*times1)
timeseries1 = sin1 + np.random.normal(0,0.2,len(times1))
zeromean1 = timeseries1 - np.mean(timeseries1)
timeseries_power1 = sum(zeromean1**2)

### Doing the FFT
fft1 = np.abs(np.fft.fft(zeromean1))**2
freqs1 = np.fft.fftfreq(zeromean1.size,times1[1]-times1[0])

### 1/N * sum of squares of FFT powers
fft_power1 = sum(fft1)/len(fft1)

print(fft_power1/timeseries_power1)

### Doing the pulse profile
phases1, phasebins1, summed_prof1 = Lv2_phase.pulse_profile(freq1,times1,sin1,0,101)

### Doing FFT of the pulse profile
fft_phase1 = np.abs(np.fft.fft(summed_prof1))**2
freqs_phase1 = np.fft.fftfreq(summed_prof1.size,(phasebins1[1]-phasebins1[0])*1/freq1)


plt.figure(1)
plt.plot(freqs1,fft1)
#plt.plot(times1,sin1,'r')
plt.figure(2)
plt.step(phasebins1[:-1],summed_prof1,'r')
plt.figure(3)
plt.hist(np.log10(fft1),bins=100,log=True)

#plt.plot(times1,sin1,'b-')
#plt.plot(times1,timeseries1,'rx')


#plt.figure(1)
#plt.plot(freqs_phase1,fft_phase1,'r')

################################################################################
################################################################################
################################################################################

### Constructing the time series
freq2 = 0.05 #0.01 Hz
times2 = np.linspace(0,32768,65540)
sin2 = np.sin(2*np.pi*freq2*times2)
timeseries2 = sin2 + np.random.normal(0,0.2,len(times2))
zeromean2 = timeseries2 - np.mean(timeseries2)
timeseries_power2 = sum(zeromean2**2)

### Doing the FFT
fft2 = np.abs(np.fft.fft(zeromean2))**2
freqs2 = np.fft.fftfreq(zeromean2.size,times2[1]-times2[0])

### 1/N * sum of squares of FFT powers
fft_power2 = sum(fft2)/len(fft2)

print(fft_power2/timeseries_power2)

### Doing the pulse profile
phases2, phasebins2, summed_prof2 = Lv2_phase.pulse_profile(freq2,times2,sin2,0,101)

### Doing FFT of the pulse profile
fft_phase2 = np.abs(np.fft.fft(summed_prof2))**2
freqs_phase2 = np.fft.fftfreq(summed_prof2.size,(phasebins2[1]-phasebins2[0])*1/freq2)


plt.figure(4)
plt.plot(freqs2,fft2)
#plt.plot(times2,sin2,'b')
plt.figure(5)
plt.step(phasebins2[:-1],summed_prof2,'b')
plt.figure(6)
plt.hist(np.log10(fft2),bins=100,log=True)

#plt.figure(2)
#plt.plot(freqs_phase2,fft_phase2,'b')

plt.show()
