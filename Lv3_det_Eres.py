#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 2 7:01pm 2020

Determining the energy resolution of detectors by fitting model Gaussians with the
instrument response folded in.

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
import time

input_E = np.array([0.5,0.6,0.7,0.8,0.9] + list(np.arange(1.0,12.5,0.5)))
pi_ev = input_E * 100
fwhm = []
fwhm_err = []

for i in range(len(pi_ev)):

    def gauss(x,a,sig,constant):
        return a*np.exp( -(x-input_E[i])**2/(2*sig**2) ) + constant

    #testspec = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/gauss_' + str(int(pi_ev[i])).zfill(5) + '.txt'
    testspec = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/getEres/gauss_mos1_' + str(int(pi_ev[i])).zfill(4) + '.txt'
    print(testspec)
    E,E_err,flux,flux_err,model = np.genfromtxt(testspec,usecols=(0,1,2,3,4),skip_header=3,skip_footer=1,unpack=True)

    pguess = np.array([100,0.05,0.01])
    popt,pcov = curve_fit(gauss,E,model,p0=pguess)

    print('Gaussian width is ' + str(round(popt[1]*1000,2)) + ' eV; FWHM is ' + str(round(popt[1]*1000*2.355,2)) + ' +- ' + str(round(np.sqrt(np.diag(pcov)[1])*1000*2.355,2)) + ' eV, and divided by 3, this is ' + str(round(popt[1]*1000*2.355/3,2))  + 'eV')
    fwhm.append(popt[1]*1000*2.355)
    fwhm_err.append(np.sqrt(np.diag(pcov)[1])*1000*2.355)
    #plt.plot(E,model,'b-')
    #plt.plot(E,gauss(E,popt[0],popt[1],popt[2]),'r-')
    #plt.show()

timeend = time.time()


def sqrtfunc(x,a,b,c):
    return a*(x+b)**c

popt,pcov = curve_fit(sqrtfunc,input_E,fwhm,sigma=fwhm_err,p0=[10,3.3,0.5])
print(popt)
print(np.sqrt(np.diag(pcov)))
#for i in range(len(input_E)):
#    print(input_E[i],sqrtfunc(input_E,popt[0],popt[1],popt[2],popt[3])[i]/3)

N = 8.7
fano = 0.114
w = 3.71

nicer_spie = 2.35*w*(N**2 + input_E*1000*fano/w)**(1/2) #for NICER

plt.plot(input_E,fwhm,'rx-')
plt.plot(input_E,sqrtfunc(input_E,*popt),'b-')
plt.annotate(str(round(popt[0],2))+'(E+'+str(round(popt[1],2))+')^'+str(round(popt[2],2)),(2,160))
plt.xlabel('Energy, E (keV)',fontsize=12)
plt.ylabel('FWHM (eV)',fontsize=12)
plt.legend(('Measured','curve_fit'),fontsize=12)
plt.show()
