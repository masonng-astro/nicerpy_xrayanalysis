#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Jan 5th

@author: masonng


"""
from __future__ import division, print_function
import datetime
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
import Lv0_call_eventcl,Lv1_data_bin
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy import signal

timestart = time.time()

print(datetime.datetime.now())

class Complex:
     def __init__(self, realpart, imagpart):
         self.r = realpart
         self.i = imagpart

x = Complex(4.0,3.0)
print(x.r,x.i)

obsid = '0034070101'
bary = True
par_list = ['PI','TIME','PI_FAST']
tbin_size=1
t1 = 0
t2 = 809

truncated_t, truncated_counts = Lv1_data_bin.binning_t(obsid,bary,par_list,tbin_size,t1,t2)

data_dict = Lv0_call_eventcl.get_eventcl(obsid,bary,par_list)

times = data_dict['TIME']
counts = np.ones(len(times))

shifted_t = times-times[0]
t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])*1/tbin_size+1)
summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the time values in the data

#for i in range(len(counts)):
#    print(i,summed_data[i]-truncated_counts[i])
"""
import Lv0_call_mkf
datadict = Lv0_call_mkf.get_mkf('1034070102',['TIME','ANG_DIST'])
times = datadict['TIME']
shifted_t = times-times[0]
angdist = datadict['ANG_DIST']
plt.plot(shifted_t,angdist)
plt.xlim([0,625])
plt.show()
"""

obsids = ['0034070101','0034070102','0034070103','0034070104','1034070101','1034070102','1034070103','1034070104','1034070105','1034070106']
"""
import Lv0_call_ufa
from scipy.stats import norm
for i in range(len(obsids)):
    plt.figure(1)
    datadict = Lv0_call_ufa.get_ufa(obsids[i],'7',['TIME','PI'])
    pis = datadict['PI']
    times = datadict['TIME']
    pis_trunc = pis[(pis>=5)&(pis<=15)]
    plt.hist(pis_trunc,bins=20,range=(5,15),density=True)
    plt.title(obsids[i],fontsize=15)

    mu,std = norm.fit(pis_trunc)
    print(mu)
    x = np.linspace(5,15,1001)
    p = norm.pdf(x,mu,std)
    plt.plot(x,p,'k')

    plt.show()
"""
#plt.figure(2)
#event = '/Users/masonng/Documents/MIT/Research/ni1034070106_0mpu7_ufa.evt'
#event = fits.open(event)

#pis = event[1].data['PI']
#plt.hist(pis,bins=20,range=(5,15))

workingdir = '/Users/masonng/Documents/MIT/Research/nicer-data/0034070101_test_xsel/xti/event_cl/'
fitsfile = 'trymodel1.fits'

E,c1,c2,c3,c4 = np.loadtxt(workingdir+fitsfile,skiprows=3,usecols=(0,1,2,3,4),unpack=True)
plt.semilogy(E,c2,'r-')
plt.semilogy(E,c3,'k-')
plt.semilogy(E,c4,'c-')
plt.ylim([1e-4,1])

plt.show()

timeend = time.time()

print(timeend-timestart)
