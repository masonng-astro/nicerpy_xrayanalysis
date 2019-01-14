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

print(len(truncated_t))
print(len(t_bins))

#for i in range(len(counts)):
#    print(i,summed_data[i]-truncated_counts[i])

timeend = time.time()

print(timeend-timestart)
