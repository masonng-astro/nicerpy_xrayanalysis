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
import matplotlib.gridspec as gridspec
import Lv0_dirs,Lv0_call_eventcl,Lv1_data_bin
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy import signal

timestart = time.time()

Lv0_dirs.global_par()

print(datetime.datetime.now())

"""

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

"""

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
"""
E,c1,c2,c3,c4 = np.loadtxt(workingdir+fitsfile,skiprows=3,usecols=(0,1,2,3,4),unpack=True)
plt.semilogy(E,c2,'r-')
plt.semilogy(E,c3,'k-')
plt.semilogy(E,c4,'c-')
plt.ylim([1e-4,1])

plt.show()
"""
"""
import xspec
tryfile = '/Users/masonng/Documents/MIT/Research/nicer-data/0034070101/xti/event_cl/0034070101_spec.pi'
#tryfile = '0034070101_spec.pi'
s = xspec.Spectrum(tryfile)
"""

"""
renorm = '/Users/masonng/Documents/MIT/Research/nicer-data/crabcor.dat'
test_spec = '/Users/masonng/Documents/MIT/Research/nicer-data/0034070101/xti/event_cl/0034070101_spec.pi'
test_renormspec = '/Users/masonng/Documents/MIT/Research/nicer-data/0034070101/xti/event_cl/0034070101_renormspec.pi'

event = fits.open(test_spec)
event_renorm = fits.open(test_renormspec)

eventchan = event[1].data['CHANNEL']*10/1000
eventcount = event[1].data['COUNTS']

chan,ratio = np.loadtxt(renorm,usecols=(0,1),unpack=True)

renormchan = event_renorm[1].data['CHANNEL']*10/1000
renormcount = event_renorm[1].data['COUNTS']

#for i in range(100,500):
#    print(chan[i],eventcount[i],ratio[i],renormcount[i])

#plt.plot(chan,ratio)
#plt.plot(eventchan[100:500],eventcount[100:500],'rx')
#plt.plot(renormchan[100:500],renormcount[100:500],'bx')
#plt.plot(eventchan,eventcount,'rx')
#plt.plot(renormchan,renormcount,'bx')
#plt.legend(('Original','Renorm'),loc='best')
#plt.show()
"""

"""
event = Lv0_dirs.NICER_DATADIR + '0034070101/xti/event_cl/0034070101_truncated.pi'
event = fits.open(event)

counts = event[1].data['COUNTS']
print(sum(counts))

datadict = Lv0_call_eventcl.get_eventcl('0034070101',True,['PI','PI_FAST','TIME'])
times = datadict['TIME']
pi = datadict['PI']
T = times[-1]-times[0]
counts = np.ones(len(times))
print(min(pi),max(pi))
print(sum(counts))
print(sum(counts)/T)
print(T)
"""

"""
datadict = Lv0_call_eventcl.get_eventcl('1034070104',True,['PI','PI_FAST','TIME','EVENT_FLAGS'])
times = datadict['TIME']
flags = datadict['EVENT_FLAGS']
counts = np.ones(len(times))

shifted_t = times-times[0]
times_trunc = times[(shifted_t>=11700)&(shifted_t<=11900)] #55281
flags_trunc = flags[(shifted_t>=11700)&(shifted_t<=11900)]
t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1]+1))
#print(len(times_trunc)) #55281
#for i in range(23200,23500):
#    print(shifted_t[i],flags_trunc[i])

binned_counts, edges, binno = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins)
#plt.plot(t_bins[:-1],binned_counts,'r-')
#plt.show()
"""


x = [np.linspace(0,100,10001),np.linspace(100,200,10001),np.linspace(200,300,10001),np.linspace(300,400,10001),np.linspace(400,500,10001),np.linspace(500,600,10001)]
y = [4*x[0]+34,4*x[1]+34,4*x[2]+34,4*x[3]+34,4*x[4]+34,4*x[5]+34]
"""
fig,axes = plt.subplots(1,1,sharex=True)
axes[0,0].plot(x,y)

plt.show()
"""

"""
#Creates four polar axes, and accesses them through the returned array
fig,(p10,p20,p30,p40,p50,p60) = plt.subplots(6,1,figsize=(12.8,4.0))
gs = gridspec.GridSpec(6,1)
p10 = plt.subplot(gs[0])
p20 = plt.subplot(gs[1])
p30 = plt.subplot(gs[2])
p40 = plt.subplot(gs[3])
p50 = plt.subplot(gs[4])
p60 = plt.subplot(gs[5])
plotting = [p10,p20,p30,p40,p50,p60]

for i in range(6):
    plt.subplot(gs[i]).plot(x[i],y[i])

plt.show()

"""

"""
import subprocess
import os
#subprocess.check_output(['/Users/masonng/Documents/MIT/Research/Nicer-Logistics/heasoft-6.25/x86_64-apple-darwin17.7.0/headas-init.sh'],shell=True)
#subprocess.call(['/bin/bash','-i','-c','heainit','nicerversion'],shell=True)
subprocess.check_call('nicerversion')#,env=os.environ)
timeend = time.time()
"""

"""
import subprocess
from astropy.io.fits import getdata,update

#event = '/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe/GTI_1000s/GTI1_breakup/ni1034090111_nicersoft_bary_GTI1_1000s.evt'
#newevent = '/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe/GTI_1000s/GTI1_breakup/ni1034090111_nicersoft_bary_GTI1_1000s_short.evt'
#subprocess.check_call(['cp',event,newevent])
#new_list = np.rec.array([(0,1000)])
event = '/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe/ni1034090111_nicersoft_bary_stitched.evt'
data,hdr = getdata(event,1,header=True)
times = data['TIME']
counts = np.ones(len(times))
shifted_t = times-times[0]
t_bins = np.linspace(0,shifted_t[-1],len(shifted_t)+1)
summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins)

plt.plot(t_bins[:-1],summed_data,'rx')
plt.show()
#update(newevent,new_list,2,header=hdr)
"""
#header = fits.Header()

"""
eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe_old/GTI_1000s/ni1034090111_nicersoft_bary_GTI1_1000s.dat'
contents = np.fromfile(eventfile,dtype='<f',count=-1)

print(len(contents),type(contents))
print(np.max(contents))

"""

"""
import subprocess
subprocess.check_call(['nicerversion'])
import sh
import glob

dir = '/Volumes/Samsung_T5/1060060127/xti/hk/'
hk_files = glob.glob(dir+'*.gz')

for i in range(len(hk_files)):
    sh.gunzip(hk_files[i])

#subprocess.check_call(['scp', '-r', 'masonng@ciri:/nfs/ciri/nicer/decrypted/1060060127','/Volumes/Samsung_T5'])
#subprocess.check_call(['gunzip','/Volumes/Samsung_T5/1060060127/xti/hk/*.gz'])
#subprocess.check_call(['psrpipe.py'])
#subprocess.check_call(['barycorr'])
#subprocess.check_call(['nicerfits2presto.py'])
#subprocess.check_call(['realfft'])
#subprocess.check_call(['accelsearch'])
#subprocess.check_call(['prepfold'])
"""

import subprocess

infile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe_old/GTI_1000s/ni1034090111_nicersoft_bary_GTI1_1000s_63.86Hz_Cand.pfd.ps'
outfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe_old/GTI_1000s/ni1034090111_nicersoft_bary_GTI1_1000s_63.86Hz_Cand.pfd.pdf'
subprocess.check_call(['ps2pdf',infile,outfile])


timeend = time.time()

print(timeend-timestart)
