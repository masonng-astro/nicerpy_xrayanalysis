#!/usr/bin/env python
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
#import binary_psr
import matplotlib.gridspec as gridspec
import subprocess
import glob
import Lv0_dirs,Lv0_call_eventcl,Lv1_data_bin,Lv1_data_gtis,Lv2_phase
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

#obsids = ['0034070101','0034070102','0034070103','0034070104','1034070101','1034070102','1034070103','1034070104','1034070105','1034070106']
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

#workingdir = '/Users/masonng/Documents/MIT/Research/nicer-data/0034070101_test_xsel/xti/event_cl/'
#fitsfile = 'trymodel1.fits'
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

#x = [np.linspace(0,100,10001),np.linspace(100,200,10001),np.linspace(200,300,10001),np.linspace(300,400,10001),np.linspace(400,500,10001),np.linspace(500,600,10001)]
#y = [4*x[0]+34,4*x[1]+34,4*x[2]+34,4*x[3]+34,4*x[4]+34,4*x[5]+34]

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

"""
import subprocess

infile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe_old/GTI_1000s/ni1034090111_nicersoft_bary_GTI1_1000s_63.86Hz_Cand.pfd.ps'
outfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe_old/GTI_1000s/ni1034090111_nicersoft_bary_GTI1_1000s_63.86Hz_Cand.pfd.pdf'
subprocess.check_call(['ps2pdf',infile,outfile])
"""

"""
eventfile = Lv0_dirs.NICERSOFT_DATADIR + '1060060127_pipe/test.dat'
bins = np.fromfile(eventfile,dtype='<f',count=-1)
ts = np.linspace(0,100,200000)
tbins = np.linspace(0,100,101)
sumcounts,bin_edges,binnumber = stats.binned_statistic(ts,bins,statistic='sum',bins=tbins)

plt.plot(tbins[:-1],sumcounts,'r-')
plt.show()
print(len(bins[bins>0])/len(bins)*100)

"""

"""
obsdir = Lv0_dirs.NICERSOFT_DATADIR + '0034070101_pipe/'
test_dat = obsdir + 'test.dat'

bins = np.fromfile(test_dat,dtype='<f',count=-1) #reads the binary file ; converts to little endian, count=-1 means grab everything
bins_with_data = len(bins[bins>0]) #number of bins with data (NOT the ones with averaged count rate!)
average_count_rate = sum(bins)/len(bins)

segment_length = 1000
tbin = 0.00025
print(len(bins))
no_desired_bins = segment_length/tbin #number of desired bins, to go up to the desired segment length
no_padded = int(no_desired_bins-len(bins)) #number of bins *needed* for original segment to have desired segment length
padding = np.ones(no_padded,dtype=np.float32)*average_count_rate #generate the array of (averaged) counts needed to pad the original segment
#padding = np.zeros(no_padded,dtype=np.float32)
#padding = np.zeros(no_padded,dtype=np.float32)
new_bins = np.array(list(bins) + list(padding)) #the new set of bins where it contains the original segment, in addition to the padded bins (with counts = average of original segment)
print(len(new_bins))
new_bins.tofile('test3.dat')

import subprocess
subprocess.check_output(['mv','test3.dat',obsdir])

test3_dat = obsdir + 'test3.dat'
bins = np.fromfile(test3_dat,dtype='<f',count=-1)

print(len(bins))

times = np.linspace(0,len(bins)*tbin,len(bins))
plt.plot(times,bins,'r-')
plt.show()

def pad_binary(obsid,tbin,segment_length):

    To pad the binary file so that it will be as long as the desired segment length.
    The value to pad with for each time bin, is the average count rate in THAT segment!

    obsid - Observation ID of the object of interest (10-digit str)
    tbin - size of the bins in time
    segment_length - length of the individual segments

    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsdir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    dat_files = glob.glob(obsdir+'*.dat')
    for i in range(len(dat_files)):
        bins = np.fromfile(dat_files[i],dtype='<f',count=-1) #reads the binary file ; converts to little endian, count=-1 means grab everything
        bins_with_data = len(bins[bins>0]) #number of bins with data (NOT the ones with averaged count rate!)
        average_count_rate = sum(bins)/len(bins)

        no_desired_bins = segment_length/tbin #number of desired bins, to go up to the desired segment length
        no_padded = int(no_desired_bins-len(bins)) #number of bins *needed* for original segment to have desired segment length
        padding = np.ones(no_padded)*average_count_rate #generate the array of (averaged) counts needed to pad the original segment
        new_bins = list(bins) + list(padding) #the new set of bins where it contains the original segment, in addition to the padded bins (with counts = average of original segment)

        final_dat = open(dat_files[i],"wb")
        final_dat.write(new_bins)

    return
"""

"""
tbin = 0.00025
duty_cycle_bin = 1
gtino = 0
segment_length = 1000

binned_data = np.fromfile('/Volumes/Samsung_T5/NICERsoft_outputs/1034090111_pipe/ni1034090111_nicersoft_bary_GTI'+str(gtino)+'_'+str(segment_length)+'s.dat',dtype='<f',count=-1)
dat_times = np.arange(0,tbin*len(binned_data),tbin)

duty_cycle_times = np.arange(0,tbin*len(binned_data)+tbin,duty_cycle_bin)

summed_data, binedges, binnumber = stats.binned_statistic(dat_times,binned_data,statistic='sum',bins=duty_cycle_times)
print(len(summed_data))
print(len(summed_data[(summed_data[i]!=stats.mode(summed_data))&(sum(summed_data[i:i+5])!=5*stats.mode(summed_data) for i in range(len(summed_data)))&(summed_data>0)]))
print(len(summed_data[(summed_data!=stats.mode(summed_data)[0][0])&(summed_data>0)]))
plt.plot(duty_cycle_times[:-1],summed_data,'rx')
plt.axhline(y=stats.mode(summed_data)[0][0])
plt.show()
"""

# mkgti.py --gtiname testgti.gti 109451762 109451962
# niextract-events /Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/ni0034070101_nicersoft_bary.evt testtrunc.evt timefile=testgti.gti
# nicerfits2presto.py --dt=0.00025 testtrunc.evt

"""
bins = np.fromfile('testtrunc.dat',dtype='<f',count=-1)
print(sum(bins[:800000]))
print(sum(bins[800000:]))
"""

"""
from PyAstronomy.pyasl import foldAt
import matplotlib.pylab as plt
import numpy as np

# Generate some data ...
time = np.random.random(1000) * 100.
flux = 0.05 * np.sin(time*(2.*np.pi/21.5) + 15)
# ... and add some noise
flux += np.random.normal(0, 0.02, len(flux))

# Obtain the phases with respect to some
# reference point (in this case T0=217.4)
phases = foldAt(time, 21.5, T0=217.4)

# Sort with respect to phase
# First, get the order of indices ...
sortIndi = np.argsort(phases)
# ... and, second, rearrange the arrays.
phases = list(phases[sortIndi]) + list(phases[sortIndi]+1)
flux = list(flux[sortIndi])*2

phase_bins = np.linspace(0,2,51)
summed_profile, bin_edges, binnumber = stats.binned_statistic(phases,flux,statistic='sum',bins=phase_bins)

# Plot the result
plt.figure(1)
plt.plot(phases, flux, 'bp')
plt.axvline(x=1,alpha=0.5)
plt.figure(2)
plt.plot(phase_bins[:-1],summed_profile,'r')
plt.show()
"""

"""
#### TESTING semicolons; using multiple commands at once! Need shell=True
import os
import subprocess
#subprocess.check_call(['cd','/Volumes/Samsung_T5/NICERsoft_outputs/testdirha/',';', 'mv','testtest2.txt','renamedYAY.txt',';','cd','..',';', 'mv','testtest1.txt','renamedAGAIN.txt'],shell=True)
subprocess.Popen('cd /Volumes/Samsung_T5/NICERsoft_outputs/testdirha/ ; mv testtest2.txt renamedYAY.txt ; cd .. ; mv testtest1.txt renamedAGAIN.txt',shell=True)
"""

"""
input_filename = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe_old/ni0034070101_nicersoft_bary_ACCEL_200'
input_file = open(input_filename,'r').read().split('\n')
print(len(input_file),type(input_file))
a = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
b = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

print(len(a))
print(len(b))
print(input_file[0],len(input_file[0]),input_file[0]==a)
print(input_file[17],len(input_file[17]),input_file[17]==b)
print(np.where(np.array(input_file)==a)[0][0])
print(np.where(np.array(input_file)==b)[0][0])
"""

#import subprocess
#logfile = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe_old/logfile.txt'
#log = open(logfile,'a')
#subprocess.Popen('cd /Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe_old/ ; prepfold -double -events -noxwin -n 50 -accelcand 1 -accelfile ni0034070101_nicersoft_bary_800-1200_ACCEL_0.cand ni0034070101_nicersoft_bary.events',stdout=log,shell=True)

"""
obsid = '1060060170'
base_folder = '/Volumes/Samsung_T5/NICERsoft_outputs/'+obsid+'_pipe/'
event = base_folder + 'ni'+obsid+'_nicersoft_bary.evt'
event = fits.open(event)
counts = len(event[1].data['TIME'])
gtis = event[2].data
total_gti = sum([ (gtis[i][1] - gtis[i][0]) for i in range(len(gtis))])
print('Number of counts: '+str(counts))
print(total_gti)
print('Exposure time is ' + str(gtis[-1][1]-gtis[0][0]) + ' s')
"""

"""
obsid = '1060060170'
base_folder = '/Volumes/Samsung_T5/NICERsoft_outputs/'+obsid+'_pipe/'
event = base_folder + 'ni'+obsid+'_nicersoft_bary.evt'
event = fits.open(event)
times = event[1].data['TIME']

import binary_psr
timea = binary_psr.binary_psr("/Volumes/Samsung_T5/NICERsoft_outputs/J1231-1411.par").demodulate_TOAs(times)
for i in range(100):
    print(times[i],timea[i])
"""

"""
basefolder = '/Volumes/Samsung_T5/NICERsoft_outputs/1034070104_pipe_old/'
event = basefolder + 'cleanfilt.evt'
event = fits.open(event)
times = event[1].data['TIME']

times_zero = times - times[0]
counts = np.ones(len(times))

tbins = np.linspace(0,len(times_zero),len(times_zero)+1)
summed_counts,binedges,binnumber = stats.binned_statistic(times_zero,counts,statistic='sum',bins=tbins)
plt.plot(tbins[:-1],summed_counts,'r-')
plt.title('Cen X-3 ; ObsID 1034070104',fontsize=12)
plt.xlabel('Elapsed time (s)',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.show()
"""

### Testing FFT from PRESTO and FFT from manual method
"""
cenx3_data = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/ni0034070101_nicersoft_bary.dat'
cenx3_fft = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/ni0034070101_nicersoft_bary.fft'

raw_data = np.fromfile(cenx3_data,dtype='<f',count=-1)

freqs = np.fft.fftfreq(raw_data.size,0.00025)
N = len(freqs)
use_freqs = freqs[1:int(N/2)]

new_raw_data = raw_data - np.mean(raw_data)
my_fft = np.fft.fft(new_raw_data)
my_ps = 2/sum(raw_data)*np.abs(my_fft)**2
use_my_ps = my_ps[1:int(N/2)]
print(use_my_ps[:5])

fft_data = np.fromfile(cenx3_fft,dtype='complex64',count=-1)
fft_ps = 2/sum(raw_data)*np.abs(fft_data)**2
#fft_ps = signal.detrend(fft_ps,type='constant')
use_fft_ps = fft_ps[1:int(N/2)]
print(use_fft_ps[:5])

print('Mean of my power spectrum: ' + str(np.mean(use_my_ps[use_freqs>10])))
print('Mean of PRESTO power spectrum: ' + str(np.mean(use_fft_ps[use_freqs>10])))

plt.figure(1)
plt.plot(use_freqs,use_my_ps)
plt.figure(2)
plt.plot(use_freqs,use_fft_ps)
plt.show()
"""

"""
test_data = '/Volumes/Samsung_T5/NICERsoft_outputs/1060020113_pipe/ni1060020113_nicersoft_bary.evt'
event = fits.open(test_data)
times = event[1].data['TIME']

times_zero = times - times[0]
counts = np.ones(len(times_zero))

tbin_size = 0.00025
startt = 0
endt = int(times_zero[-1])
t_bins = np.linspace(startt,endt,(endt-startt)*1/tbin_size+1)
summed_data, bin_edges, binnumber = stats.binned_statistic(times_zero,counts,statistic='sum',bins=t_bins)

no_bins = (endt-startt)/tbin_size
nicersoft_t_bins = np.arange(no_bins+1,dtype=np.float)*tbin_size
sums,edges = np.histogram(times_zero,bins=nicersoft_t_bins)
dat = np.array(sums,np.float32)

print(t_bins[20],summed_data[20])
print(nicersoft_t_bins[20],sums[20])
print('--')
print(t_bins[-20:],summed_data[-20:])
print(nicersoft_t_bins[-20:],sums[-20:])

plt.plot(t_bins[:100000],summed_data[:100000],'rx-')
plt.plot(nicersoft_t_bins[:100000],sums[:100000],'bx-')
plt.show()
"""

"""
import Lv3_average_ps_segments
#### FROM PRESTO
testdata = '/Volumes/Samsung_T5/NICERsoft_outputs/1060020113_pipe/accelsearch_1000s/ni1060020113_nicersoft_bary_GTI0_1000s.dat'
presto_bin_data = np.fromfile(testdata,dtype='<f',count=-1)
binsize = 1000/len(presto_bin_data)
print(binsize)
presto_t_bins = np.arange(0,1000,binsize)

plt.figure(1)
plt.plot(presto_t_bins,presto_bin_data,'rx-')
#print(presto_t_bins[:20],presto_bin_data[:20])
#print(presto_t_bins[-20:],presto_bin_data[-20:])
print('--')

#### FROM MY Lv3_AVERAGE_PS_SEGMENTS!
#test_data = '/Volumes/Samsung_T5/NICERsoft_outputs/1060020113_pipe/accelsearch_1000s/ni1060020113_nicersoft_bary_GTI0_1000s.evt'
test_data = '/Volumes/Samsung_T5/NICERsoft_outputs/1060020113_pipe/ni1060020113_nicersoft_bary.evt'
event = fits.open(test_data)
times = event[1].data['TIME']

times_zero = times - event[2].data[0][0]
counts = np.ones(len(times_zero))

time_bins = np.arange(0,int(times_zero[-1]),1)
summed_threshold,bin_edges,binnumber = stats.binned_statistic(times_zero,counts,statistic='sum',bins=time_bins)

#binned_t,binned_data = Lv3_average_ps_segments.binned_data('1060020113',['TIME','PI','PI_FAST'],0.00025)

#truncated_t = binned_t[(binned_t>=0)&(binned_t<=1000)]
#truncated_counts = binned_data[(binned_t>=0)&(binned_t<=1000)]


#plt.figure(2)
plt.plot(time_bins[:-1],summed_threshold,'bx-')
#print(truncated_t[:20],truncated_counts[:20])
#print(truncated_t[-20:],truncated_counts[-20:])

plt.show()
"""

"""
obsids = ['12002501' + str(i+1).zfill(2) for i in range(26)]
for i in range(len(obsids)):
    eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/' + obsids[i] + '_pipe/ni' + obsids[i] + '_nicersoft_bary.evt'
    event = fits.open(eventfile)
    gtis = event[2].data
    print('Observation duration for ' + obsids[i] + ': ' + str(gtis[-1][1]-gtis[0][0]))
"""

"""
obsid = '1200250101'
basefile = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
header_card = fits.open(basefile)[0].header
date_obs = header_card['DATE-OBS']
date_end = header_card['DATE-END']
tstart = header_card['TSTART']
tend = header_card['TSTOP']
print(date_obs,date_end,tstart,tend)
"""

"""
import Lv0_call_nicersoft_eventcl
import Lv3_average_ps_segments

data_dict = Lv0_call_nicersoft_eventcl.get_eventcl('0034070101',[True,'','','','',''],['PI','TIME'])
gtis = Lv0_call_nicersoft_eventcl.open_fits('0034070101',[True,'','','','',''])[2].data

times = data_dict['TIME']
shifted_t = times - gtis[0][0]
counts = np.ones(len(times))
startt = shifted_t[0]
endt = shifted_t[-1]
tbin_size = 0.00025

print('Binning started.')
t_bins = np.arange(startt,endt,tbin_size)
t_bins_paul = np.arange(int((endt-startt)/tbin_size)+1)*tbin_size
print(len(t_bins))
print(len(t_bins_paul))
print(t_bins[:5],t_bins_paul[:5])

summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins) #binning the counts in the data
sums,edges = np.histogram(shifted_t,bins=t_bins_paul)
print('Binning finished.')

print(summed_data[:5],summed_data[-5:],sums[:5],sums[-5:])
#print(t_bins[399995:400005])
#print(summed_data[399995:400005])

plt.figure(1)
#plt.step(t_bins[:-1],summed_data,'r-')
plt.step(t_bins_paul[:-1],sums,'r-')
print(t_bins[:5])
dat_files = Lv3_average_ps_segments.presto_dat('0034070101',100)
fft_files = Lv3_average_ps_segments.presto_FFT('0034070101',100)

counter = 0
for i in range(len(dat_files)):
    counts = np.fromfile(dat_files[i],dtype='<f',count=-1)
    #t_data = np.arange(shifted_t[0],shifted_t[0]+100,tbin_size)+i*100
    t_data = np.arange(0,100,tbin_size)+i*100
    print(t_data[:5],t_data[-5:])
    print(counts[:5],counts[-5:])
    plt.step(t_data,counts,'b-')


for i in range(399995,400005):
    print(t_bins[i],summed_data[i],t_bins_paul[i],sums[i])

#print(sum(summed_data[400000:800000]))
#test2 = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/ni0034070101_nicersoft_bary_GTI1_100s.dat'
#counts2 = np.fromfile(test2,dtype='<f',count=-1)
#print(counts2[:5])
#print(sum(counts2))

for i in range(9):
    actual_gti = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/100s_GTI'+str(i)+'.gti'
    actual_gti = fits.open(actual_gti)
    actual_gti = actual_gti[1].data
#    print(actual_gti[0][0],actual_gti[0][1])

testgti = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/test.gti'
testgti = fits.open(testgti)
testgti = testgti[1].data
print(testgti[0][0],testgti[0][1],type(testgti[0][0]),type(testgti[0][1]))

gti_paul = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/accelsearch_100s/ni0034070101_nicersoft_bary_GTI0_100s.evt'
gti_paul = fits.open(gti_paul)
gti_paul = gti_paul[2].data
print(gti_paul[0][0])
print(type(gti_paul[0][0]))

gti_paul = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/accelsearch_100s/ni0034070101_nicersoft_bary_GTI1_100s.evt'
gti_paul = fits.open(gti_paul)
gti_paul = gti_paul[2].data
print(gti_paul[0][0])
print(type(gti_paul[0][0]))

gti_orig = '/Volumes/Samsung_T5/NICER-data/0034070101/xti/event_cl/ni0034070101_0mpu7_cl_bary.evt'
gti_orig = fits.open(gti_orig)
gti_orig = gti_orig[2].data
print(gti_orig[0][0])
print(type(gti_orig[0][0]))
"""

"""
import Lv3_average_ps_segments

tbin_size = 0.00025
t_bins_data,counts_data = Lv3_average_ps_segments.binned_data('0034070101',['PI','TIME'],0.00025)
segments = np.arange(0,1000,100)

dat_files = sorted(Lv3_average_ps_segments.presto_dat('0034070101',100))

for i in range(len(segments)-1):
    dat_file_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
    dat_file_times = np.arange((segments[i+1]-segments[i])/tbin_size)*tbin_size+i*100
    plt.plot(t_bins_data,counts_data,'r')
    plt.plot(dat_file_times,dat_file_data,'b')
    plt.xlim([segments[i],segments[i+1]])

    plt.show()
"""

"""
To test whether presto's fft and np's fft were the same
test_dir = '/Users/masonng/Documents/MIT/Research/nicer-data/testing_prestofft/'

binary_dat = np.fromfile(test_dir+'ni0034070101_nicersoft_bary.dat',dtype='<f',count=-1)
binary_fft = np.fft.fft(binary_dat)

presto_fft = np.fromfile(test_dir+'ni0034070101_nicersoft_bary.fft',dtype='complex64',count=-1)

freqs = np.fft.fftfreq(binary_dat.size,0.00025)
N = len(freqs)

binary_ps = 2/sum(binary_dat) * np.abs(binary_fft)**2
presto_ps = 2/sum(binary_dat) * np.abs(presto_fft)**2

plt.plot(freqs[1:int(N/2)],binary_ps[1:int(N/2)],'r-')
plt.plot(freqs[1:int(N/2)],presto_ps[1:],'b-')
plt.show()
"""

"""
test_dir = '/Users/masonng/Documents/MIT/Research/nicer-data/testing_prestofft/accelsearch_100s/'
dat_files = sorted(glob.glob(test_dir+'*.dat'))
fft_files = sorted(glob.glob(test_dir+'*.fft'))

for i in range(len(dat_files)):
    binary_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
    no_photons = sum(binary_data)
    freqs = np.fft.fftfreq(binary_data.size,0.00025)
    N = len(freqs)
    binary_ps = 2/no_photons * np.abs(np.fft.fft(binary_data))**2

    fft_data = np.fromfile(fft_files[i],dtype='complex64',count=-1)
    fft_ps = 2/no_photons * np.abs(fft_data)**2

    #plt.plot(freqs[1:int(N/2)],binary_ps[1:int(N/2)],'r')
    #plt.plot(freqs[1:int(N/2)],fft_ps[1:],'b')
    #plt.show()
"""

"""
import Lv3_average_ps_segments

fft_dict = Lv3_average_ps_segments.segments_FFT('0034070101',['PI','TIME'],0.00025,100,1)
dict_keys = sorted(fft_dict.keys())

fft_files = sorted(Lv3_average_ps_segments.presto_FFT('0034070101',100))
dat_files = sorted(Lv3_average_ps_segments.presto_dat('0034070101',100))

for i in range(len(dict_keys)):
    dat_file_data = np.fromfile(dat_files[i],dtype='<f',count=-1)
    fft_data = np.fromfile(fft_files[i],dtype='complex64',count=-1)
    print('Length of fft data: ' + str(len(fft_data)) + ' , for dat data: ' + str(len(dat_file_data)))
    ps_data = 2/sum(dat_file_data) * np.abs(fft_data)**2

    freqs = np.fft.fftfreq(ps_data.size,0.00025)
    N = len(freqs)

    fft_file_freq = fft_dict[dict_keys[i]][0]
    fft_file_ps = fft_dict[dict_keys[i]][1]

    fft_presto_freq = freqs[1:int(N/2)]
    fft_presto_ps = ps_data[1:int(N/2)]

    print(len(fft_file_freq),len(fft_presto_freq),len(fft_file_ps),len(fft_presto_ps))

    plt.plot(fft_file_freq,fft_file_ps,'r')
    plt.plot(fft_presto_freq,fft_presto_ps,'b')

    plt.show()
"""

"""
P_orb = np.linspace(0,1)*86400
a_c = 0.12/(2*np.pi)*P_orb
plt.figure(1)
plt.plot(P_orb,a_c,'rx-')
plt.figure(2)
plt.plot(P_orb,a_c*3e8,'bx-')
plt.axhline(y=1.5e11,alpha=0.5,lw=0.5)
plt.show()
"""

"""
### Now test writing to new fits file...
basedir = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/accelsearch_100s/testfits/'
oldfile = basedir+'ni0034070101_nicersoft_bary_GTI0_100s.evt'
newfile = basedir+'ni0034070101_nicersoft_bary_GTI0_100s_demod.evt'
subprocess.check_call(['cp',oldfile,newfile])
with fits.open(newfile,mode='update') as hdu5:
    hdu5[1].data['TIME'] = hdu5[1].data['TIME'] + 150000
    hdu5.flush()
"""

"""
test_dir = '/Volumes/Samsung_T5/NICERsoft_outputs/trymerge/'
evt1 = test_dir + 'ni2060060363_nicersoft_bary.evt'
evt2 = test_dir + 'ni2060060364_nicersoft_bary.evt'
evt3 = test_dir + 'ni2060060365_nicersoft_bary.evt'

subprocess.check_call(['ftmerge',evt1+','+evt2+','+evt3,test_dir+'merged1.evt','copyall=YES','clobber=YES'])
"""

"""
test_dir = '/Volumes/Samsung_T5/NICERsoft_outputs/trymerge/testdirectories/'
print(sorted(glob.glob(test_dir+'merged*')))
"""

"""
test_dir = '/Volumes/Samsung_T5/NICERsoft_outputs/1200250101_pipe/'
testfile = test_dir + 'ni1200250101_nicersoft_bary_E250-1200.events'
openfile = np.fromfile(testfile,dtype='<f',count=-1)
print(len(openfile))
"""

"""
DOES NOT WORK THE WAY I WANT IT TO...
def testfunction():
    print('hi')
    if __name__ == "__main__":
        print('whatthe')

def testfunction2():
    print('lolllll')
    testfunction()

#testfunction()
testfunction2()
"""

"""
mergeddir = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000002/accelsearch_05000s/'
data_files = glob.glob(mergeddir + '*0.dat')
lc = np.zeros(20000000)
for i in tqdm(range(len(data_files))):
    dat_file = np.fromfile(data_files[i],dtype='<f',count=-1)
    lc += dat_file

lc.tofile(mergeddir+'merged_lc.dat')
"""

"""
test_parfile = '/Volumes/Samsung_T5/NICERsoft_outputs/testB1957+20.par'

T_asc = 48196.0635242 + 0.04*0.3819666069
par_contents = open(test_parfile,'r')
contents = par_contents.read()
contents = contents.split('\n')
par_contents.close()

newstring = 'TASC ' + str(T_asc) + '    #  TASC     Epoch of ascending node passage (MJD)'
par_contents = open(test_parfile,'w')
for j in range(len(contents)):
    if j != 13:
        par_contents.write(contents[j]+'\n')
    else:
        par_contents.write(newstring + '\n')
par_contents.close()
"""
"""
pyfile = 'Lv3_average_ps_main.py'
pyfile_contents = open(pyfile,'r')
contents = pyfile_contents.read().split('\n')
print(contents[116])
"""

"""
testspectrum = '/Volumes/Samsung_T5/NGC300_ULX/testdata.txt'
E,E_error,counts,counts_error = np.genfromtxt(testspectrum,skip_header=3,usecols=(0,1,2,3),unpack=True)

plt.errorbar(E,counts,xerr=E_error,yerr=counts_error,fmt='+')
plt.xlim([0.3,12])
plt.xscale('log')
plt.show()
"""

"""
par_file = '/Volumes/Samsung_T5/NICERsoft_outputs/B1957+20.par'

oldfile = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000005/merged000005_nicersoft_bary.evt' #old event FITS file
newfile = oldfile[:-4]+'_demod.evt' #new event FITS file, to be demodulated
subprocess.check_call(['cp',oldfile,newfile])
#really to prevent myself from repeating the demodulation multiple times if I run the function again...
with fits.open(newfile,mode='update') as fitsfile_demod:
    MJDREFI = fitsfile_demod[1].header['MJDREFI'] #integer for MJD reference
    MJDREFF = fitsfile_demod[1].header['MJDREFF'] #float decimal for MJD reference

    times = fitsfile_demod[1].data['TIME'] #original time series
    gtis_start = fitsfile_demod[2].data['START'] #original GTI start times
    gtis_stop = fitsfile_demod[2].data['STOP'] #original GTI end times

    times_MJD = MJDREFI + MJDREFF + times/86400 #converting METs to MJD
    gtis_start_MJD = MJDREFI + MJDREFF + gtis_start/86400 #converting GTIs in METs to MJD
    gtis_stop_MJD = MJDREFI + MJDREFF + gtis_stop/86400 #converting GTIs in METs to MJD

    try:
        times_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(times_MJD) #demodulated event times
        gtis_start_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(gtis_start_MJD) #demodulated GTI start times
        gtis_stop_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(gtis_stop_MJD) #demodulated GTI end times

        fitsfile_demod[1].data['TIME'] = (times_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs
        fitsfile_demod[2].data['START'] = (gtis_start_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs
        fitsfile_demod[2].data['STOP'] = (gtis_stop_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs

        fitsfile_demod.flush()

    except ValueError:
        pass
"""

"""
### just reproduce light curve for NGC300 ULX - what Ron sent me

cms = '/Volumes/Samsung_T5/NGC300_ULX/n300_ulx.bgsub_cl50_RGcms.ffphot'
err = '/Volumes/Samsung_T5/NGC300_ULX/n300_ulx.bgsub_cl50_RGerr_norm.ffphot'
norm = '/Volumes/Samsung_T5/NGC300_ULX/n300_ulx.bgsub_cl50_RGnorm.ffphot'

A_band,B_band,inband,MJD = np.genfromtxt(norm,usecols=(5,6,9,11),unpack=True)
A_err,B_err,inband_err,MJD_err = np.genfromtxt(err,usecols=(5,6,9,11),unpack=True)

color_BA = B_band/A_band
color_BA_err = np.sqrt( (B_err/A_band)**2 + (B_band*A_err/A_band**2)**2)

fig,(ax1,ax2) = plt.subplots(2,1)

ax1.errorbar(x=MJD,y=inband,yerr=inband_err,fmt='^',mfc='none')
ax1.axhline(y=0,ls='--',lw=0.5,alpha=0.5)
ax1.set_ylim([-2,4])

ax2.errorbar(x=MJD,y=color_BA,yerr=color_BA_err,fmt='^',mfc='none')
ax2.set_ylim([0,2])

plt.show()
"""


"""
obsids = ['12002501' + str(i).zfill(2) for i in range(1,27)]
for i in range(len(obsids)):
    event = fits.open(Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/ni' + obsids[i] + '_nicersoft_bary.evt')
    times = event[1].data['TIME']
    print(obsids[i] + ': ' + str(times[0]) + ' - ' + str(times[-1]))
"""

"""
x = np.linspace(0,1000,10000)
y = 2 * np.sin(2*np.pi*0.05*x) + np.random.normal(0,1,10000)

y2 = np.array(list(y) + list(np.zeros(10000)))
y5 = np.array(list(y) + list(np.zeros(40000)))
y50 = np.array(list(y) + list(np.zeros(490000)))

ps2 = 2/10000 * np.abs(np.fft.fft(y2))**2
freqs2 = np.fft.fftfreq(len(y2),0.1)
N2 = len(freqs2)

ps5 = 2/10000 * np.abs(np.fft.fft(y5))**2
freqs5 = np.fft.fftfreq(len(y5),0.1)
N5 = len(freqs5)

ps50 = 2/10000 * np.abs(np.fft.fft(y50))**2
freqs50 = np.fft.fftfreq(len(y50),0.1)
N50 = len(freqs50)

plt.plot(freqs2[1:int(N2/2)],ps2[1:int(N2/2)],'rx-')
plt.plot(freqs5[1:int(N5/2)],ps5[1:int(N5/2)],'bx-')
plt.plot(freqs50[1:int(N50/2)],ps50[1:int(N50/2)],'kx-')

plt.axhline(y=2)
plt.axhline(y=0.5*np.max(ps2),color='r')
plt.axhline(y=0.5*np.max(ps5),color='b')
plt.axhline(y=0.5*np.max(ps50),color='k')

plt.legend(('OS=2','OS=5','OS=50','P=2'),loc='best')
plt.show()
"""

"""
c = 299792458 #speed of light in m/s
inclin = 60 * np.pi/180 #inclination in radians
f_0 = 271.453 #pulsation frequency for J1231-1411
P_orb = 1.86 * 86400 #orbital period in seconds
G = 6.674e-11 #gravitational constant
e = 0 #eccentricity ; pretty small Laplace parameters, so just let it be 0
m_2 = 0.220 * 1.989e30 #"median" mass of companion star
m_1 = 1.4 * 1.989e30 #assumed mass of neutron star
a = 2.0426 * c/np.sin(inclin)

K_1 = np.sqrt(G/(1-e**2)) * m_2 * np.sin(inclin)/np.sqrt(a*(m_1+m_2))
f_dot = f_0 * K_1/c * 2*np.pi/P_orb
print(f_dot*500**2)
"""

"""
b1957_20_data = np.array(['ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_07/1030180101.tar.gz', 'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_07/1030180102.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_07/1030180103.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_07/1030180104.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_07/1030180105.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_07/1030180106.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_08/1030180107.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_08/1030180108.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_08/1030180109.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_08/1030180110.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_08/1030180111.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_08/1030180112.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_08/1030180113.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_10/1030180114.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_10/1030180115.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_10/1030180116.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_10/1030180117.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_10/1030180118.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_11/1030180119.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_11/1030180120.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_11/1030180121.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_11/1030180122.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_11/1030180123.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_11/1030180124.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_11/1030180125.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180126.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180127.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180128.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180129.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180130.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180131.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180132.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180133.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180134.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2017_12/1030180135.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_02/1030180136.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_02/1030180137.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_02/1030180138.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_02/1030180139.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_02/1030180140.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_02/1030180141.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_02/1030180142.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180143.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180144.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180145.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180146.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180147.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180148.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180149.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180150.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180151.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180152.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180153.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180154.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180155.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180156.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180157.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180158.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180159.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_03/1030180160.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_04/1030180161.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_04/1030180162.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_04/1030180163.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_04/1030180164.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_04/1030180165.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180166.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180167.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180168.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180169.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180170.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180171.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180172.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180173.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180174.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180175.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180176.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180177.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_07/1030180178.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180179.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180180.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180181.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180182.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180183.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180184.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180185.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180186.tar.gz',
'ftp://legacy.gsfc.nasa.gov/nicer/data/obs/2018_09/1030180187.tar.gz'])

#for i in range(len(b1957_20_data)):
#    subprocess.check_call(['wget',b1957_20_data[i]])
"""

"""
from PyAstronomy.pyasl import foldAt

gtis = Lv1_data_gtis.raw_gtis('1013010105',True)

MJDREFI = Lv0_call_eventcl.open_fits('1013010105',True)[1].header['MJDREFI']
MJDREFF = Lv0_call_eventcl.open_fits('1013010105',True)[1].header['MJDREFF']
datadict = Lv0_call_eventcl.get_eventcl('1013010105',True,['TIME','PI'])
T = sum([ (gtis[i][1]-gtis[i][0]) for i in range(len(gtis)) ])

times = datadict['TIME']
#times = MJDREFI + MJDREFF + times/86400.0 #to convert to MJD

#phase,summed_profile = Lv2_phase.pulse_folding(times,T,48442.5,29.639575,-3.77535E-10,1.1147E-20,100)
phase,summed_profile = Lv2_phase.pulse_folding(times,T,times[0],29.639575,-3.77535E-10,1.1147E-20,100)
plt.step(phase[:-1],summed_profile,'r')

plt.show()
"""

"""
basedir = '/Volumes/Samsung_T5/NGC300_ULX/'
testfits = basedir + 'test.pha'
event = fits.open(testfits)
exposure = 17949.31
chans = event[1].data['CHANNEL']
counts = event[1].data['COUNTS']
grouping = event[1].data['GROUPING']
counter = 0
for i in range(len(chans)):
    counter+=counts[i]
#    if grouping[i] == 1:
    print(chans[i],counts[i],counter,grouping[i])
"""

"""
spectra_all = '/Volumes/Samsung_T5/NGC300_ULX/get_alldata.txt'
spectra = np.array(open(spectra_all,'r').read().split('\n'))

separator = list(np.where(spectra=='NO NO NO NO')[0])
#separator.insert(0,2)
#separator.insert(len(separator),len(spectra)-1)

energy_data = [float(spectra[i].split(' ')[0]) for i in range(3,separator[0])]
energy_unc_data = [float(spectra[i].split(' ')[1]) for i in range(3,separator[0])]
flux_data = [float(spectra[i].split(' ')[2]) for i in range(3,separator[0])]
flux_unc_data = [float(spectra[i].split(' ')[3]) for i in range(3,separator[0])]

energy_bg_counts = [float(spectra[i].split(' ')[0]) for i in range(separator[0]+1,separator[1])]
energy_unc_bg_counts = [float(spectra[i].split(' ')[1]) for i in range(separator[0]+1,separator[1])]
flux_bg_counts = [float(spectra[i].split(' ')[2]) for i in range(separator[0]+1,separator[1])]
flux_unc_bg_counts = [float(spectra[i].split(' ')[3]) for i in range(separator[0]+1,separator[1])]

energy_bg_rate = [float(spectra[i].split(' ')[0]) for i in range(separator[1]+1,len(spectra)-1)]
energy_unc_bg_rate = [float(spectra[i].split(' ')[1]) for i in range(separator[1]+1,len(spectra)-1)]
flux_bg_rate = [float(spectra[i].split(' ')[2]) for i in range(separator[1]+1,len(spectra)-1)]
flux_unc_bg_rate = [float(spectra[i].split(' ')[3]) for i in range(separator[1]+1,len(spectra)-1)]

plt.errorbar(x=energy_data,y=flux_data,xerr=energy_unc_data,yerr=flux_unc_data,fmt='+-')
plt.errorbar(x=energy_bg_counts,y=flux_bg_counts,xerr=energy_unc_bg_counts,yerr=flux_unc_bg_counts,fmt='+-')
plt.errorbar(x=energy_bg_rate,y=flux_bg_rate,xerr=energy_unc_bg_rate,yerr=flux_unc_bg_rate,fmt='+-')
plt.legend(('Data','BG (counts)','BG (rate)'),loc='best')
plt.show()
"""

"""
writefile = open('/Volumes/Samsung_T5/NICERsoft_outputs/J1231_merge_filenames.txt','w')
obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113'] + [str(i) for i in range(1060060101,1060060200)] + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
for i in range(len(obsids)):
    writefile.write('/Volumes/Samsung_T5/NICER-data/' + obsids[i] + '/' + '\n')
writefile.close()
"""

"""
writefile = open('/Volumes/Samsung_T5/NICERsoft_outputs/J1231_merge_cleanfilt.txt','w')
obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113'] + [str(i) for i in range(1060060101,1060060200)] + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
for i in range(len(obsids)):
    writefile.write('/Volumes/Samsung_T5/NICERsoft_outputs/' + obsids[i] + '_pipe/ni' + obsids[i] + '_nicersoft_bary.evt' + '\n')
writefile.close()
"""

"""
import pint.toa as toa
import pint.models as models
#t = toa.get_TOAs("/Volumes/Samsung_T5/NICERsoft_outputs/1060060293_pipe/ni1060060293_nicersoft_bary.evt")
m = models.get_model("/Volumes/Samsung_T5/NICERsoft_outputs/1060060293_pipe/J1959+2048_test.par")
print(m.as_parfile())
"""

"""
#vtfile = '/Volumes/Samsung_T5/NICER-data/1013010105/xti/event_cl/ni1013010105_0mpu7_cl_bary.evt'
evtfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1060060282_pipe/accelsearch_500s/ni1060060282_nicersoft_bary_demod.evt'
#evtfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1050230107_pipe/ni1050230107_nicersoft_bary_demod.evt'
#evtfile = '/Volumes/Samsung_T5/NICERsoft_outputs/1060060293_pipe/cleanfilt.evt'
event = fits.open(evtfile)
times = event[1].data['TIME']
#phase = event[1].data['PULSE_PHASE']
#print(phase[0:100])

gtis = event[2].data
T = sum([ (gtis[i][1]-gtis[i][0]) for i in range(len(gtis)) ])
#print(T)

zerotime_crab = times[0]/86400 + 56658 + 0.000777592592592593
zerotime_j1231 = times[0]/86400 + 56658 + 0.000777592592592593
#phase_bins,profile = Lv2_phase.pulse_folding(times,T,zerotime_crab,29.639575,-3.77535E-10,1.1147E-20,50)
phase_bins,profile = Lv2_phase.pulse_folding(times,T,zerotime_j1231,271.45301962438478,-1.6670366322227105621e-15,0,50)
#phase_bins,profile = Lv2_phase.pulse_folding(times,T,58211.6,182.06580377,1.4E-12,0,50)
plt.figure(1)
plt.step(phase_bins[:-1],profile,'r')
plt.figure(2)
plt.plot(phase_bins[:-1],profile,'r-')
plt.show()
"""

"""
data = 'data'
bgsub_overplot = 'grp_bg_overplot.xcm'
base_folder = '/Volumes/Samsung_T5/NGC300_ULX/'
bgsub_files = sorted(glob.glob(base_folder + 'grp*_bg_*pha'))

data_list = []
for i in range(len(bgsub_files)):
    data_list.append(str(i+1)+':'+str(i+1))
    data_list.append(bgsub_files[i])

write_file = open(bgsub_overplot,'w')
data_string = 'data ' + ' '.join([data_list[i] for i in range(len(data_list))])
write_file.write(data_string + '\n')
write_file.write('ignore 0.0-0.285,12.01-**' + '\n')
write_file.close()
"""

"""
ffphot = "/Volumes/Samsung_T5/NGC300_ULX/n300_ulx.bgsub_cl50_RGnorm.ffphot"
exptime = np.genfromtxt(ffphot,usecols=(1),unpack=True)
print(sum(exptime))
"""

"""
startdir = '/Volumes/Samsung_T5/NICERsoft_outputs/'
total_exp = 0
for i in range(1030180149,1030180188):
    event = fits.open(startdir + str(i) + '_pipe/cleanfilt.evt')
    exposure = event[1].header['EXPOSURE']
    total_exp += exposure
    print(str(i),exposure)
"""

"""
eventslist = '/Volumes/Samsung_T5/NICERsoft_outputs/B1957_events_ELL1.txt'
#orbslist = '/Volumes/Samsung_T5/NICERsoft_outputs/B1957_orbs.txt'

events = open(eventslist,'w')
#orbs = open(orbslist,'w')
for i in range(1030180101,1030180188):
    events.write('/Volumes/Samsung_T5/NICERsoft_outputs/' + str(i) + '_pipe/cleanfilt_phase_ELL1.evt' + '\n')
#    orbs.write('/Volumes/Samsung_T5/NICERsoft_outputs/' + str(i) + '_pipe/ni' + str(i) + '.orb' + '\n')
events.close()
#orbs.close()
"""


"""
import pint.models as models

#parfile = '/Volumes/Samsung_T5/NICERsoft_outputs/J1959+2048_NUPPI_forPaul.par'
parfile = '/Volumes/Samsung_T5/NICERsoft_outputs/J1231-1411_paul.par'
m = models.get_model(parfile)
print(m.as_parfile())
"""

"""
jsgrp_files = glob.glob('/Volumes/Samsung_T5/NGC300_ULX/jsgrp*cl*pha')
for i in range(len(jsgrp_files)):
    event = fits.open(jsgrp_files[i])
    print(jsgrp_files[i],event[1].header['BACKFILE'])
"""

"""
#mjds = ['58239','58244','58249','58254','58259','58264','58269','58274',
#        '58289','58309','58314','58324','58329','58334','58339','58389',
#        '58399','58449','58454','58459','58464','58484','58489','58504',
#        '58509'] #for 0.4-5 keV
mjds = ['58239','58244','58249','58254','58259','58264','58269','58274',
        '58289','58309','58314','58324','58329','58389','58449']

jsgrp_files = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/jsgrp_58*_cl_*pha'))
jsgrp_files_spectra = [jsgrp_files[i] for i in range(len(jsgrp_files)) if jsgrp_files[i][-21:-16] in mjds]
xcm_file = '/Volumes/Samsung_T5/NGC300_ULX/jsgrp_cl_spectralfit_0050_0200.xcm'
writing = open(xcm_file,'w')
command1 = 'data '
for i in range(len(mjds)):
    writing.write(command1 + str(i+1)+':'+str(i+1) + ' ' + jsgrp_files_spectra[i] + ' ')
writing.write('\n')
writing.write('ignore **:0.0-0.485,2.01-**')
writing.close()
"""

"""
ell1_ex = '/Volumes/Samsung_T5/NICERsoft_outputs/1030180121_pipe/cleanfilt_phase_ELL1.evt'
bt_ex = '/Volumes/Samsung_T5/NICERsoft_outputs/1030180121_pipe/cleanfilt_phase_BT.evt'
ell1_phase = fits.open(ell1_ex)[1].data['pulse_phase']
bt_phase = fits.open(bt_ex)[1].data['pulse_phase']
#for i in range(len(ell1_phase)):
#    print(np.abs(ell1_phase[i]-bt_phase[i])/ell1_phase[i]*100)
print(len(ell1_phase),len(bt_phase))
"""

"""
from pint.eventstats import sf_z2m,z2m,sig2sigma
#merged_event = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000009/merged000009_cut.evt' #assuming I do a 2.5c/s rate cut
merged_event = '/Volumes/Samsung_T5/NICERsoft_outputs/B1957_merged_BT/B1957_merged_BT_cut.evt'
#10/24 for B1957 - lmao no. 40 harmonics needed.... ya...
phases = fits.open(merged_event)[1].data['PULSE_PHASE']

m = 100
z_vals = z2m(phases,m=m)
probs = sf_z2m(z_vals)
significances = sig2sigma(probs)
print(significances)
"""

"""
parfile = '/Users/masonng/Downloads/test.par'
contents = open(parfile,'r').read().split('\n')
for i in range(len(contents)):
    print('A1' in contents[i])
"""

"""
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)

F0_txt = '/Volumes/Samsung_T5/NICERsoft_outputs/J1231-1411_F0.txt'
freq,harm1,harm2 = np.genfromtxt(F0_txt,dtype='str',skip_footer=1,usecols=(1,4,5),unpack=True)
freq = [float(freq[i][:-1]) for i in range(len(freq))]
harm1 = [float(harm1[i][1:]) for i in range(len(harm1))]
harm2 = [float(harm2[i][:-1]) for i in range(len(harm2))]
plt.plot(freq,harm1,'rx-')
plt.plot(freq,harm2,'bx-')
plt.xlabel('Frequency',fontsize=12)
plt.ylabel('Significance (sigma)',fontsize=12)
plt.legend(('1st Harmonic','2nd Harmonic'),loc='best')
plt.show()
"""

"""
from pint.eventstats import sigma2sig

print(sigma2sig(14.5))
print(sigma2sig(8.25))
"""

"""
xmm_newton = '/Volumes/Samsung_T5/NICERsoft_outputs/B1957+20_XMM/P0204910201PNS003TIEVLI0000.FTZ'
mike = '/Volumes/Samsung_T5/NICERsoft_outputs/B1957+20_XMM/pn.evt9.psrj1959+2048.bary.phase.fits'

event = fits.open(mike)
pis = event[1].data['PI']
phases = event[1].data['PULSE_PHASE']
filtered = pis[(pis>=500)&(pis<=4500)]
phases_filtered = phases[(pis>=500)&(pis<=4500)]

profbins = np.linspace(0.0,1.0,25+1,endpoint=True)
plt.hist(phases_filtered,bins=profbins)
plt.show()
print(len(filtered))
print(len(pis))
"""

"""
test_lumin = '/Volumes/Samsung_T5/NGC300_ULX/spectral_fit_0040-0500/tbnew-powerlaw_lumin.txt'
contents = np.array(open(test_lumin,'r').read().split('\n'))
lumin_lines = [float(contents[i].split(' ')[2]) for i in range(2,len(contents)-1,3)] #lines containing luminosity
print(lumin_lines)
print(len(lumin_lines))
"""

"""
mjds_used = np.array([58239,58244,58249,58254,58259,58264,58269,58274,58289,
             58309,58314,58324,58329,58334,58339,58389,58399,58449,
             58454,58459,58464,58484,58489,58504,58509])

#21
mjds_used = np.array([58239,58244,58249,58254,58259,58264,58269,58274,
            58309,58314,58324,58329,58334,58339,58389,58449,
            58454,58484,58489,58504,58509])

update_hid = '/Volumes/Samsung_T5/NGC300_ULX/soft_color_HID.txt'
mjd,soft,soft_err,intens,intens_err = np.genfromtxt(update_hid,skip_header=1,unpack=True,usecols=(0,1,2,3,4))

new_mjd = [mjd[i] for i in range(len(mjd)) if mjd[i] in mjds_used]
new_soft = [soft[i] for i in range(len(mjd)) if mjd[i] in mjds_used]
new_soft_err = [soft_err[i] for i in range(len(mjd)) if mjd[i] in mjds_used]
new_intens = [intens[i] for i in range(len(mjd)) if mjd[i] in mjds_used]
new_intens_err = [intens_err[i] for i in range(len(mjd)) if mjd[i] in mjds_used]

#plt.errorbar(x=soft,y=intens,xerr=soft_err,yerr=intens_err,fmt='rx')
plt.errorbar(x=new_soft,y=new_intens,xerr=new_soft_err,yerr=new_intens_err,fmt='kx')
plt.errorbar(x=new_soft[0],y=new_intens[0],xerr=new_soft_err[0],yerr=new_intens_err[0],fmt='rx',markersize=10)
plt.errorbar(x=new_soft[-1],y=new_intens[-1],xerr=new_soft_err[-1],yerr=new_intens_err[-1],fmt='bx',markersize=10)
plt.legend(('Data','Sample 1','Sample 2'),loc='best',fontsize=12)
#plt.xlim([-0.1,2.9])
#plt.ylim([-0.9,1.6])
plt.xlabel('Soft Color: (1-2 keV)/(0.4-1 keV)',fontsize=12)
plt.ylabel('Intensity in counts/s (0.4-12 keV)',fontsize=12)
plt.show()
"""


"""
uftxt = '/Volumes/Samsung_T5/NGC300_ULX/spectral_fit_0040-0500_25_original/tbnew-powerlaw_ufspectra.txt'
print(np.genfromtxt(uftxt,skip_header=3,skip_footer=1,usecols=(tuple(range(4+1))),unpack=True))
"""

"""
cleanfile = '/Volumes/Samsung_T5/NICERsoft_outputs/B1957+20_XMM/SAS_DataAnalysisThread/P0204910201PNS003TIEVLI0000_clean.FTZ'
pis = fits.open(cleanfile)[1].data['PI']
filtered_pis = pis[(pis>=500)&(pis<=4500)]
print(len(filtered_pis))
"""

"""
exposure = []
phafiles = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/jsgrp_58*pha'))
for i in range(len(phafiles)):
    event = fits.open(phafiles[i])
    exptime = event[1].header['EXPOSURE']
    exposure.append(exptime)
    print(phafiles[i] + ': ' + str(exptime))

used = np.array([17949.31,39294.18,2934.0,8829.248,6688.78,6437.754,6156.472,
                6765.789,6232.0,8644.451,2744.0,4581.0,7116.36,6893.0,5685.0,
                7315.0,6258.0,5094.0,7354.0,10082.0,7088.0])
print(len(used),sum(used)/21)
"""

"""
import Lv3_analyze_xspec_pars
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

powerlaw_vals = Lv3_analyze_xspec_pars.xspec_par('tbnew-powerlaw','0040','0500')
MJDs_21 = ['58239','58244','58249','58254','58259','58264','58269','58274',
            '58309','58314','58324','58329','58334','58339','58389','58449',
            '58454','58484','58489','58504','58509']

x_data = np.array(MJDs_21,dtype='float')
y_data = powerlaw_vals['powerlaw-PhoIndex']
y_data_err = powerlaw_vals['powerlaw-PhoIndex_unc']

plt.figure(figsize=(16,9))
plt.errorbar(x=x_data,y=y_data,yerr=y_data_err,fmt='kx-')
plt.errorbar(x=x_data[0],y=y_data[0],yerr=y_data_err[0],fmt='rx-',markersize=10)
plt.errorbar(x=x_data[-1],y=y_data[-1],yerr=y_data_err[-1],fmt='bx-',markersize=10)
plt.legend(('XSPEC Fit','Sample 1','Sample 2'),fontsize=12,loc='best')
plt.xlabel('Time (MJD)',fontsize=12)
plt.ylabel(r'powerlaw - $\Gamma$ ',fontsize=12)
plt.show()
"""


"""
MJDs_21 = ['58239','58244','58249','58254','58259','58264','58269','58274',
            '58309','58314','58324','58329','58334','58339','58389','58449',
            '58454','58484','58489','58504','58509']

powerlaw = sorted(glob.glob(Lv0_dirs.NGC300+'spectral_fit_0040-0500/indiv_ratio/tbnew-powerlaw_*.txt'))
powerlaw_gauss = sorted(glob.glob(Lv0_dirs.NGC300+'spectral_fit_0040-0500/indiv_ratio/tbnew-powerlaw-gauss_*.txt'))
powerlaw_laor = sorted(glob.glob(Lv0_dirs.NGC300+'spectral_fit_0040-0500/indiv_ratio/tbnew-powerlaw-laor_*.txt'))
for i in range(len(powerlaw)):
    pl_e,pl_e_err,pl_r,pl_r_err = np.genfromtxt(powerlaw[i],usecols=(0,1,2,3),unpack=True)
    pl_gauss_e,pl_gauss_e_err,pl_gauss_r,pl_gauss_r_err = np.genfromtxt(powerlaw_gauss[i],usecols=(0,1,2,3),unpack=True)
    pl_laor_e,pl_laor_e_err,pl_laor_r,pl_laor_r_err = np.genfromtxt(powerlaw_laor[i],usecols=(0,1,2,3),unpack=True)

    plt.figure(figsize=(16,9))
    plt.errorbar(x=pl_e,y=pl_r,xerr=pl_e_err,yerr=pl_r_err,fmt='+',alpha=0.5)
    plt.errorbar(x=pl_gauss_e,y=pl_gauss_r,xerr=pl_gauss_e_err,yerr=pl_gauss_r_err,fmt='+',alpha=0.5)
    #plt.errorbar(x=pl_laor_e,y=pl_laor_r,xerr=pl_laor_e_err,yerr=pl_laor_r_err,fmt='+',alpha=0.5)
    plt.legend(('powerlaw','powerlaw-gauss','powerlaw-laor'),fontsize=12)
    plt.axhline(y=1,lw=0.5,alpha=0.5)
    plt.xlabel('Energy, E (keV)',fontsize=12)
    plt.ylabel('ratio (data/model)',fontsize=12)
    plt.title('MJD ' + MJDs_21[i],fontsize=12)

    plt.show()
"""

"""
powerlaw = sorted(glob.glob(Lv0_dirs.NGC300+'spectral_fit_0040-0500/indiv_ratio/tbnew-powerlaw_*.txt'))
powerlaw_gauss = sorted(glob.glob(Lv0_dirs.NGC300+'spectral_fit_0040-0500/indiv_ratio/tbnew-powerlaw-gauss_*.txt'))
powerlaw_laor = sorted(glob.glob(Lv0_dirs.NGC300+'spectral_fit_0040-0500/indiv_ratio/tbnew-powerlaw-laor_*.txt'))
for i in range(len(powerlaw)):
    pl_e,pl_e_err,pl_chi= np.genfromtxt(powerlaw[i],usecols=(0,1,2),unpack=True)
    pl_gauss_e,pl_gauss_e_err,pl_gauss_chi = np.genfromtxt(powerlaw_gauss[i],usecols=(0,1,2),unpack=True)
    pl_laor_e,pl_laor_e_err,pl_laor_chi = np.genfromtxt(powerlaw_laor[i],usecols=(0,1,2),unpack=True)

    plt.figure(figsize=(16,9))
    plt.axhline(y=1,lw=0.5,alpha=0.5)
    plt.step(x=pl_e,y=pl_chi)
    plt.step(x=pl_gauss_e,y=pl_gauss_chi)
    plt.step(x=pl_laor_e,y=pl_laor_chi)
    plt.legend(('y=0','powerlaw','powerlaw-gauss','powerlaw_laor'),fontsize=12)
    plt.xlabel('Energy, E (keV)',fontsize=12)
    plt.ylabel('ratio (data/model)',fontsize=12)
    plt.title(MJDs_21[i],fontsize=12)

    plt.show()
"""

"""
mjd58239 = '/Volumes/Samsung_T5/NGC300_ULX/spectral_fit_0040-0500/indiv_spectra/tbnew-powerlaw_ufspectra_58239.txt'
mjd58504 = '/Volumes/Samsung_T5/NGC300_ULX/spectral_fit_0040-0500/indiv_spectra/tbnew-powerlaw_ufspectra_58504.txt'

day1_E,day1_E_err,day1_flux,day1_flux_err,day1_model = np.genfromtxt(mjd58239,usecols=(0,1,2,3,4),unpack=True)
day2_E,day2_E_err,day2_flux,day2_flux_err,day2_model = np.genfromtxt(mjd58504,usecols=(0,1,2,3,4),unpack=True)

plt.figure(figsize=(16,9))
plt.yscale('log')
plt.xscale('log')

plt.errorbar(x=day1_E,y=day1_flux,xerr=day1_E_err,yerr=day1_flux_err,fmt='+',color='r')
plt.errorbar(x=day2_E,y=day2_flux,xerr=day2_E_err,yerr=day2_flux_err,fmt='+',color='b')

plt.errorbar(x=day1_E,y=day1_model,xerr=day1_E_err,fmt='-',color='k')
plt.errorbar(x=day2_E,y=day2_model,xerr=day2_E_err,fmt='-',color='k')

plt.legend(('Sample 1','Sample 2','powerlaw'),loc='best',fontsize=12)
plt.xlabel('Energy, E (keV)',fontsize=12)
plt.ylabel('Flux (photons/cm^2/s/keV)',fontsize=12)
plt.xlim([0.4,5])

#plt.annotate('For 21 spectra, chi-squared is 3910/2729 (1.43)',(0.5,1e-5))

plt.show()
"""


timeend = time.time()

print(str(timeend-timestart) + ' seconds')
