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
from presto import binary_psr
import matplotlib.gridspec as gridspec
import subprocess
from astropy.wcs import WCS,utils
#import pint.toa as toa
import Lv0_fits2dict,Lv2_dj_lsp
import glob
from PyAstronomy.pyasl import foldAt
import peakutils
import pathlib
import Lv0_dirs,Lv1_data_bin,Lv2_phase,Lv3_detection_level
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy import signal
from scipy.optimize import curve_fit
#from astropy.visualization import quantity_support
#from stingray.pulse.pulsar import pulse_phase,phase_exposure,fold_events
#from stingray.pulse.search import plot_profile

#quantity_support()

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

"""
def GaussSum(x,*p):
    n=int(len(p)/3)
    A=p[:n] #amplitude
    w=p[n:2*n] #std
    c=p[2*n:3*n] #centroid
    const = p[-1]
    return sum([ A[i]*np.exp(-(x-c[i])**2/(2*w[i]**2))/np.sqrt(2*np.pi*w[i]**2) for i in range(n)]) + const

x = np.linspace(1,1000,10001)
y = GaussSum(x,[5,10,3,5,100,200,50])
print(len(x),len(y))
plt.plot(x,y,'rx-')
plt.show()


efsearch_results = '/Volumes/Samsung_T5/NICER-data/rxj0209/rxj0209_100segs.fes'
seg1_period = fits.open(efsearch_results)[6].data['PERIOD']
seg1_chisq = fits.open(efsearch_results)[6].data['CHISQRD1']
seg1_error = fits.open(efsearch_results)[6].data['ERROR1']
mindist = 0.0001

peak_indices = peakutils.indexes(seg1_chisq,thres=0,min_dist=mindist)
N_gauss = len(peak_indices)

amp_guess = [seg1_chisq[peak_indices[i]] - np.min(seg1_chisq) for i in range(N_gauss)]
std_guess = [mindist*5] * N_gauss
mean_guess = [seg1_period[peak_indices[i]] for i in range(N_gauss)]
pguess = amp_guess + std_guess + mean_guess + [np.mean(seg1_chisq)]

popt,pcov = curve_fit(f=GaussSum,xdata=seg1_period,ydata=seg1_chisq,sigma=seg1_error,p0=pguess)
for i in range(len(popt)):
    print(popt[i],np.sqrt(np.diag(pcov))[i])
plt.plot(seg1_period,GaussSum(seg1_period,*popt),'b-')

plt.errorbar(x=seg1_period,y=seg1_chisq,yerr=seg1_error)
for i in range(len(peak_indices)):
    plt.plot(seg1_period[peak_indices[i]],seg1_chisq[peak_indices[i]],'rx')
    plt.text(seg1_period[peak_indices[i]],1.1*np.max(seg1_chisq),str(i) + '|',fontsize=8)
plt.show()
"""

"""
print(Lv3_detection_level.signal_significance(10/0.0005,1,1,31.633))
print(Lv3_detection_level.power_for_sigma(3,10/0.0005,1,1))
"""

"""
for i in range(1200250101,1200250127): #for all of the AT2018cow ObsIDs:
    gtis = fits.open(Lv0_dirs.NICERSOFT_DATADIR + str(i) + '_pipe/ni' + str(i) + '_nicersoft_bary.evt')[2].data
    print('ObsID: ' + str(i) + ' ; Start time in GTI: ' + str(gtis[0][0]) + ' ; End time in GTI: ' + str(gtis[-1][1]))
"""

"""
r_sun = 6.96e8 #radius of sun in m
P = 88*86400 #orbital period of mercury
a = 0.387 * 1.5e11 #orbital radius in m
e = 0.206 #eccentricity of mercury orbit
M_sun = 1.989e30 #mass of sun
G = 6.674e-11 #grav. constant

omega_sun = 2*np.pi/(25*86400)
omega_b = np.sqrt(G*M_sun/r_sun**3)

omega_prec = 6*np.pi*r_sun**2/(P*(a*(1-e**2))**2) * (omega_sun/omega_b)
print(omega_prec)
print(omega_prec*180/np.pi*3600 * 86400 * 365 * 100)
"""

"""
gtifolder = '/Volumes/Samsung_T5/NICER-data/xtej1739/accelsearch_GTIs/GTI'
for i in range(1,93):
    gti = gtifolder+str(i)+'.gti'
    gti_fits = fits.open(gti)[1].data
    print(i,gti_fits[0][0]-193073992.57127872,gti_fits[0][1]-193073992.57127872)
"""

"""
eventfile = '/Users/masonng/Documents/MIT/Research/xtej1739_burst8/2002251704_bary.evt'
gtifolder = '/Users/masonng/Documents/MIT/Research/xtej1739_burst8/accelsearch_GTIs/'
#gtistart = fits.open('/Users/masonng/Documents/MIT/Research/xtej1739_burst8/2002251704_bary.evt')[2].data[0][0]

#for i,starttime in enumerate(np.arange(gtistart+100,gtistart+205,5)):
    #subprocess.run(['mkgti.py','--gtiname',gtifolder+'GTI'+str(i+1)+'.gti',str(starttime),str(starttime+10)])
#    subprocess.run(['niextract-events',eventfile,gtifolder+'2002251704_bary_GTI'+str(i+1)+'.evt','timefile='+gtifolder+'GTI'+str(i+1)+'.gti'])

def exp_func(x,a,b,c,d):
    return a*np.exp(-b*x+c) + d #a for amplitude, b for decay constant, and c for a constant background

times = fits.open(eventfile)[1].data['TIME']
times_t = times-times[0]
counts = np.ones(len(times_t))

t_bins = np.linspace(0,np.ceil(times_t[-1]),np.ceil(times_t[-1])+1)
summed_data, bin_edges, binnumber = stats.binned_statistic(times_t,counts,statistic='sum',bins=t_bins)

plt.plot(t_bins[:-1],summed_data,'r-')

#popt,pcov = curve_fit(exp_func,t_bins[:-1],summed_data,p0=[350,10,50])
#print(popt,pcov)

plt.plot(t_bins[:-1],exp_func(t_bins[:-1],20,0.02,2,50),'b-')

for i in range(1,22):
    gtis = fits.open(gtifolder+'GTI'+str(i)+'.gti')[1].data
    start = gtis[0][0]-times[0]
    end = gtis[0][1]-times[0]
    plt.axvline(x=start)
    plt.axvline(x=end)

plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.axvline(x=150,color='k')
plt.axvline(x=160,color='k')
"""

"""
#plt.show()
timezero = 193073992.641447
artifacts = [617550,618550,645580,646500,673350,674375,741020,741290,974690,975200,985830,986350,1013740,1014200,1041630,1041990,
            1489670,1490520,1529070,1529860,1589670,1590900,1595745,1596900,1612070,1612660,1617800,1619200,1656460,1658240,
            1718650,1719220,1774380,1775140,1841220,1842320,1880190,1881360,1919150,1920400,1969270,1970120,2047220,2048070,
            2052800,2053640,2080650,2081500,2152530,2153700,2181020,2181850,2253630,2255040]

"""

"""
eventfile = Lv0_dirs.NICER_DATADIR + 'xtej1739_burst27+/2002131540_bary.evt'
gtis = fits.open(eventfile)[2].data
times = fits.open(eventfile)[1].data['TIME']
segment_length = 10000

Tobs_start = gtis[0][0] #MJD for the observation start time
Tobs_end = gtis[-1][1] #MJD for the observation end time

segment_times = np.arange(Tobs_start,Tobs_end+segment_length,segment_length) #array of time values, starting
binned_counts, bin_edges, binnumber = stats.binned_statistic(times,np.ones(len(times)),statistic='sum',bins=segment_times)
"""

#for i in range(len(binned_counts)):
#    print(binned_counts[i],bin_edges[i])

"""
gtis = list(fits.open(eventfile)[2].data)
rm_artifacts = 193073992.641447008+np.array(open( str(pathlib.Path(eventfile).parent)+'/rm_artifacts.txt','r').read().split('\n')[:-1],dtype=np.float64)

N = len(gtis)

gtis_remove = []
for i in range(len(gtis)):
    for j in range(0,len(rm_artifacts)-1,2):
        if (gtis[i][0] >= rm_artifacts[j]) and (gtis[i][1] <= rm_artifacts[j+1]):
            gtis_remove.append(gtis[i])
        #    break

new_gtis = []
for i in range(len(gtis)):
    if gtis[i] not in np.array(gtis_remove):
        new_gtis.append(gtis[i])
print(len(new_gtis))
"""

#eventfile = Lv0_dirs.NICER_DATADIR + 'xtej1739_burst27+/2002131540.evt'
#t = toa.get_TOAs(eventfile, usepickle=True)

"""
eventfile = '/Volumes/Samsung_T5/NICER-data/1034100102/xti/event_cl/ni1034100102_0mpu7_cl_bary.evt'
import pint.event_toas as event_toas
t = event_toas.load_event_TOAs(eventfile,'nicer')
print(t[0],t[1])
tt = t.table
tt.show_in_browser()
"""

"""
txtfile = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000017/S500_W5_T30_E30-200.txt'
f,ps = np.genfromtxt(txtfile,usecols=(0,1),unpack=True)
print(f)
"""

"""
for i in range(101,128):
    print('1034100'+str(i))

merged_files = [Lv0_dirs.NICERSOFT_DATADIR+'merged_events/merged0000'+str(i)+'/merged0000'+str(i)+'_nicersoft_bary.evt' for i in range(15,19)]
for i in range(len(merged_files)):
    gtis = fits.open(merged_files[i])[2].data
    print(len(gtis))
    differences = np.array( [ (gtis[j][1]-gtis[j][0]) for j in range(len(gtis)) ])
    print(sum(differences))
"""

"""
filenames = sorted(glob.glob('/Volumes/Samsung_T5/NICER-data/xtej1812/accelsearch_64s/*.evt'))
for i in range(len(filenames)):
    print(filenames[i],len(fits.open(filenames[i])[1].data))
"""

"""
tbin_size = 16

eventfile = Lv0_dirs.NICER_DATADIR + 'xtej1701/2003261259_filt_bary.evt'
times = fits.open(eventfile)[1].data['TIME']
trunc = times - times[0]
t_bins = np.linspace(0,np.ceil(trunc[-1]),np.ceil(trunc[-1])*1/tbin_size+1)
summed_data, bin_edges, binnumber = stats.binned_statistic(trunc,np.ones(len(times)),statistic='sum',bins=t_bins)

spectra_GTIs = np.array([1,2,3,4,5,6,7,8,9,10,11,12,20]) #corresponding GTI no. that had spectra extracted
gti_start,gti_stop = np.loadtxt('/Volumes/Samsung_T5/NICER-data/xtej1701/bunchedgtis.txt',usecols=(0,1),unpack=True)

gti_used_start = np.array([gti_start[i-1] for i in spectra_GTIs]) - times[0]
gti_used_stop = np.array([gti_stop[i-1] for i in spectra_GTIs]) - times[0]
for i in range(len(gti_used_start)):
    print(gti_used_stop[i]-gti_used_start[i])
gti_used_centroid = gti_used_start + (gti_used_stop-gti_used_start)/2
"""

"""
fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

time_bins = t_bins[:-1]
count_bins = summed_data/tbin_size
ax1.plot(time_bins[count_bins>0],count_bins[count_bins>0],'rx') #plot the light curve
for i in range(len(gti_used_start)):
    ax1.axvline(x=gti_used_start[i],alpha=0.5,lw=0.5)
    ax1.axvline(x=gti_used_stop[i],alpha=0.5,lw=0.5)

ax1.set_ylabel('Counts/s',fontsize=12)

#nH fixed
gamma_fixed = np.array([1.6687,1.5988,1.5699,1.7324,1.7758,1.8787,1.8873,1.8097,1.8023,2.1124,2.0908,2.2772,2.1631])
gamma_fixed_err = np.array([0.011,0.010,0.011,0.011,0.011,0.014,0.012,0.024,0.023,0.041,0.015,0.045,0.022])
norm_fixed = np.array([0.168,0.168,0.152,0.164,0.161,0.159,0.148,0.153,0.143,0.091,0.102,0.092,0.083])
norm_fixed_err = np.array([0.0025,0.0023,0.0021,0.0024,0.0024,0.0028,0.0023,0.0045,0.0040,0.0042,0.0019,0.0046,0.0021])

ax2.errorbar(x=gti_used_centroid,y=gamma_fixed,yerr=gamma_fixed_err)
ax2.set_ylabel('Gamma',fontsize=12)

ax3.errorbar(x=gti_used_centroid,y=norm_fixed,yerr=norm_fixed_err)
ax3.set_ylabel('norm',fontsize=12)
ax3.set_xlabel('Time from first event (s)',fontsize=12)

plt.show()
"""

"""
fig,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex=True)

time_bins = t_bins[:-1]
count_bins = summed_data/tbin_size
ax1.plot(time_bins[count_bins>0],count_bins[count_bins>0],'rx') #plot the light curve
for i in range(len(gti_used_start)):
    ax1.axvline(x=gti_used_start[i],alpha=0.5,lw=0.5)
    ax1.axvline(x=gti_used_stop[i],alpha=0.5,lw=0.5)

ax1.set_ylabel('Counts/s',fontsize=12)

#nH not fixed
nH_not_fixed = np.array([4.501,4.574,4.311,4.735,4.700,4.695,4.680,4.981,4.940,5.534,4.508,4.777,4.797])
nH_not_fixed_err = np.array([0.048,0.044,0.044,0.048,0.048,0.062,0.051,0.122,0.111,0.215,0.063,0.199,0.098])
gamma_not_fixed = np.array([1.633,1.586,1.478,1.769,1.803,1.904,1.908,1.933,1.912,2.481,2.053,2.343,2.230])
gamma_not_fixed_err = np.array([0.018,0.016,0.016,0.018,0.019,0.025,0.021,0.048,0.044,0.094,0.027,0.094,0.043])
norm_not_fixed = np.array([0.158,0.164,0.130,0.174,0.168,0.166,0.153,0.187,0.170,0.159,0.096,0.102,0.092])
norm_not_fixed_err = np.array([0.004,0.004,0.003,0.005,0.005,0.006,0.005,0.013,0.011,0.021,0.004,0.013,0.006])

ax2.errorbar(x=gti_used_centroid,y=nH_not_fixed,yerr=nH_not_fixed_err)
ax2.set_ylabel('nH',fontsize=12)

ax3.errorbar(x=gti_used_centroid,y=gamma_not_fixed,yerr=gamma_not_fixed_err)
ax3.set_ylabel('Gamma',fontsize=12)

ax4.errorbar(x=gti_used_centroid,y=norm_not_fixed,yerr=norm_not_fixed_err)
ax4.set_ylabel('norm',fontsize=12)

ax4.set_xlabel('Time from first event (s)',fontsize=12)

plt.show()
"""

"""
tbin_size = 16
eventfile = Lv0_dirs.NICER_DATADIR + 'xtej1739_mostrec/2002131540_filt_bary.evt'
times = fits.open(eventfile)[1].data['TIME']
trunc = times - times[0]
t_bins = np.linspace(0,np.ceil(trunc[-1]),np.ceil(trunc[-1])*1/tbin_size+1)
summed_data, bin_edges, binnumber = stats.binned_statistic(trunc,np.ones(len(times)),statistic='sum',bins=t_bins)

spectra_GTIs = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,36,37,38,39,
                        46,60,61,62,63,64,65,66,67,77,78,83,84,85,86,91,92,93,94,97,98,
                        103,106,107,108,110,117,120,121,122,123,125,126,127,135,136,
                        139,146,150,152,155,156,157,158,161,162,167,168,169,170,171,
                        172,179,188,197,198,199,200,201,203,204,207,208,213,214,238,
                        239,240,241,242,244]) #corresponding GTI no. that had spectra extracted
gti_start,gti_stop = np.loadtxt('/Volumes/Samsung_T5/NICER-data/xtej1739_mostrec/bunchedgtis.txt',usecols=(0,1),unpack=True)

gti_used_start = np.array([gti_start[i-1] for i in spectra_GTIs]) - times[0]
gti_used_stop = np.array([gti_stop[i-1] for i in spectra_GTIs]) - times[0]
#for i in range(len(gti_used_start)):
#    print(gti_used_stop[i]-gti_used_start[i])
gti_used_centroid = gti_used_start + (gti_used_stop-gti_used_start)/2

contents = open('/Volumes/Samsung_T5/NICER-data/xtej1739_mostrec/accelsearch_GTIs/tbabs-powerlaw.txt','r').read().split('\n')
phoindex = []
phoindex_err = []
norm = []
norm_err = []

for i in range(len(contents)): #for each line
    if "PhoIndex" in contents[i]:
        line = [x for x in contents[i].split(' ') if x]
        phoindex.append(float(line[4]))
        phoindex_err.append(float(line[6]))
    elif "norm" in contents[i]:
        line = [x for x in contents[i].split(' ') if x]
        norm.append(float(line[4]))
        norm_err.append(float(line[6]))

fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

time_bins = t_bins[:-1]
count_bins = summed_data/tbin_size
ax1.plot(time_bins[count_bins>0],count_bins[count_bins>0],'rx') #plot the light curve
for i in range(len(gti_used_start)):
    ax1.axvline(x=gti_used_start[i],alpha=0.5,lw=0.5)
    ax1.axvline(x=gti_used_stop[i],alpha=0.5,lw=0.5)

ax1.set_ylabel('Counts/s',fontsize=12)

ax2.errorbar(x=gti_used_centroid,y=phoindex,yerr=phoindex_err,fmt='x-')
ax2.set_ylabel('Gamma',fontsize=12)

ax3.errorbar(x=gti_used_centroid,y=norm,yerr=norm_err,fmt='x-')
ax3.set_ylabel('norm',fontsize=12)
ax3.set_xlabel('Time from first event (s)',fontsize=12)

plt.show()
"""

"""
counts = 0
gtis_no = 0
for i in range(2050280101,2050280125):
    fitsfile = Lv0_dirs.NICERSOFT_DATADIR + str(i)+'_pipe/ni'+str(i)+'_nicersoft_Bary.evt'
    counts += len(fits.open(fitsfile)[1].data['TIME'])
    gtis_no += len(fits.open(fitsfile)[2].data['START'])

print(counts)
print(gtis_no)
"""

"""
evtfiles = sorted(glob.glob(Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/merged000019/accelsearch_500s/merged000019_nicersoft_bary_GTI8**_500s_E200-1200.evt'))
for i in range(len(evtfiles)):
    times = fits.open(evtfiles[i])[1].data['TIME']
    trunc = times-times[0]
    t_bins = np.linspace(0,np.ceil(trunc[-1]),np.ceil(trunc[-1])*1/1+1)
    summed_data, bin_edges, binnumber = stats.binned_statistic(trunc,np.ones(len(times)),statistic='sum',bins=t_bins)

    plt.plot(t_bins[:-1],summed_data,'rx-')
    plt.show()
"""

"""
ngc300x1_files = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX_Swift/ngc300ulx-1/xrt/event/*pc*po*cl*evt'))
exp = 0
for i in range(len(ngc300x1_files)):
    exptime = fits.open(ngc300x1_files[i])[1].header['EXPOSURE']
    exp += exptime
print(exp)
"""

"""
import Lv2_ps_method
eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge.evt'
times = fits.open(eventfile)[1].data['TIME']
truncated_times = (times-times[0])/3600
tbins = np.arange(0,np.ceil(truncated_times[-1]))

counts,bin_edges,binnumber=stats.binned_statistic(truncated_times,np.ones(len(truncated_times)),statistic='sum',bins=tbins)

nout = 1e4
f = np.linspace(0.01,1,nout)
import scipy.signal as signal
pgram = signal.lombscargle(times,np.ones(len(times)),f,normalize=True)
plt.plot(f,pgram,'r-')

#f,ps = Lv2_ps_method.manual(tbins[:-1],counts,[False],[False],True,[False,0,0])
#plt.plot(f,ps,'r-')

#plt.plot(tbins[:-1],counts,'rx')
#plt.xlabel('Time (h)',fontsize=12)
#plt.ylabel('Counts/hr', fontsize=12)
plt.show()
"""

"""
t = np.linspace(0,100,10001)
y = 0.5*np.sin(3*t) + np.random.normal(0,0.01,10001)
f = np.linspace(0.01,5,10001)

indices = np.random.randint(0,10001,100)
filt_t = t#[indices]
filt_y = y#[indices]

#filt_t = t[::100]
#filt_y = y[::100]

ffilt_t = []
ffilt_y = []

bound1 = 10
bound2 = 20
bound3 = 40
bound4 = 50
bound5 = 70
bound6 = 90

for i in range(len(filt_t)):
    if (filt_t[i] >= bound1 and filt_t[i] <= bound2) or (filt_t[i] >= bound3 and filt_t[i] <= bound4) or (filt_t[i] >= bound5 and filt_t[i] <= bound6):
        ffilt_t.append(filt_t[i])
        ffilt_y.append(filt_y[i])

#ffilt_t = filt_t
#ffilt_y = filt_y

phase_bins = np.linspace(0,2,20*2+1)

normal_pgram = signal.lombscargle(t,y,f,normalize=True)
pgram = signal.lombscargle(ffilt_t,ffilt_y,f,normalize=True)

###
normal_phases = foldAt(t,2*np.pi/3,T0=0)
index_sort = np.argsort(normal_phases)
normal_phases = list(normal_phases[index_sort]) + list(normal_phases[index_sort]+1)
normal_fluxes = list(y[index_sort]) + list(y[index_sort])

normal_profile, bin_edges, binnumber = stats.binned_statistic(normal_phases,normal_fluxes,statistic='mean',bins=phase_bins)
###

###
phases = foldAt(np.array(ffilt_t),2*np.pi/3,T0=0)
index_sort = np.argsort(phases)
phases = list(phases[index_sort]) + list(phases[index_sort]+1)
fluxes = list(np.array(ffilt_y)[index_sort]) + list(np.array(ffilt_y)[index_sort])

profile, bin_edges, binnumber = stats.binned_statistic(phases,fluxes,statistic='mean',bins=phase_bins)
###

print(normal_profile)
print(profile)

###
#sr_ph,sr_profile,sr_profile_err = fold_events(y,3/(2*np.pi),nbin=20)
###

plt.figure()
plt.plot(t,y,'bx')
plt.plot(ffilt_t,ffilt_y,'rx')
plt.figure()
plt.plot(f,pgram,'r-')
plt.plot(f,normal_pgram,'b-')
plt.legend(('Altered','Original'),loc='best')
plt.figure()
plt.step(phase_bins[:-1],normal_profile,'b-')
plt.step(phase_bins[:-1],profile,'r-')
#plt.step(sr_ph,sr_profile/1000-0.5,'k-')
plt.legend(('No gaps','Gaps','Stingray (no gaps)'),loc='best')
plt.show()
"""

"""
eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070102_pipe/ni0034070102_nicersoft_bary.evt'
times = fits.open(eventfile)[1].data['TIME']
gtis = Lv0_fits2dict.fits2dict(eventfile,2,['START','STOP'])
T = sum([ (gtis['STOP'][i]-gtis['START'][i]) for i in range(len(gtis['START'])) ])
trunc_t = times - times[0]

sr_ph,sr_profile,sr_profile_err = fold_events(trunc_t,0.2081,nbin=20,expocorr=True)
plt.step(sr_ph,sr_profile,'k-')
_ = plot_profile(sr_ph,sr_profile*T/(trunc_t[-1]))
plt.show()
"""

"""
nin = 1000
nout = 100000
x = np.linspace(0.01,10*np.pi,nin)
y = 2*np.sin(1*x)
f = np.linspace(0.01,10,nout)
pgram = signal.lombscargle(x,y,f,normalize=True)
plt.figure()
plt.plot(x,y,'bx-')
plt.figure()
plt.plot(f,pgram,'r-')
plt.show()
"""

"""
basefolder = '/Volumes/Samsung_T5/n300_ulx_2020/'
obsids = ['1034200' + str(i) for i in range(101,200)] + ['1034200' + str(i) for i in range(201,241)] + ['2034200' + str(i) for i in range(201,206)]
print(len(obsids))
"""

"""
ngc300x1 = fits.open('/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge.evt')[1].data['PI']
bins = np.array([20,100,200,300,400,500,600,700,800,900,1000])
plt.hist(ngc300x1*10/1000,bins*10/1000,log=True)
plt.xlabel('Energy (keV)',fontsize=12)
plt.ylabel('Number of events',fontsize=12)
plt.show()
"""

"""
eventlist = open('/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/eventfiles.list','r').read().split('\n')[:-1]
E_boundary = 2.0 #energy boundary in keV

mjds = []
color = []
color_err = []

for i in tqdm(range(len(eventlist))):
    eventheader = fits.open(eventlist[i])[1].header
    mjd = eventheader['MJDREFI'] + eventheader['MJDREFF'] + (eventheader['TSTART']+eventheader['TSTOP'])/(2*86400)

    pis = fits.open(eventlist[i])[1].data['PI']
    if len(pis) != 0:
        soft = len(pis[pis*10/1000<E_boundary])
        hard = len(pis[pis*10/1000>=E_boundary])
        soft_err = np.sqrt(soft)
        hard_err = np.sqrt(hard)
        if soft != 0:
            color.append(hard/soft)
            error = np.sqrt( (hard_err/soft)**2 + (hard*soft_err/soft**2)**2 )
            color_err.append(error)
            mjds.append(mjd)

mjds = np.array(mjds)
color = np.array(color)
color_err = np.array(color_err)
print(mjds[0],mjds[-1])

plt.figure()
plt.errorbar(x=mjds,y=color,yerr=color_err,fmt='rx-')
plt.xlabel('Time (MJD)',fontsize=12)
plt.ylabel('Color (H/S) - (2.0-10 keV)/(0.2-2.0 keV)',fontsize=12)
plt.show()
"""

"""
plt.figure()
circle1 = plt.Circle((0.5,0.5),radius=0.2,color='r',lw=0.5,label='NICER FOV')
plt.gcf().gca().add_artist(circle1)
plt.show()
"""

"""
energy = np.linspace(100,12000,100001)
N = 8.7
fano = 0.114
w = 3.71
fwhm = 2.35*w*(N**2+fano*energy/w)**(1/2)
plt.plot(energy,fwhm,'r-')
plt.show()
"""

"""
checkfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge_niceroverlap2_spec2.evt'
times = fits.open(checkfile)[1].data['TIME']
gtis = fits.open(checkfile)[2].data
startt = gtis['START']
stopt = gtis['STOP']
checkdict = {}

for i in range(len(times)):
    if times[i] == times[i-1]:
        print(times[i])
    checkdict[times[i]] = []
    for j in range(len(startt)):
        if (times[i]>=startt[j]) and (times[i]<=stopt[j]):
            checkdict[times[i]].append(startt[j])

dictkeys = list(checkdict.keys())
print(len(dictkeys))

#for i in range(len(dictkeys)):
#    if len(checkdict[dictkeys[i]]) == 0:
#        print(dictkeys[i])
"""

"""
import numpy.ma as ma
from reproject import reproject_interp

testimage = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/products/sw00049834027xpc_sk.img'
expimage = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/products/sw00049834027xpc_ex.img'
testimage2 = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/products/sw00049834028xpc_sk.img'

obsid = str(pathlib.Path(testimage).name)[:13]
obsid2 = str(pathlib.Path(testimage2).name)[:13]

fitsfile = fits.open(testimage)[0]
expfile = fits.open(expimage)[0]
fitsfile2 = fits.open(testimage2)[0]

wcs = WCS(fitsfile.header)
wcs2 = WCS(fitsfile2.header)

plt.figure()
plt.subplot(projection=wcs)
a = plt.imshow(fitsfile.data,vmin=0,vmax=np.nanmax(fitsfile.data),cmap='gist_heat')
plt.colorbar()

plt.figure()
plt.subplot(projection=wcs)
b = plt.imshow(fitsfile.data/expfile.data,vmin=0,vmax=np.nanmax(fitsfile.data/expfile.data),cmap='gist_heat')
plt.colorbar()

print(np.max(fitsfile.data/expfile.data))
"""

"""
plt.figure()
plt.subplot(projection=wcs)
b = plt.imshow(fitsfile2.data,vmin=0,vmax=np.max(fitsfile2.data),cmap='gist_heat')


getarray = np.array(a.get_array())
getarray2 = np.array(b.get_array())

plt.figure()
plt.imshow(getarray+getarray2,vmin=0,vmax=np.max(getarray+getarray2),cmap='gist_heat')

plt.figure()
plt.subplot(projection=wcs)
array,footprint = reproject_interp(fitsfile2,fitsfile.header)
print(np.nanmax(array))
plt.imshow(fitsfile.data+array,vmin=0,cmap='gist_heat')
"""

#plt.show()


"""
input_E = np.array([0.5,0.6,0.7,0.8,0.9] + list(np.arange(1.0,12.5,0.5)))
pi_ev = input_E * 100
fwhm = []
fwhm_err = []

for i in range(len(pi_ev)):

    def gauss(x,a,sig,constant):
        return a*np.exp( -(x-input_E[i])**2/(2*sig**2) ) + constant

    #testspec = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/find_E_res/gauss_' + str(int(pi_ev[i])).zfill(4) + '.txt'
    testspec = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/gauss_' + str(int(pi_ev[i])).zfill(5) + '.txt'
    print(testspec)
    E,E_err,flux,flux_err,model = np.genfromtxt(testspec,usecols=(0,1,2,3,4),skip_header=3,skip_footer=1,unpack=True)

    pguess = np.array([100,0.05,0.01])
    popt,pcov = curve_fit(gauss,E,model,p0=pguess)

    print('Gaussian width is ' + str(round(popt[1]*1000,2)) + ' eV; FWHM is ' + str(round(popt[1]*1000*2.355,2)) + ' +- ' + str(round(np.sqrt(np.diag(pcov)[1])*1000*2.355,2)) + ' eV')
    fwhm.append(popt[1]*1000*2.355)
    fwhm_err.append(np.sqrt(np.diag(pcov)[1])*1000*2.355)
    #plt.plot(E,model,'b-')
    #plt.plot(E,gauss(E,popt[0],popt[1],popt[2]),'r-')
    #plt.show()

timeend = time.time()


def sqrtfunc(x,a,b,c):
    return a*(x+b)**c

#def parabola(x,a,b,c):
#    return a*(x-)

popt,pcov = curve_fit(sqrtfunc,input_E,fwhm,sigma=fwhm_err,p0=[10,3.3,0.5])
print(popt)
print(np.sqrt(np.diag(pcov)))
#for i in range(len(input_E)):
#    print(input_E[i],sqrtfunc(input_E,popt[0],popt[1],popt[2],popt[3])[i]/3)

N = 8.7
fano = 0.114
w = 3.71

nicer_spie = 2.35*w*(N**2 + input_E*1000*fano/w)**(1/2)

plt.plot(input_E,fwhm,'rx-')
plt.plot(input_E,sqrtfunc(input_E,popt[0],popt[1],popt[2]),'b-')
plt.plot(input_E,nicer_spie,'kx-')
plt.annotate(str(round(popt[0],2))+'(E+'+str(round(popt[1],2))+')^'+str(round(popt[2],2)),(2,160))
plt.xlabel('Energy, E (keV)',fontsize=12)
plt.ylabel('FWHM (eV)',fontsize=12)
plt.legend(('Measured','curve_fit','NICER SPIE'),fontsize=12)
plt.show()
"""

"""
MJDs_nicer = np.array([58239,58244,58249,58254,58259,58264,58269,58274,58309,58314,58324,58329,58334,58339,58389,58449,58454,58484,58489,58504,58509])
flux_nicer = np.array([3.93e-12,4.79e-12,3.28e-12,3.01e-12,3.00e-12,3.05e-12,3.17e-12,2.80e-12,2.36e-12,2.70e-12,2.05e-12,1.73e-12,1.95e-12,2.48e-12,1.36e-12,1.17e-12,2.21e-12,2.85e-12,1.65e-12,1.57e-12,9.21e-13])

plt.semilogy(MJDs_nicer,flux_nicer,'rx-')
#plt.axhline(y=1.36e-13,color='k',xmax=0.18)
#plt.axhline(y=1.55e-13,color='k',xmin=0.2)

#plt.axhline(y=1.83e-13,color='b',xmax=0.18)
#plt.axhline(y=1.35e-13,color='b',xmin=0.2)

#plt.axhline(y=1.37e-12,color='m',xmax=0.18)
#plt.axhline(y=3.41e-13,color='m',xmin=0.2)

#plt.legend(('NICER flux','NGC300 X-1','NGC300 X-1','NGC300 bg','NGC300 bg','NGC300 ULX-1','NGC300 ULX-1'),loc='best')

plt.axhline(y=1.42e-13,color='k')
plt.axhline(y=1.59e-13,color='b')
plt.axhline(y=8.40e-13,color='m')
plt.legend(('NICER flux','NGC300 X-1','NGC300 bg','NGC300 ULX-1'),loc='best')

plt.xlabel('Time (MJD)',fontsize=12)
plt.ylabel('Flux (ergs/s/cm^2)',fontsize=12)
plt.show()
"""

"""
mjds = ['58239','58244','58249','58254','58259','58264','58269','58274','58279','58284','58289','58294','58309','58314','58324','58329','58334','58339','58344','58349','58384','58389','58394','58399','58409','58449','58454','58459','58464','58469','58474','58479','58484','58489','58494','58499','58504','58509','58604']
textfile = '/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/fparkey_jsgrp_05d.go'
writing = open(textfile,'w')
for i in range(len(mjds)):
    writing.write('fparkey' + ' grp_' + mjds[i] + '_05d_bg_cl50.pha ' +  'jsgrp_' + mjds[i] + '_05d_cl_cl50.pha ' + 'BACKFILE' + '\n')
writing.close()
"""

"""
mjds = ['58239','58244','58249','58254','58259','58264','58269','58274','58279','58284','58289','58294','58309','58314','58324','58329','58334','58339','58344','58349','58384','58389','58394','58399','58409','58449','58454','58459','58464','58469','58474','58479','58484','58489','58494','58499','58504','58509','58604']
cl_bg_file = '/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/getspec_cl_bg.xcm'
basefolder = '/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/'
writing = open(cl_bg_file,'w')
writing.write('data ')
base_string = ''
for i in range(len(mjds)):
    base_string += str(i+1) + ':' + str(i+1) + ' ' + basefolder + 'jsgrp_' + mjds[i] + '_05d_cl_cl50.pha '
writing.write(base_string + '\n')

#for i in range(len(mjds)):
#    writing.write('backgrnd ' + str(i+1) + ' grp_' + mjds[i] + '_05d_bg_cl50.pha' + '\n')

writing.write('setplot energy' + '\n')
#writing.write('setplot background on' + '\n')
writing.write('ignore **:0.0-0.285,12.01-**')
writing.close()
"""

"""
def sinfunc(t,a,b,c,d):
    return a*np.sin(d*2*np.pi*t+b) + c

eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/lightcurve/ngc300x1_100s.lc'
fitsfile = fits.open(eventfile)[1].data

times = fitsfile['TIME'][:1000]
rate = fitsfile['RATE'][:1000]

popt,pcov = curve_fit(sinfunc,times,rate,p0=np.array([0.05,0.1,0.05,8.46e-6]))
print(popt)
print(np.sqrt(np.diag(pcov)))

plt.plot(times,rate,'rx-')
plt.plot(times,sinfunc(times,popt[0],popt[1],popt[2],popt[3]),'b-')
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.show()
"""


"""
import Lv2_ps_method

cenx3 = fits.open('/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/trylc/haha.lc')[1].data
times = cenx3['TIME']
rates = cenx3['RATE']
errors = cenx3['ERROR']
fracexp = cenx3['FRACEXP']

plt.figure()
T = times[fracexp==1][-1] - times[fracexp==1][0]
dt = T/len(times[fracexp==1])
freq = 1/dt

dt = 0.5
freq = 2
f,pxx = signal.periodogram(rates[fracexp==1],fs=freq)
print(T,dt,freq)
#pdgm_f,pdgm_ps = Lv2_ps_method.pdgm(times,rates,[False,0,10],[False,2],True,[False,5])

#plt.plot(pdgm_f,pdgm_ps,'b-')
plt.plot(f,pxx,'b-')

#plt.figure()
#plt.errorbar(x=times[fracexp==1],y=rates[fracexp==1],yerr=errors[fracexp==1],fmt='bx-')
plt.show()
"""

"""
#eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/ni0034070101_nicersoft_bary.evt'
eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge_niceroverlap_all.evt'
#eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge.evt'

times = fits.open(eventfile)[1].data['TIME']
gtis_data = fits.open(eventfile)[2].data
T = sum([ gtis_data[i]['STOP']-gtis_data[i]['START'] for i in range(len(gtis_data)) ])
#T0 = 56658 + 0.000777592592592593 + (fits.open(eventfile)[1].header['TSTART'] + (-1))/86400
T0 = 51910 + 7.428703700000000E-04 + fits.open(eventfile)[1].header['TSTART']/86400
#T0 = 58244.6660

#plt.figure()
#plt.title('NGC300x1_merge_niceroverlap_all.evt')
#plt.xlabel('Phase',fontsize=12)
#plt.ylabel('Counts/s',fontsize=12)

#phase,profile,profile_error = Lv2_phase.pulse_folding(times,T,T0,8.466777057846186e-06,0,0,20)
#phase,profile,profile_error = Lv2_phase.pulse_folding(times,T,T0,8.464749477590145e-06,0,0,20)
#phase,profile,profile_error = Lv2_phase.pulse_folding(times,T,T0,8.47122e-6,0,0,20)
### 8.461949E-6 is from chi^2 ; 8.46185E-6 is from Z^2, 8.465941706160968E-6 is from FR/RSS

#plt.step(phase[:-1],profile,'r-')
#plt.errorbar(x=phase[:-1],y=profile,yerr=profile_error,color='r',drawstyle='steps-mid')
gtis_conform = []
for i in range(len(gtis_data)):
    gtis_conform.append([gtis_data[i][0],gtis_data[i][1]])
"""

"""
expos = phase_exposure(times[0]-times[0],times[-1]-times[0],1/8.47122e-6,nbin=20,gtis=gtis_conform-times[0])
print(expos)
total_expos = np.array(list(expos) + list(expos))
#plt.figure()
plt.errorbar(x=phase[:-1],y=profile/total_expos,yerr=profile_error/total_expos,color='b',drawstyle='steps-mid')
plt.title('NGC300x1_merge_niceroverlap_all.evt, exposure-corrected')
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.legend(('Folded profile','Exposure-corrected profile'),loc='best',fontsize=12)

#print(times[0],times[-1])
#t0 = (times[0]+times[-1])/2+0.175*1/8.47122e-06
#print(t0)
t0 = times[0]

phase_sr,prof_sr,err_sr = fold_events(times,[8.47122e-06],gtis=np.array(gtis_conform),ref_time=t0,nbin=20)
phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times,[8.47122e-6],gtis=np.array(gtis_conform),ref_time=t0,expocorr=True,nbin=20)

print(prof_sr/prof_sr_expo)

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
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)

#print(prof_sr)
#print(sum(prof_sr))
#print(prof_sr_expo)
#print(sum(prof_sr_expo))

#print(sum(profile*T))
#print(sum(profile/total_expos*T))

plt.show()
"""

"""
chi2 = []
freqs = np.arange(8.0e-6,9.0e-6,0.0001e-6)
for i in tqdm(range(len(freqs))):
    phase_sr_expo,prof_sr_expo,err_sr_expo = fold_events(times,[freqs[i]],gtis=np.array(gtis_conform),ref_time=times[0],expocorr=True)
    chi2.append( Lv2_phase.get_chi2(prof_sr_expo,err_sr_expo) )

plt.figure()
plt.plot(freqs,chi2,'rx-')
plt.xlabel('Frequency (Hz)',fontsize=12)
plt.ylabel('chi^2 [ sum( (profile-mean)^2/error^2) ]',fontsize=12)
plt.show()
"""

"""
medians = np.array([8.335751931919311,8.33584176174108,8.337165146653654,8.337165146653544,
                    8.337652596774758,8.337623326928423,8.337600802469019,8.338121728007296,
                    8.338676789004637,8.338070034352355,8.338064677098994])

median_err = np.array([0.015450388314367997,0.014905023834124513,0.016690004956593195,0.011717682240614905,
                    0.01499435360625864,0.011451324497186582,0.014450108614211365,0.01536640919810147,
                    0.01590570701944077,0.01510912946652743,0.017427593377916324])

pdgm = np.array([8.335234336108602,8.335234336108602,8.337495697570686,8.337230902833188,
                8.337931466792437,8.3384201347159,8.337091556933076,8.337865199746701,
                8.338395922865838,8.337507666407173,8.338676789004637])

expected = (1/120e3)/1e-6

stddev = np.abs(medians-pdgm)/median_err
stddev_exp = np.abs(medians-expected)/median_err

print(stddev)
print(np.median(stddev),np.mean(stddev),np.std(stddev))
print(stddev_exp)
print(np.median(stddev_exp),np.mean(stddev_exp),np.std(stddev_exp))
"""

"""
basefolder = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/'
obsids = [str(i) for i in range(49834027,49834042)] + [str(i) for i in range(49834043,49834062)] + [str(i) for i in range(49834063,49834066)] + ['88810002'] + [str(i) for i in range(49834066,49834069)] + [str(i) for i in range(49834070,49834079)] + [str(i) for i in range(49834080,49834088)]
counter = 0
for i in range(len(obsids)):
    fitsfile = fits.open(basefolder + 'sw000' + obsids[i] + 'xpcw3po_bary_cl_ngc300x1_ngc300x1.evt')
    print(len(fitsfile[2].data))
    counter += len(fitsfile[2].data)

print(counter)
"""

"""
##### eclipse timings - 6/24, see the spreadsheet

#eclipse_times = np.array([547432660.5728,548500113.9923,548730242.8801,549555729.5291,
#                        550037412.9875,552279376.5862,552747356.4104,553577945.8216,557242679.2037,
#                        559373193.2472,569064024.9780]) #8.461949e-6

eclipse_times = np.array([547432660.5728,548500113.9923,548730242.8801,
                        549555729.5291,550037412.9875,552279376.5862,552747356.4104,553577945.8216,
                        554042992.9464,557235752.3213,567870701.9700]) #8.464749e-6

MJDREFI = 51910
MJDREFF = 7.428703700000000E-04
mjds = MJDREFI + MJDREFF + eclipse_times/86400

#no_cycles = np.array([1,10,12,19,23,42,46,53,84,102,184]) #8.461949e-6
no_cycles = np.array([1,10,12,19,23,42,46,53,57,84,174]) #8.464749e-6

T0 = MJDREFI + MJDREFF + 547313659.6723/86400

def linear_mod(N,x0,a):
    return x0 + a*N

def quad_mod(N,x0,a,b):
    return x0 + a*N + b*N**2

popt,pcov = curve_fit(linear_mod,no_cycles,mjds,p0=[T0,1.5])
popt_q,pcov_q = curve_fit(quad_mod,no_cycles,mjds,p0=[T0,1.5,1e-5])

print('T0 from linear: ' + str(popt[0]))
print('T0 from quadratic: ' + str(popt_q[0]))

print( (popt[0] - MJDREFI - MJDREFF) * 86400 )
print( (popt[0] - MJDREFI - MJDREFF) * 86400 - 547313659.6723 )
print( (popt_q[0] - MJDREFI - MJDREFF) * 86400 )
print( (popt_q[0] - MJDREFI - MJDREFF) * 86400 - 547313659.6723 )

print('Linear:')
print(popt)
print(np.diag(np.sqrt(pcov)))
print('Quadratic:')
print(popt_q)
print(np.diag(np.sqrt(pcov_q)))

plt.figure()
plt.plot(no_cycles,mjds,'rx-')
plt.plot(no_cycles,linear_mod(no_cycles,popt[0],popt[1]),'bx-')
plt.figure()
plt.plot(no_cycles,mjds-linear_mod(no_cycles,popt[0],popt[1]))
plt.figure()
plt.plot(no_cycles,quad_mod(no_cycles,popt_q[0],popt_q[1],popt_q[2]))
plt.figure()
plt.plot(no_cycles,mjds-quad_mod(no_cycles,popt_q[0],popt_q[1],popt_q[2]))
plt.show()
"""

"""
nsample = 100
T = 480
x = np.arange(0,86400,nsample)
y = 0.05*np.sin(2*np.pi/T*x) + np.random.normal(loc=0,scale=0.05,size=len(x))
#plt.plot(x,y,'r-')

omega,psd,prob3,prob4,prob5 = Lv2_dj_lsp.lsp(x,y)
nu_reg = omega/(2.0*np.pi)
freq = omega/(2*np.pi)

plt.axhline(y=prob3,lw=0.5,alpha=0.5)
plt.axhline(y=prob4,lw=0.5,alpha=0.5)
plt.axhline(y=prob5,lw=0.5,alpha=0.5)

plt.plot(freq,psd,'rx-')

plt.show()
"""

"""
eventfiles = open('/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/niceroverlap_all.list','r').read().split('\n')[:-1]
expos = []
for i in range(len(eventfiles)):
    gtis = fits.open(eventfiles[i])[2].data
    total_gti = sum( [gtis[i][1] - gtis[i][0] for i in range(len(gtis)) ])
    expos.append(total_gti)

expos = np.array(expos)
print(np.mean(expos))
print(np.median(expos))
print(np.std(expos))
print(np.sum(expos))
plt.hist(expos)
plt.xlabel('Exposure times (s)',fontsize=12)
plt.ylabel('Number',fontsize=12)
plt.show()
"""

"""
popt,pcov = curve_fit(linear_mod,no_cycles,mjds,p0=[1.5])
popt_q,pcov_q = curve_fit(quad_mod,no_cycles,mjds,p0=[1.5,1e-5])
print('Linear:')
print(popt)
print(np.diag(np.sqrt(pcov)))
print('Quadratic:')
print(popt_q)
print(np.diag(np.sqrt(pcov_q)))

plt.figure()
plt.plot(no_cycles,mjds,'rx-')
plt.plot(no_cycles,linear_mod(no_cycles,popt[0]),'bx-')
plt.figure()
plt.plot(no_cycles,mjds-linear_mod(no_cycles,popt[0]))
plt.figure()
plt.plot(no_cycles,quad_mod(no_cycles,popt_q[0],popt_q[1]))
plt.figure()
plt.plot(no_cycles,mjds-quad_mod(no_cycles,popt_q[0],popt_q[1]))
plt.show()
"""

"""
def phase_exposure(start_time, stop_time, period, nbin=16, gtis=None):

    Calculate the exposure on each phase of a pulse profile.

    Parameters
    ----------
    start_time, stop_time : float
        Starting and stopping time (or phase if ``period==1``)
    period : float
        The pulse period (if 1, equivalent to phases)

    Other parameters
    ----------------
    nbin : int, optional, default 16
        The number of bins in the profile
    gtis : [[gti00, gti01], [gti10, gti11], ...], optional, default None
        Good Time Intervals

    Returns
    -------
    expo : array of floats
        The normalized exposure of each bin in the pulse profile (1 is the
        highest exposure, 0 the lowest)

    if gtis is None:
        gtis = np.array([[start_time, stop_time]])

    # Use precise floating points -------------

    start_time = np.longdouble(start_time)
    stop_time = np.longdouble(stop_time)
    period = np.longdouble(period)
    gtis = np.array(gtis, dtype=np.longdouble)
    # -----------------------------------------
    #print('GTIs from longdouble: ' + str(sum([gtis[i][1]-gtis[i][0] for i in range(len(gtis))])))

    expo = np.zeros(nbin)
    phs = np.linspace(0, 1, nbin + 1)
    phs = np.array(list(zip(phs[0:-1], phs[1:])))

    # Discard gtis outside [start, stop]
    good = np.logical_and(gtis[:, 0] < stop_time, gtis[:, 1] > start_time)
    gtis = gtis[good]

    for g in gtis:
        g0 = g[0]
        g1 = g[1]
        if g0 < start_time:
            # If the start of the fold is inside a gti, start from there
            g0 = start_time
        if g1 > stop_time:
            # If the end of the fold is inside a gti, end there
            g1 = stop_time
        length = g1 - g0
        # How many periods inside this length?
        nraw = length / period
        # How many integer periods?
        nper = nraw.astype(int)

        # First raw exposure: the number of periods
        expo += nper / nbin

        # FRACTIONAL PART =================
        # What remains is additional exposure for part of the profile.
        start_phase = np.fmod(g0 / period, 1)
        end_phase = nraw - nper + start_phase
        limits = [[start_phase, end_phase]]
        # start_phase is always < 1. end_phase not always. In this case...
        if end_phase > 1:
            limits = [[0, end_phase - 1], [start_phase, 1]]

        for l in limits:
            l0 = l[0]
            l1 = l[1]
            # Discards bins untouched by these limits
            goodbins = np.logical_and(phs[:, 0] <= l1, phs[:, 1] >= l0)
            idxs = np.arange(len(phs), dtype=int)[goodbins]
            for i in idxs:
                start = np.max([phs[i, 0], l0])
                stop = np.min([phs[i, 1], l1])
                w = stop - start
                expo[i] += w

    #print('Expo*period: ' + str(sum(expo)*period))
    return expo / np.max(expo)
"""

"""
#eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/0034070101_pipe/ni0034070101_nicersoft_bary.evt'
eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge_niceroverlap_all.evt'

times = fits.open(eventfile)[1].data['TIME']
gtis_data = fits.open(eventfile)[2].data
T = sum([ gtis_data[i]['STOP']-gtis_data[i]['START'] for i in range(len(gtis_data)) ])
#T0 = 56658 + 0.000777592592592593 + (fits.open(eventfile)[1].header['TSTART'] + (-1))/86400
T0 = 51910 + 7.428703700000000E-04 + fits.open(eventfile)[1].header['TSTART']/86400

gtis_conform = []
for i in range(len(gtis_data)):
    gtis_conform.append([gtis_data[i][0],gtis_data[i][1]])

#print('T: ' + str(T))

expos = phase_exposure(times[0],times[-1],1/8.4666e-6,nbin=20,gtis=gtis_conform)
print(expos)
total_expos = np.array(list(expos) + list(expos))
"""

"""
t1 = 5
t2 = 35
GTIs1 = np.array([ [3,6], [10,19], [20,25], [26,32], [33,36] ])
nbin = 20
expos = phase_exposure(t1,t2,3,nbin=20,gtis=GTIs1)
print(expos)

gti_phases = pulse_phase(GTIs1,[1/3],to_1=False)
start_phase,stop_phase = pulse_phase(np.array([t1,t2]),[1/3],to_1=False)
expo_norm = phase_exposure(start_phase,stop_phase,1,20,gtis=gti_phases)
print(expo_norm)
"""


"""
f1 = [8.465068630012935e-06, 8.46813510807163e-06, 8.464302689232871e-06, 8.466605940469421e-06, 8.467369167046e-06, 8.468898333913252e-06, 8.467369167046e-06, 8.464302689232871e-06, 8.467366452353124e-06, 8.465834570792998e-06, 8.465065916058498e-06, 8.466597797127494e-06, 8.46277894733641e-06, 8.467363737661991e-06, 8.470432931148515e-06, 8.467363737661991e-06, 8.464305402943485e-06, 8.465834570792998e-06, 8.466600511573061e-06, 8.465834570792998e-06, 8.467363737661991e-06, 8.465831856592996e-06, 8.47119615625344e-06, 8.470432931148515e-06, 8.471193440334479e-06, 8.467363737661991e-06, 8.46812967819649e-06, 8.468895618730987e-06, 8.465065916058498e-06, 8.467363737661991e-06, 8.468901049097257e-06, 8.46584271340345e-06, 8.466597797127494e-06, 8.466597797127494e-06, 8.465065916058498e-06, 8.472725321403475e-06, 8.467366452353124e-06, 8.463547602323434e-06, 8.465839999198226e-06, 8.468898333913252e-06, 8.469664274693314e-06, 8.467366452353124e-06, 8.470432931148515e-06, 8.465831856592996e-06, 8.468132393133189e-06, 8.464305402943485e-06, 8.465831856592996e-06, 8.464305402943485e-06, 8.466600511573061e-06, 8.467363737661991e-06]

f2 = [8.468895618730987e-06, 8.468895618730987e-06, 8.46507677188669e-06, 8.46812967819649e-06, 8.466597797127494e-06, 8.465834570792998e-06, 8.468132393133189e-06, 8.464302689232871e-06, 8.4696697055542e-06, 8.465068630012935e-06, 8.466597797127494e-06, 8.464302689232871e-06, 8.462768094455006e-06, 8.466597797127494e-06, 8.463536748452808e-06, 8.465079485848087e-06, 8.468132393133189e-06, 8.469664274693314e-06, 8.469661559265485e-06, 8.466597797127494e-06, 8.465068630012935e-06, 8.46431083036993e-06, 8.47119615625344e-06, 8.466600511573061e-06, 8.464302689232871e-06, 8.465071343969115e-06, 8.47503129445256e-06, 8.463534034989504e-06, 8.46812967819649e-06, 8.469661559265485e-06, 8.46123892611262e-06, 8.465834570792998e-06, 8.46507677188669e-06, 8.466605940469421e-06, 8.468903764283003e-06]

f3 = [8.46812967819649e-06, 8.468132393133189e-06, 8.465068630012935e-06, 8.466597797127494e-06, 8.465834570792998e-06, 8.465065916058498e-06, 8.468895618730987e-06, 8.462776234113448e-06, 8.468901049097257e-06, 8.464302689232871e-06, 8.46812967819649e-06, 8.462002153920508e-06, 8.468132393133189e-06, 8.46660322602037e-06, 8.463539461917857e-06]

total_f = f1+f2+f3
print(len(total_f))
print(np.std(np.array(total_f)))
print(np.mean(np.array(total_f)))
"""

"""
import numpy as np
import numpy.polynomial.polynomial as npoly
from scipy import optimize
import matplotlib.pyplot as plt
np.random.seed(2017)

def f(breakpoints, x, y, fcache):
    breakpoints = tuple(map(int, sorted(breakpoints)))
    if breakpoints not in fcache:
        total_error = 0
        for f, xi, yi in find_best_piecewise_polynomial(breakpoints, x, y):
            total_error += ((f(xi) - yi)**2).sum()
        fcache[breakpoints] = total_error
    # print('{} --> {}'.format(breakpoints, fcache[breakpoints]))
    return fcache[breakpoints]

def find_best_piecewise_polynomial(breakpoints, x, y):
    breakpoints = tuple(map(int, sorted(breakpoints)))
    xs = np.split(x, breakpoints)
    ys = np.split(y, breakpoints)
    result = []
    for xi, yi in zip(xs, ys):
        if len(xi) < 2: continue
        coefs = npoly.polyfit(xi, yi, 1)
        f = npoly.Polynomial(coefs)
        result.append([f, xi, yi])
    return result

x = np.arange(0.45,1.95,0.05)
y = np.array([0.00105478,0.00108443,0.00104594,0.00126743,0.00125505,0.00110672,
        0.00106399,0.00094255,0.00111247,0.00120255,0.00122203,0.00125356,
        0.0010629,0.00050578,0.00029983,0.00033695,0.00036172,0.00035548,
        0.00056006,0.00064167,0.00105478,0.00108443,0.00104594,0.00126743,
        0.00125505,0.00110672,0.00106399,0.00094255,0.00111247,0.00120255])

num_breakpoints = 4
breakpoints = optimize.brute(
    f, [slice(1, len(x), 1)]*num_breakpoints, args=(x, y, {}), finish=None)

plt.scatter(x, y, c='blue', s=50)
for f, xi, yi in find_best_piecewise_polynomial(breakpoints, x, y):
    x_interval = np.array([xi.min(), xi.max()])
    print('y = {:35s}, if x in [{}, {}]'.format(str(f), *x_interval))
    plt.plot(x_interval, f(x_interval), 'ro-')

plt.show()
"""

"""
xmm_lc1 = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_lccorr.lc'
xmm_lc2 = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_lccorr.lc'

times1 = fits.open(xmm_lc1)[1].data['TIME']
rates1 = fits.open(xmm_lc1)[1].data['RATE']
times2 = fits.open(xmm_lc2)[1].data['TIME']
rates2 = fits.open(xmm_lc2)[1].data['RATE']

total_time = np.array(list(fits.open(xmm_lc1)[1].header['MJDREF'] + times1) + list(fits.open(xmm_lc2)[1].header['MJDREF'] + times2) )
total_rates = np.array(list(rates1) + list(rates2))

#print(sum(100*total_rates))


trunc_times = total_time-total_time[0]
"""
#t_bins = np.arange(trunc_times[0],trunc_times[-1]+100,100)
#summed_data, bin_edges, binnumber = stats.binned_statistic(trunc_times,total_rates,statistic='mean',bins=t_bins) #binning the time values in the data

#freqs = np.linspace(1e-7,1e-4,10001)
"""
plt.figure()
plt.plot(total_time,total_rates,'r-')
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Rate (counts/s)',fontsize=12)
plt.title('XMM Newton light curve for NGC 300 X-1 (0791010101 and 0791010301)')
plt.show()


plt.figure()

import scipy.signal as signal
pgram = signal.lombscargle(total_time, total_rates, freqs, normalize=True)
plt.plot(freqs,pgram,'rx-')

plt.figure()
pgram2 = signal.lombscargle(t_bins[:-1],summed_data,freqs,normalize=True)
plt.plot(freqs,pgram2,'bx-')
print(freqs[:10],pgram2[:10])
"""

"""
plt.figure()
omega,psd,prob3,prob4,prob5 = Lv2_dj_lsp.lsp(trunc_times,total_rates)
freq = omega/(2*np.pi)

plt.plot(freq,psd,'rx-')
plt.show()
"""

"""
xmm_evt1 = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_bary_ngc300x1.evt'
xmm_evt1_bg = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_bary_ngc300bg.evt'
xmm_evt2 = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_bary_ngc300x1.evt'
xmm_evt2_bg = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_bary_ngc300bg.evt'

times1 = fits.open(xmm_evt1)[1].data['TIME']
counts1 = fits.open(xmm_evt1)[1].data['COUNTS']
times2 = fits.open(xmm_evt2)[1].data['TIME']
counts2 = fits.open(xmm_evt2)[1].data['COUNTS']
exp1 = fits.open(xmm_evt1)[1].header['EXPOSURE']
exp2 = fits.open(xmm_evt2)[1].header['EXPOSURE']

trunc_times1 = times1-times1[0]
trunc_times2 = times2-times1[0]

counts_bg1 = fits.open(xmm_evt1_bg)[1].data['COUNTS']
counts_bg2 = fits.open(xmm_evt2_bg)[1].data['COUNTS']

tbins = np.arange(trunc_times1[0],trunc_times2[-1]+100,100)

summed_data1, bin_edges, binnumber = stats.binned_statistic(trunc_times1,counts1,statistic='sum',bins=tbins)
summed_data2, bin_edges, binnumber = stats.binned_statistic(trunc_times2,counts2,statistic='sum',bins=tbins)

bg_summed_data1, bin_edges, binnumber = stats.binned_statistic(trunc_times1,counts_bg1,statistic='sum',bins=tbins)
bg_summed_data2, bin_edges, binnumber = stats.binned_statistic(trunc_times2,counts_bg2,statistic='sum',bins=tbins)

total_sum = summed_data1 + summed_data2
total_bg = bg_summed_data1 + bg_summed_data2



plt.figure()
plt.plot(tbins[:-1],total_sum,'k-')

plt.figure()
omega,psd,prob3,prob4,prob5 = Lv2_dj_lsp.lsp(tbins[:-1],total_sum)
freq = omega/(2*np.pi)

plt.plot(freq,psd,'rx-')
plt.show()
"""

"""
print(sum(total_sum),sum(total_bg))

plt.plot(tbins[:-1],total_sum,'k-')
plt.plot(tbins[:-1],total_bg*1/9,'r-')
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Counts/s',fontsize=12)
plt.title('Exposure time for 0791010101: ' + str(exp1) + ' ; exposure time for 0791010301: ' + str(exp2))
plt.show()
"""

"""
completeness = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
data = np.array([43.6,43.6,45.5,78.1,84.8,84.2,84.3,84.2,82.8,79.0,78.3])
sim1 = np.array([213.51,213.51,228.35,228.86,230.04,228.82,225.64,223.16,221.71,216.18,210.67])
sim2 = np.array([207.05,207.05,223.88,224.41,225.56,224.25,221.81,219.14,217.22,211.45,206.46])
sim3 = np.array([141.57,141.57,145.70,145.72,143.23,142.83,141.18,143.39,144.69,141.67,135.36])
sim4 = np.array([61.80,61.80,64.46,63.96,65.15,65.55,66.46,65.46,68.81,67.69,66.76])
sim5 = np.array([28.03,28.03,29.29,27.29,26.17,25.40,22.50,24.02,23.08,22.09,23.01])
sim6 = np.array([123.64,123.64,143.11,142.797,142.802,144.12,140.16,140.24,137.53,133.64,128.01])

plt.plot(completeness,data,'x-')
plt.plot(completeness,sim1,'x-')
plt.plot(completeness,sim2,'x-')
plt.plot(completeness,sim3,'x-')
plt.plot(completeness,sim4,'x-')
plt.plot(completeness,sim5,'x-')
plt.plot(completeness,sim6,'x-')

plt.xlabel('Completeness (fraction)',fontsize=12)
plt.ylabel('Normalized Power from LS periodogram',fontsize=12)
plt.legend(('Data','Sim 1','Sim 2','Sim 3','Sim 4','Sim 5','Sim 6'),loc='best')

plt.title('Sim. 1) No noise' + '\n' + 'Sim. 2) Noise level ~ 10% of amplitude of data' + '\n' + 'Sim. 3) Noise level ~ half of amplitude of data' + '\n' + 'Sim. 4) Noise level ~ same amplitude of data' + '\n' + 'Sim. 5) Noise level ~ 2x amplitude of data' + '\n' + 'Sim. 6) Low amplitude, and noise level ~ half amplitude of data')
plt.show()
"""

"""
freq,freqdot = np.meshgrid(np.arange(8.25e-6,8.7e-6,1e-8),np.arange(1e-18,9e-17,1e-20))
a = np.meshgrid(np.arange(8.25e-6,8.7e-6,1e-8),np.arange(1e-18,9e-17,1e-20))
freq_swift,freqdot_swift,chi2_swift = np.genfromtxt('/Volumes/Samsung_T5/NGC300_XMMdata/chi2_swift_freq_freqdot.txt',usecols=(0,1,2),unpack=True)
N_freq = len(np.arange(8.25e-6,8.7e-6,1e-8))
N_freqdot = len(np.arange(1e-18,9e-17,1e-20))
chi2_reshape = np.reshape(chi2_swift,(N_freqdot,N_freq),order='F')
"""

"""
print('0791010101 PN')
data1a = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_pn_clean.pha')
print(sum(data1a[1].data['COUNTS'][30:1000]))
bg1a = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_pn_bg.pha')
print(sum(bg1a[1].data['COUNTS'][30:1000]))

print('0791010101 MOS1')
data1b = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_mos1_clean.pha')
print(sum(data1b[1].data['COUNTS'][30:1000]))
bg1b = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_mos1_bg.pha')
print(sum(bg1b[1].data['COUNTS'][30:1000]))

print('0791010101 MOS2')
data1c = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_mos2_clean.pha')
print(sum(data1c[1].data['COUNTS'][30:1000]))
bg1c = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_mos2_bg.pha')
print(sum(bg1c[1].data['COUNTS'][30:1000]))

print('0791010301 PN')
data2a = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_pn_clean.pha')
print(sum(data2a[1].data['COUNTS'][30:1000]))
bg2a = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_pn_bg.pha')
print(sum(bg2a[1].data['COUNTS'][30:1000]))

print('0791010301 MOS1')
data2b = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_mos1_clean.pha')
print(sum(data2b[1].data['COUNTS'][30:1000]))
bg2b = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_mos1_bg.pha')
print(sum(bg2b[1].data['COUNTS'][30:1000]))

print('0791010301 MOS2')
data2c = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_mos2_clean.pha')
print(sum(data2c[1].data['COUNTS'][30:1000]))
bg2c = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_mos2_bg.pha')
print(sum(bg2c[1].data['COUNTS'][30:1000]))
"""

"""
eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/at2018cow_2020/at2018cow_filtered.evt'
eventfiles = glob.glob('/Volumes/Samsung_T5/NICERsoft_outputs/at2018cow_2020/accelsearch_GTIs/*.evt')
gtis_highsnr = fits.open(eventfile)[5].data
niextract_folder = '/Volumes/Samsung_T5/NICERsoft_outputs/at2018cow_2020/accelsearch_GTIs/'

ACCEL_files = glob.glob('/Volumes/Samsung_T5/NICERsoft_outputs/at2018cow_2020/accelsearch_GTIs/*ACCEL_100')
eventfiles_ACCEL = [ACCEL_files[i][:-10] + '.evt' for i in range(len(ACCEL_files))]
gti_no = [float(ACCEL_files[i][-14:-10]) for i in range(len(ACCEL_files))]

#with PdfPages('/Volumes/Samsung_T5/NICERsoft_outputs/at2018cow_2020/at2018cow_accelsearch.pdf') as pdf:

#header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
#header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

for i in tqdm(range(len(ACCEL_files))): #for each successful ACCEL_zmax file:
    freqs = []
    sigmas = []
    accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
    index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
    index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
    no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates
    for j in range(no_cands):
        freqs.append(float(accel_textfile[j+3].split()[6][:-3]))
        sigmas.append(float(accel_textfile[j+3].split()[1]))
    #times = np.ones(no_cands)*(fits.open(eventfiles_ACCEL[i])[1].header['TSTART']+fits.open(eventfiles_ACCEL[i])[1].header['TSTOP'])/2
    times = np.ones(no_cands)*gti_no[i]
    plt.scatter(x=times,y=freqs,c=sigmas,marker='o',cmap='gist_heat',vmin=2,vmax=4,edgecolors='k')

plt.colorbar().set_label('Significance (sigma)')

import mplcursors
mplcursors.cursor(hover=True)
plt.xlabel('Time',fontsize=12)
plt.ylabel('Frequency',fontsize=12)
#pdf.savefig()
#plt.close()
plt.show()
"""

"""
pn_file = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_mos1_grp.pha')[1].data
pn_chans = pn_file['CHANNEL']
pn_group = pn_file['GROUPING']

pn_chans_1 = pn_chans[pn_group==1]
for i in range(len(pn_chans)-1):
    if pn_chans_1[i] <= 2500:
        if i == 0:
            print(pn_chans_1[i],pn_chans_1[i]*5)
        else:
            print(pn_chans_1[i],pn_chans_1[i]*5,(pn_chans_1[i+1]-pn_chans_1[i])*5)
"""

"""
energies = np.array(list(np.arange(0.5,1.05,0.1)) + list(np.arange(1.5,12.5,0.5)))
fwhm = np.array([126.06,109.68,107.23,107.73,109.08,110.81,121.43,138.13,147.88,
                158.82,168.35,178.80,189.79,199.75,209.74,220.28,229.68,239.02,
                248.66,257.29,265.68,274.18,281.65,288.76,295.80,301.66,307.05,312.19])

for i in range(len(energies)):
    print(energies[i],fwhm[i]/3)

def sqrtfunc(x,a,b,c):
    return a*(x+b)**c


#popt,pcov = curve_fit(sqrtfunc,energies,fwhm,p0=[10,1.0,1.2],maxfev=50000)
#print(popt)
#print(np.sqrt(np.diag(pcov)))
#for i in range(len(energies)):
#    print(input_E[1:][i],sqrtfunc(energies,popt[0],popt[1],popt[2],popt[3])[i]/3)

plt.plot(energies,fwhm,'rx-')
#plt.plot(energies,sqrtfunc(energies,*popt),'bx-')
plt.annotate('48.2(E+2.5)^0.71',(4,300))
plt.xlabel('Energy (keV)',fontsize=12)
plt.ylabel('FWHM (eV)',fontsize=12)
plt.legend(('Measured','curve_fit'))
plt.show()
"""

"""
from scipy.interpolate import interp1d

##### FOR PN:
phases = np.arange(0.025,2,0.05)
pn_prof = np.array([1604.0,2164.0,1937.4058099752326,2676.0,3300.0,3839.637730914639,4038.0,4056.0,3540.7859071725998,2140.5, 3804.0,3668.75965331829,3826.5,3006.0,
1034.213403644357,1404.0,483.0,703.8586733495172,474.2092322833178,487.5]*2)
pn_inter = interp1d(phases,pn_prof,kind='cubic')
pn_prof_err = np.array([40.049968789001575,46.51881339845203,55.65052783180561,89.5991071384085,99.49874371066204,107.4471256562774,110.06361796706489,110.30865786510145,100.06852620303427,
56.66348030257233,75.53806987208498,74.28709253420669,75.76113779504647,67.14908785679812,39.400392759523484,45.89117562233503,26.916537667389523,32.50117373800949,
26.67635348746897,27.04163456597996]*2)

##### FOR MOS1:
mos1_prof = np.array([448.0,572.0,645.7875788096308,864.0,1152.0,1032.3053915630078,1167.0,1236.0,978.8970623790145,676.5,1014.0,1073.182141931018,1099.5,943.5,
444.50831287810894,507.0,318.0,210.0,208.69091755823553,178.58161485859577]*2)
mos1_inter = interp1d(phases,mos1_prof,kind='cubic')
mos1_prof_err = np.array([21.166010488516726,23.916521486202797,29.055356803825024,50.911688245431435,58.78775382679629,55.73925190095661,59.169248769948084,60.89334939055334,
48.883786528584615,31.855140872392965,38.99999999999997,40.190987067477344,40.610959112042714,37.61980861195333,25.8365208870395,27.57716446627533,21.840329667841537,
17.748239349298878, 17.700932599765082, 16.37054979377456]*2)

##### FOR MOS2:
mos2_prof = np.array([336.0,651.0,588.6128641386009,855.3569611778419,1308.0,1011.0,1218.0,1116.0,1323.0,804.5623109294818,1038.0,1129.4531069776747,1100.2580003973176,1054.5,
759.8367288864442,486.4282601354391,439.5,171.0,304.56707598060035,150.0]*2)
mos2_inter = interp1d(phases,mos2_prof,kind='cubic')
mos2_prof_err = np.array([18.33030277982336,25.514701644346147,24.9179694998494,44.77142958345511,62.64183905346333,55.07267925205741,60.44832503882968,57.86190456595775,63.0,36.57101413315826,
39.458839313897684,41.21432922222213,40.63895649552821,39.77122075068852,33.77886277372683,27.023792229746615,25.675864152935514,16.015617378046993,21.376418084077713,15.0]*2)

spacing = 0.05/20
phases_long = np.arange(0.025,1.975+spacing,spacing)

plt.figure()
#plt.errorbar(x=phases,y=pn_prof,yerr=pn_prof_err,color='r',drawstyle='steps-mid')
plt.plot(phases_long,pn_inter(phases_long),'r--')
#plt.errorbar(x=phases,y=mos1_prof,yerr=mos1_prof_err,color='b',drawstyle='steps-mid')
plt.plot(phases_long,mos1_inter(phases_long),'b--')
#plt.errorbar(x=phases,y=mos2_prof,yerr=mos2_prof_err,color='k',drawstyle='steps-mid')
plt.plot(phases_long,mos2_inter(phases_long),'k--')
plt.legend(('PN','MOS1','MOS2'),loc='best',fontsize=12)
plt.title('Exposure-corrected (using Stingray fold_events)',fontsize=12)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts',fontsize=12)
plt.show()
"""

#from stingray.pulse.pulsar import pulse_phase
#a = np.array([2,8,9,11,12,14,16,17,18,19,20])
#print(pulse_phase(a,1/3))
#b = np.array([2,8,9,11,12,14,16,17,18,19,20])-9
#print(pulse_phase(b,1/3))

#bgfile = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/combined_bg_pn.pha')
#ounts = bgfile[1].data['COUNTS']
#print(sum(counts))

"""
T_orb = 1/8.4712e-6

pn_expos_frac = np.array([0.05251764,0.07677285,0.14237649,0.13152822,0.04864319,
                        0.03305885,0.04322354,0.04902803,0.05,0.05,0.05,0.06141165,
                        0.1,0.09972594,0.1,0.08189309,0.03509542,0.06472928,
                        0.04661179,0.05677649])
0.00-0.05, 0.05-0.10, 0.10-0.15, 0.15-0.20, 0.20-0.25, 0.25-0.30, 0.30-0.35
0.35-0.40, 0.40-0.45, 0.45-0.50, 0.50-0.55, 0.55-0.60, 0.60-0.65, 0.65-0.70
0.70-0.75, 0.75-0.80, 0.80-0.85, 0.85-0.90, 0.90-0.95, 0.95-1.00
pn_expos = pn_expos_frac * T_orb
#pn_off_eclipse = pn_expos[4:-3]
#pn_on_eclipse = np.array(list(pn_expos[:4]) + list(pn_expos[-3:]))
#print(sum(pn_off_eclipse))
#print(sum(pn_on_eclipse))
pn_off = sum(pn_expos[5:-4])
pn_ingress = sum(pn_expos[-4:-2])
pn_on = sum(np.array(list(pn_expos[-2:]) + list(pn_expos[:2])))
pn_egress = sum(pn_expos[2:5])
print(pn_off,pn_ingress,pn_on,pn_egress)

mos1_expos_frac = np.array([0.07237615,0.09413837,0.15,0.13729413,0.06624048,0.03560002,
                            0.04915294,0.04981703,0.05,0.05,0.05169412,0.07558012,
                            0.1,0.09883275,0.1,0.09143779,0.04672711,0.07289424,
                            0.07510179,0.07400623])

mos1_expos = mos1_expos_frac * T_orb
#mos1_off_eclipse = mos1_expos[4:-3]
#mos1_on_eclipse = np.array(list(mos1_expos[:4]) + list(mos1_expos[-3:]))
#print(sum(mos1_off_eclipse))
#print(sum(mos1_on_eclipse))
mos1_off = sum(mos1_expos[5:-4])
mos1_ingress = sum(mos1_expos[-4:-2])
mos1_on = sum(np.array(list(mos1_expos[-2:]) + list(mos1_expos[:2])))
mos1_egress = sum(mos1_expos[2:5])
print(mos1_off,mos1_ingress,mos1_on,mos1_egress)

mos2_expos_frac = np.array([0.08557638,0.10630647,0.15,0.14237648,0.08618251,0.04152942,
                            0.05,0.04986786,0.05,0.05,0.05423529,0.08545808,0.1,
                            0.09966964,0.1,0.09576473,0.0643135,0.07628246,
                            0.09402377,0.07879993])

mos2_expos = mos2_expos_frac * T_orb
#mos2_off_eclipse = mos2_expos[4:-3]
#mos2_on_eclipse = np.array(list(mos2_expos[:4]) + list(mos2_expos[-3:]))
#print(sum(mos2_off_eclipse))
#print(sum(mos2_on_eclipse))
mos2_off = sum(mos2_expos[5:-4])
mos2_ingress = sum(mos2_expos[-4:-2])
mos2_on = sum(np.array(list(mos2_expos[-2:]) + list(mos2_expos[:2])))
mos2_egress = sum(mos2_expos[2:5])
print(mos2_off,mos2_ingress,mos2_on,mos2_egress)
"""

"""
pn_spectra = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/ngc300x1_pn_*pha'))
mos1_spectra = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/ngc300x1_mos1_*pha'))
mos2_spectra = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/ngc300x1_mos2_*pha'))

change_expos = open('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/change_expos.txt','w')
for i in range(len(pn_spectra)):
    change_expos.write('fparkey ' + str(pn_expos[i]) + ' ' + str(pn_spectra[i]) + ' EXPOSURE' + '\n')
for i in range(len(mos1_spectra)):
    change_expos.write('fparkey ' + str(mos1_expos[i]) + ' ' + str(mos1_spectra[i]) + ' EXPOSURE' + '\n')
for i in range(len(mos2_spectra)):
    change_expos.write('fparkey ' + str(mos2_expos[i]) + ' ' + str(mos2_spectra[i]) + ' EXPOSURE' + '\n')
change_expos.close()
"""

"""
cleaned_event = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_pn_clean.evt'
cleaned_spectra = '/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_pn_clean.pha'

plt.figure()
print(cleaned_event)
stdgti_indices = range(52,64)
stdgti_names = np.array(['stdgti' + str(i) for i in range(1,13)])
for i in range(len(stdgti_indices)):
    eventgti = fits.open(cleaned_event)[stdgti_indices[i]].data
    gti_sum = sum( [eventgti[j][1] - eventgti[j][0] for j in range(len(eventgti))] )
    for j in range(len(eventgti)):
        plt.plot(eventgti[j][0]-eventgti[0][0],i,'rx')
        plt.plot(eventgti[j][1]-eventgti[0][0],i,'rx')
    print(stdgti_names[i] + ' has sum of ' + str(gti_sum))

#plt.xscale('log')

plt.figure()
gti_indices = range(64,76)
gti_names = np.array(['gti0003','gti0103','gti0203','gti0303','gti0403','gti0503','gti0603','gti0703','gti0803','gti0903','gti1003','gti1103'])
for i in range(len(gti_indices)):
    eventgti = fits.open(cleaned_event)[gti_indices[i]].data
    gti_sum = sum( [eventgti[j][1] - eventgti[j][0] for j in range(len(eventgti))] )
    for j in range(len(eventgti)):
        plt.plot(eventgti[j][0]-eventgti[0][0],i,'rx')
        plt.plot(eventgti[j][1]-eventgti[0][0],i,'rx')
    print(gti_names[i] + ' has sum of ' + str(gti_sum))

#plt.xscale('log')


plt.show()

print(cleaned_spectra)
specgti_indices = [2] + list(range(4,15))
gti_names = np.array(['gti0003','gti0103','gti0203','gti0303','gti0403','gti0503','gti0603','gti0703','gti0803','gti0903','gti1003','gti1103'])
for i in range(len(gti_indices)):
    eventgti = fits.open(cleaned_spectra)[specgti_indices[i]].data
    gti_sum = sum( [eventgti[j][1] - eventgti[j][0] for j in range(len(eventgti))] )
    print(gti_names[i] + ' has sum of ' + str(gti_sum))
"""

"""
eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000023/merged000023_nicersoft_bary_phase.evt'
gtis = fits.open(eventfile)[2].data
print( sum([gtis[i][1]-gtis[i][0] for i in range(len(gtis))]) )
"""

"""
results_file = '/Volumes/Samsung_T5/NICER-data/at2019wey/from_dragonfly/last_line.txt'
contents = open(results_file,'r').read().split('\n')
freq_sig = np.array([contents[i] for i in range(1,len(contents),3)])
freqs = np.array([np.float(freq_sig[i].split()[0]) for i in range(len(freq_sig))])
plt.hist(freqs,bins=700)
plt.xlabel('Frequency (Hz)',fontsize=12)
plt.plot()
plt.show()
"""

"""
obs1_pn = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_pn_bg.pha')[9].data
print( sum([obs1_pn[i][1]-obs1_pn[i][0] for i in range(len(obs1_pn))]) )

obs2_pn = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_pn_bg.pha')[9].data
print( sum([obs2_pn[i][1]-obs2_pn[i][0] for i in range(len(obs2_pn))]) )

obs1_mos1 = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_mos1_bg.pha')[2].data
print( sum([obs1_mos1[i][1]-obs1_mos1[i][0] for i in range(len(obs1_mos1))]) )

obs2_mos1 = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_mos1_bg.pha')[2].data
print( sum([obs2_mos1[i][1]-obs2_mos1[i][0] for i in range(len(obs2_mos1))]) )

obs1_mos2 = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010101/PROC/xmm_0791010101_mos2_bg.pha')[2].data
print( sum([obs1_mos2[i][1]-obs1_mos2[i][0] for i in range(len(obs1_mos2))]) )

obs2_mos2 = fits.open('/Volumes/Samsung_T5/NGC300_XMMdata/0791010301/PROC/xmm_0791010301_mos2_bg.pha')[2].data
print( sum([obs2_mos2[i][1]-obs2_mos2[i][0] for i in range(len(obs2_mos2))]) )
"""

"""
folder = '/Volumes/Samsung_T5/n300_ulx_2020/spectra_10d/'
fparkey = open(folder + 'fparkey_jsgrp_10d.txt','w')
jsgrp_files = sorted(glob.glob(folder+'jsgrp_*_cl_*pha'))
for i in range(len(jsgrp_files)):
    bg_file = str(pathlib.Path(jsgrp_files[i]).name)[2:16] + 'bg_cl50.pha'
    fparkey.write('fparkey ' + bg_file + ' ' + jsgrp_files[i] + ' BACKFILE' + '\n')
fparkey.close()
"""

"""
folder = '/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/'
load_data = open(folder+'bgsub_load_all.xcm','w')
jsgrp_files = sorted(glob.glob(folder+'grp_*_bgsub_*pha'))
load_data.write('data ')
for i in range(len(jsgrp_files)):
    load_data.write(str(i+1)+':'+str(i+1) + ' ' + str(pathlib.Path(jsgrp_files[i]).name) + ' ')
load_data.write('\n')
load_data.write('ignore **:0.0-0.295,12.01-**' + '\n')
load_data.write('setplot energy' + '\n')
load_data.write('setplot background' + '\n')
load_data.close()
"""


"""
spec = '/Volumes/Samsung_T5/NGC300_XMMdata/jsgrp_ngc300x1_mos2_020-085.pha'
counts = sum(fits.open(spec)[1].data['COUNTS'])
print(counts)
"""

"""
eventfiles = sorted(glob.glob('/Volumes/Samsung_T5/NICER-data/at2019wey/accelsearch_GTIs/2008041401_filt_bary_*E0050-0200.evt'))
exps = []
for i in range(len(eventfiles)):
    exps.append(fits.open(eventfiles[i])[1].header['EXPOSURE'])
exps = np.array(exps)
print(min(exps),max(exps),np.median(exps),np.std(exps))

print(len(exps))
print(len(exps[exps>64]))
print(min(exps[exps>64]),max(exps[exps>64]),np.median(exps[exps>64]),np.std(exps[exps>64]))

plt.hist(exps,bins=100)
plt.show()
"""

"""
eventfile = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000026/merged000026_nicersoft_bary_center_phase.evt'
times = fits.open(eventfile)[1].data['TIME'] #getting array of times
gtis_data = fits.open(eventfile)[2].data #getting GTIs
T = sum([ gtis_data[i]['STOP']-gtis_data[i]['START'] for i in range(len(gtis_data)) ]) #exposure time
T0_MJD = fits.open(eventfile)[1].header['MJDREFI'] + fits.open(eventfile)[1].header['MJDREFF'] + fits.open(eventfile)[1].header['TSTART']/86400 #SWIFT
gtis_conform = []
for i in range(len(gtis_data)):
    gtis_conform.append([gtis_data[i][0],gtis_data[i][1]]) #conform to the input that Stingray uses

nbins = 20
freq = 1/(0.38196674159633911059*86400)
freqdot = 0#-2.3091587842511194363e-11/(0.38196674159633911059*86400)**2
freqdotdot = 0

expos = Lv2_phase.phase_exposure(times[0]-times[0],times[-1]-times[0],1/freq,nbin=nbins,gtis=np.array(gtis_conform)-times[0])
total_expos = np.array(list(expos) + list(expos))

print('Original expos:')
print(expos)

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

"""
mergedfile = '/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000026/merged000026_nicersoft_bary_center_phase.evt'
orbphase = fits.open(mergedfile)[1].data['ORBIT_PHASE']
exposure = fits.open(mergedfile)[1].header['EXPOSURE']
plt.hist(orbphase,bins=20)
plt.xlabel('Phase',fontsize=12)
plt.ylabel('Counts',fontsize=12)
plt.show()
"""

"""ad_all
binsize = '10d'
sourcefile = open('/Volumes/Samsung_T5/n300_ulx_2020/spectra_' + binsize + '/3C50x_BACKFILE.go','w')
jsgrp = sorted(glob.glob('/Volumes/Samsung_T5/n300_ulx_2020/spectra_' + binsize + '/jsgrp_58*3C50x*pha'))
for i in range(len(jsgrp)):
    mjd = jsgrp[i][-24:-19]
    backfile = 'grp_' + mjd + '_' + binsize + '_x_cl50.pha'
    sourcefile.write('fparkey ' + backfile + ' ' + str(jsgrp[i]) + '[1] BACKFILE' + '\n')
sourcefile.close()
"""

"""
cl_files = sorted(glob.glob(Lv0_dirs.NGC300_2020 + 'cl50/pha/cl50_*pha'))
xbg_files = sorted(glob.glob(Lv0_dirs.NGC300_2020 + 'bg_cl50/xpha/xbg_cl50_*pha'))
xbgsub_command = open(Lv0_dirs.NGC300_2020 + 'xbgsub_generate.go','w')
for i in range(len(xbg_files)): #because there are 387 xbg files and 396 cl files
    gtino = xbg_files[i][-8:-4]
    xbgfile = xbg_files[i]
    clfile = Lv0_dirs.NGC300_2020 + 'cl50/pha/cl50_' + gtino + '.pha'
    outfile = Lv0_dirs.NGC300_2020 + 'bgsub_cl50/xpha/xbgsub_cl50_' + gtino + '.pha'

    cl_exp = fits.open(clfile)[1].header['EXPOSURE']
    xbg_exp = fits.open(xbgfile)[1].header['EXPOSURE']

    xbgsub_command.write('mathpha ' + '"(' + "'" + clfile + "'-'" + xbgfile + "')" + '"' + " R outfil='" + outfile + "' exposure=" + str(cl_exp+xbg_exp) + ' errmeth=gaussian properr=yes ncomments=0 areascal=NULL clobber=YES')
    xbgsub_command.write('\n')
xbgsub_command.close()
"""

"""
bgsub_soft, bgsub_soft_err, bgsub_softintens, bgsub_softintens_err = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'soft_intens_bgsub.txt',dtype='float64',usecols=(0,1,2,3),unpack=True)
xbgsub_soft, xbgsub_soft_err, xbgsub_softintens, xbgsub_softintens_err = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'soft_intens_xbgsub.txt',dtype='float64',usecols=(0,1,2,3),unpack=True)

bgsub_mjd, bgsub_intens, bgsub_intens_err = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'time_intens_bgsub.txt',dtype='float64',usecols=(0,1,2),unpack=True)
xbgsub_mjd, xbgsub_intens, xbgsub_intens_err = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'time_intens_xbgsub.txt',dtype='float64',usecols=(0,1,2),unpack=True)

plt.figure(1)
plt.errorbar(x=bgsub_soft,y=bgsub_softintens,xerr=bgsub_soft_err,yerr=bgsub_softintens_err,color='r',fmt='x',label='bgsub')
plt.errorbar(x=xbgsub_soft,y=xbgsub_softintens,xerr=xbgsub_soft_err,yerr=xbgsub_softintens_err,color='b',fmt='x',label='xbgsub')
plt.legend(('bgsub','xbgsub'),loc='best',fontsize=12)
plt.axhline(y=0,lw=0.5,alpha=0.5)
plt.xlabel('Soft Color: (1-2 keV)/(0.4-1 keV)',fontsize=12)
plt.ylabel('Intensity in counts/s (0.4-12 keV)',fontsize=12)

plt.figure(2)
plt.errorbar(x=bgsub_mjd,y=bgsub_intens,yerr=bgsub_intens_err,color='r',fmt='x',label='bgsub')
plt.errorbar(x=xbgsub_mjd,y=xbgsub_intens,yerr=xbgsub_intens_err,color='b',fmt='x',label='xbgsub')
plt.legend(('bgsub','xbgsub'),loc='best',fontsize=12)
plt.axhline(y=0,lw=0.5,alpha=0.5)
plt.xlabel('Time (MJD)',fontsize=12)
plt.ylabel('Intensity in counts/s (0.4-12 keV)',fontsize=12)

plt.show()
"""

"""
plt.figure()

init_back_array = sorted(glob.glob(Lv0_dirs.NGC300_2020 + 'bg_cl50/pha/*.pha'))
Porb = (1/8.4712e-6)/86400
T0 = 58239.3498
for i in range(len(init_back_array)):
    init_back = init_back_array[i]
    gti_no = init_back[-8:-4]
    init_cl = Lv0_dirs.NGC300_2020 + 'cl50/pha/cl50_' + gti_no + '.pha'

    backfile = str(pathlib.Path(init_back).name)
    comb_spec_folder = Lv0_dirs.NGC300_2020 + 'bg_cl50/xpha/'

    tstart = fits.open(init_cl)[1].header['TSTART']
    tstop = fits.open(init_cl)[1].header['TSTOP']
    centroid_t = fits.open(init_cl)[1].header['MJDREFI'] + fits.open(init_cl)[1].header['MJDREFF'] + ((tstart+tstop)/2)/86400
    exp_time = fits.open(init_cl)[1].header['EXPOSURE']

#first assuming the timing model, determine how much time (in d) has passed
    orb_phase = (centroid_t-T0)%Porb/Porb
    plt.errorbar(x=orb_phase,xerr=exp_time/(Porb*86400),y=0.0012 + (0.0003-0.0012)/380*i)

a,b,c = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'swift_shifted_folded_curve.txt',usecols=(0,1,2),unpack=True)
plt.errorbar(x=a,y=b,yerr=c,color='r',drawstyle='steps-mid')
plt.xlabel('Phase')
plt.ylabel('Rate (Counts/s)')
plt.axvline(x=0.2,lw=0.5,alpha=0.5,color='c')
plt.axvline(x=0.85,lw=0.5,alpha=0.5,color='c')
plt.annotate('off-eclipse',(0.4,0.0006))

plt.show()
"""

"""
orig_file = fits.open('/Volumes/Samsung_T5/NICER-data/at2019wey/2008041401_filt_bary.evt')[1].data['TIME']
first_file = fits.open('/Volumes/Samsung_T5/NICERsoft_outputs/merged_events/merged000027/merged000027_nicersoft_bary.evt')[1].data['TIME']
mjdrefi = fits.open('/Volumes/Samsung_T5/NICER-data/at2019wey/2008041401_filt_bary.evt')[1].header['MJDREFI']
mjdreff = fits.open('/Volumes/Samsung_T5/NICER-data/at2019wey/2008041401_filt_bary.evt')[1].header['MJDREFF']

print(mjdrefi + mjdreff + orig_file[-1]/86400)
print(mjdrefi + mjdreff + first_file[0]/86400)
"""

"""
gtis = sorted(glob.glob('/Volumes/Samsung_T5/NICER-data/at2019wey/gtis/00064s_*gti'))
putting = []
for i in range(len(gtis)):
    gtiname = str(pathlib.Path(gtis[i]).name)
    if float(gtiname[-10:-4]) > 73450:
        putting.append(gtiname)

print(len(putting))
print(putting[0],putting[-1])
"""

"""
jsgrp = sorted(glob.glob('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_58*_3C50_*pha'))
for i in range(len(jsgrp)):
    input = jsgrp[i]
    output = jsgrp[i][:-13] + '3C50xs_cl50.pha'
    subprocess.run(['cp',input,output])
"""

"""
fparkey = open(Lv0_dirs.NGC300_2020 + 'spectra_05d/3C50xs_exposure.go','w')
spectra = sorted(glob.glob(Lv0_dirs.NGC300_2020 + 'spectra_05d/jsgrp_*_05d_3C50xs_cl50.pha'))
for i in range(len(spectra)):
    mjd = str(pathlib.Path(spectra[i]).name)[-25:-20]
    bg_file = Lv0_dirs.NGC300_2020 + 'spectra_05d/grp_' + mjd + '_05d_xsbg_cl50.pha'
    exposure = fits.open(spectra[i])[1].header['EXPOSURE']
    fparkey.write('fparkey ' + str(exposure) + ' ' + bg_file + ' EXPOSURE' + '\n')
fparkey.close()
"""

"""
period = np.array([1.37,0.42,1.1,31.6])
pdot = np.array([-2e-10,-3.5e-11,-8e-10,-5.56e-7])
freq = 1/period #order of M82 X-2, NGC 7793 P13, NGC 5907 ULX 1, NGC 300 ULX-1
fdot = -pdot/period**2
print(freq)
print(fdot)

plt.plot(freq,fdot,'x')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Frequency (Hz)',fontsize=12)
plt.ylabel('Frequency Derivative (Hz/s)',fontsize=12)
plt.xlim([0,1000])
plt.ylim([1e-18,1e-7])
plt.show()
"""

"""
xs_mjd,xs_s1,xs_s2,xs_a,xs_b,xs_c,xs_d,xs_in = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.xsbgsub_cl50_g2020norm_05d.fffphot',usecols=(0,1,2,3,4,5,6,7),unpack=True)
err_mjd,err_s1,err_s2,err_a,err_b,err_c,err_d,err_in = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.xbgsub_cl50_g2020err_norm_05d.fffphot',usecols=(0,1,2,3,4,5,6,7),unpack=True)
for i in range(len(xs_mjd)):
    print(xs_mjd[i],round(xs_s1[i]+err_s1[i],4),round(xs_s2[i]+err_s2[i],4),round(xs_a[i]+err_a[i],4),round(xs_b[i]+err_b[i],4),round(xs_c[i]+err_c[i],4),round(xs_d[i]+err_d[i],4),round(xs_in[i]+err_in[i],4))
"""

"""
eventfile = Lv0_dirs.NICER_DATADIR + 'at2019wey/final.evt'
fitsfile = fits.open(eventfile)[1]
print(fitsfile.header['MJDREFI'] + fitsfile.header['MJDREFF'] + fitsfile.data['TIME'][-1]/86400)
"""

"""
gtis = sorted(glob.glob(Lv0_dirs.NICER_DATADIR + 'at2019wey/accelsearch_GTIs/GTI*gti'))
gti_len = []
for i in range(len(gtis)):
    gti = fits.open(gtis[i])[1].data
    gti_len.append(gti['STOP'][0]-gti['START'][0])

print(np.min(gti_len),np.max(gti_len))
gti_len = np.array(gti_len)
print(sum(gti_len))
filtered_gtis = gti_len[gti_len>64]
print(len(gti_len))
print(len(filtered_gtis))
print(np.median(filtered_gtis))
"""

"""
from pint.eventstats import sf_z2m,z2m,sig2sigma
z_vals = 144.37
probs = sf_z2m(z_vals)
print(probs)
significances = sig2sigma(probs)
print(significances)
"""

mjd = np.array([59149.0232,59149.0877,59149.2968,59149.6840,59149.7330,59150.3138,59150.3784,59150.7193])
freqs = np.array([376.043646,376.056458,376.056901,376.045625,376.050446,376.041895,376.045996,376.05253])
plt.plot(mjd,freqs,'rx-')
plt.show()

timeend = time.time()

print(str(timeend-timestart) + ' seconds')
