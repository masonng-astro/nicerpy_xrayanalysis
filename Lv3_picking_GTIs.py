#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues May 14 2:30pm 2019

End goal:

Probably need to do 2 files for each segment:
1. original .evt file which has the proper timestamp information etc. with the appropriate GTI card
2. refined .evt file that has a GTI card with just (0,1000) for the GTI!

Overarching steps:

1) Use Python to obtain the timestamps that corresponding to T=0, 1000s, 2000s, of data
2) Use niextract-events to get *.evt
3) Also on Python - make a new copy of the file first, then EDIT that copy such that the GTI card has just (0,1000)

"""
from __future__ import division, print_function
import numpy as np
from scipy import stats, signal
from tqdm import tqdm
import matplotlib.pyplot as plt
import subprocess
from astropy.io import fits
from astropy.io.fits import getheader,getdata,update
import Lv0_dirs,Lv2_ps_method,Lv2_phase

Lv0_dirs.global_par()

def get_gtis(obsid):
    """
    Get a data dictionary from the NICERsoft pipe folders for a desired ObsID -
    data was pre-processed by NICERsoft!
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data

    return gtis

def change_GTI_card(obsid):
    """
    Change the GTI card such that the duration of the observation is the total
    length of the GTIs! So NOT timestamps.

    5/14/19: Do it manually for now - so print out start time (to check) and the end
    time (and manually enter the number with "fv"). Remember to MAKE A COPY!
    5/15/19:
    """
    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    event = fits.open(event)
    gtis = event[2].data

    T = sum([gtis[i][1]-gtis[i][0] for i in range(len(gtis))])

    startt = gtis[0][0]
    endt = gtis[0][0] + T

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    newevent = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary_short.evt'
    subprocess.check_call(['cp',event,newevent])
    new_list = np.rec.array([(startt,endt)])
    data,hdr = getdata(event,2,header=True)
    update(newevent,new_list,2,header=hdr)

    return startt,endt

#print(change_GTI_card('1034090111'))

def get_timestamps(obsid,desired_length):
    gtis = get_gtis(obsid)
    T = sum([gtis[i][1]-gtis[i][0] for i in range(len(gtis))])
    gti_1000s = [gtis[0][0]]

    deltat = 0
    boundary = gtis[0][0] #initialize boundary
    for i in range(len(gtis)-1):
        if deltat + (gtis[i][1]-boundary) < desired_length:
            deltat += (gtis[i][1]-boundary)
            boundary = gtis[i+1][0]
        else:
            boundary = (desired_length-deltat) + gtis[i][0]
            if i == len(gtis)-2:
                gti_1000s += [boundary]
            else:
                gti_1000s += [boundary,boundary]
            deltat = 0

    gti_1000s_tuple = []
    for i in range(0,len(gti_1000s)-1,2):
        gti_1000s_tuple.append((gti_1000s[i],gti_1000s[i+1]))

    return gti_1000s_tuple

#segment_bounds = np.array(get_timestamps('1034090111',1000))

def change_data_stitched(obsid):
    """
    Change the GTI card such that the duration of the observation is the total
    length of the GTIs! Alter timestamps by an offset too!

    This is for the stitched files!
    """
    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    newevent = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary_stitched.evt'

    data1,hdr1 = getdata(event,1,header=True)
    data2,hdr2 = getdata(event,2,header=True)

    old_times = data1['TIME']
    print(len(old_times))
    gtis = get_gtis(obsid)
    T = gtis[-1][1] - gtis[0][0]

    new_times = list(old_times[(old_times>=gtis[0][0])&(old_times<=gtis[0][1])])
    for i in range(1,len(gtis)):
        segment = np.array(old_times[(old_times>=gtis[i][0])&(old_times<=gtis[i][1])])
        segment_offset = segment - (gtis[i][0]-gtis[i-1][1])
        segment_offset = list(segment_offset)
        new_times += segment_offset

    print(len(new_times))
    new_times = new_times[::1000]
    print(len(new_times))
    subprocess.check_call(['cp',event,newevent])
    new_gti = np.rec.array([(0,T)],names='START,STOP',formats='D,D')
    new_data = np.rec.array([new_times],names='TIME')

    update(newevent,new_gti,2,header=hdr2)
    update(newevent,new_data,1,header=hdr1)

    return


def change_data_stitched(obsid):
    """
    Change the GTI card such that the duration of the observation is the total
    length of the GTIs! Alter timestamps by an offset too!

    This is for the stitched files!
    Columns: TIME, RAWX, RAWY, PHA, PHA_FAST, DET_ID, DEADTIME, EVENT_FLAGS,
    TICK, MPU_A_TEMP, MPU_UNDER_COUNT, PI_FAST, PI, PI_RATIO
    """
    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    newevent = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary_stitched.evt'

    prihdr = getheader(event,0)
    data1,hdr1 = getdata(event,1,header=True)
    data2,hdr2 = getdata(event,2,header=True)

    old_times = data1['TIME']
    rawx = data1['RAWX']
    rawy = data1['RAWY']
    pha = data1['PHA']
    pha_fast = data1['PHA_FAST']
    det_id = data1['DET_ID']
    deadtime = data1['DEADTIME']
    event_flags = data1['EVENT_FLAGS']
    tick = data1['TICK']
    mpu_a_temp = data1['MPU_A_TEMP']
    mpu_under_count = data1['MPU_UNDER_COUNT']
    pi_fast = data1['PI_FAST']
    pi = data1['PI']
    pi_ratio = data1['PI_RATIO']

    gtis = get_gtis(obsid)
    T = gtis[-1][1] - gtis[0][0]

    new_times = list(old_times[(old_times>=gtis[0][0])&(old_times<=gtis[0][1])])
    for i in range(1,len(gtis)):
        segment = np.array(old_times[(old_times>=gtis[i][0])&(old_times<=gtis[i][1])])
        segment_offset = segment - (gtis[i][0]-gtis[i-1][1])
        segment_offset = list(segment_offset)
        new_times += segment_offset

    prihdu = fits.PrimaryHDU(header=prihdr)

    tablehdu1 = fits.BinTableHDU.from_columns([
            fits.Column(name='TIME',format='1D',unit='s',array=new_times),
            fits.Column(name='RAWX',format='1B',unit='pixel',array=rawx),
            fits.Column(name='RAWY',format='1B',unit='pixel',array=rawy),
            fits.Column(name='PHA',format='1I',unit='chan',array=pha),
            fits.Column(name='PHA_FAST',format='1I',unit='chan',array=pha_fast),
            fits.Column(name='DET_ID',format='1B',unit='',array=det_id),
            fits.Column(name='DEADTIME',format='1B',unit='s',array=deadtime),
            fits.Column(name='EVENT_FLAGS',format='8X',unit='',array=event_flags),
            fits.Column(name='TICK',format='1K',unit='',array=tick),
            fits.Column(name='MPU_A_TEMP',format='I',unit='Celsius',array=mpu_a_temp),
            fits.Column(name='MPU_UNDER_COUNT',format='J',unit='',array=mpu_under_count),
            fits.Column(name='PI_FAST',format='1I',unit='CHAN',array=pi_fast),
            fits.Column(name='PI',format='1I',unit='CHAN',array=pi),
            fits.Column(name='PI_RATIO',format='1E',unit='',array=pi_ratio)])

    tablehdu2 = fits.BinTableHDU.from_columns([
            fits.Column(name='START',format='1D',unit='s',array=np.array([0])),
            fits.Column(name='STOP',format='1D',unit='s',array=np.array([T])) ])

    new_hdu = fits.HDUList([prihdu,tablehdu1,tablehdu2])
    new_hdu.writeto(newevent,overwrite=True)

    return

change_data_stitched('1034090111')

def change_GTI_stitched(obsid,gtino,desired_length):
    """
    Change the GTI card such that the duration of the observation is the total
    length of the GTIs! Alter timestamps by an offset too!

    This is for the stitched files!
    """
#    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary_GTI' + str(gtino) + '_' +str(desired_length) + 's_stitched.evt'
#    newevent = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary_GTI' + str(gtino) + '_' +str(desired_length) + 's_stitched_short.evt'

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary_GTI' + str(gtino) + '_' +str(desired_length) + 's_stitched.evt'
    newevent = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary_GTI' + str(gtino) + '_' +str(desired_length) + 's_stitched_short.evt'

    data1,hdr1 = getdata(event,1,header=True)
    data2,hdr2 = getdata(event,2,header=True)

    old_times = data1['TIME']
    gtis = get_gtis(obsid)

    new_times = old_times[(old_times>=gtis[0][0])&(old_times<=gtis[0][1])]
    for i in range(1,len(gtis)):
        segment = np.array(old_times[(old_times>=gtis[i][0])&(old_times<=gtis[i][1])])
        segment_offset = segment - (gtis[i][0]-gtis[i-1][1])
        segment_offset = list(segment_offset)
        new_times += segment_offset

    data1['TIME'] = new_times

    subprocess.check_call(['cp',event,newevent])
    new_gti = np.rec.array([(0,desired_length)],names='START,STOP',units='s,s')
    new_data = np.rec.array([new_times],names='TIME',units='s,s')

    update(newevent,new_gti,2,header=hdr2)
    update(newevent,new_data,1,header=hdr1)

    return

"""
obsid = '1034090111'
event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
from astropy.io.fits import getdata,update

data1,hdr1 = getdata(event,1,header=True)
data2,hdr2 = getdata(event,2,header=True)
"""

#for i in range(1,29):
#    change_GTI_stitched('1034090111',i,1000)
################################## TEST #######################################
"""
event = Lv0_dirs.NICERSOFT_DATADIR + '1034090111_pipe/ni1034090111_nicersoft_bary.evt'
event = fits.open(event)
times = event[1].data['TIME']
counts = np.ones(len(times))

shifted_t = times-times[0]
tbin_size = 0.1
t_bins = np.linspace(shifted_t[0],shifted_t[-1],(shifted_t[-1]-shifted_t[0])*1/tbin_size+1)

summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins)

plt.plot(t_bins[:-1],summed_data,'rx')
plt.show()
"""
