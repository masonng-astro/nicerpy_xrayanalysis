#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thurs Jun 11 7:24pm 2020

Creating light curves and periodograms from Swift data

"""
from __future__ import division, print_function
import numpy as np
import numpy.ma as ma
from astropy.io import fits
import matplotlib.pyplot as plt
import Lv0_dirs,Lv1_data_gtis,Lv2_presto_subroutines,Lv2_mkdir
from matplotlib.backends.backend_pdf import PdfPages
import os
from tqdm import tqdm
import subprocess
import pathlib
import glob

Lv0_dirs.global_par() #obtaining the global parameters

def get_ra_dec(eventfile):
    """
    Obtain the RA_OBJ and DEC_OBJ corresponding to the observation!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    """
    event = fits.open(eventfile)
    event_header = event[1].header

    return event_header['RA_OBJ'], event_header['DEC_OBJ']

def time_order(eventlist,mjd1,mjd2):
    """
    Takes as input, a list of event files, and outputs a subset of event files,
    which are ordered in time, and are contained within mjd1 and mj2

    eventlist - list of event files
    mjd1 - lower bound on MJD (i.e., earliest)
    mjd2 - upper bound on MJD (i.e., latest)
    """
    if type(eventlist) != np.ndarray and type(eventlist) != list:
        raise TypeError('eventfile should be an array or a list!')

    start_times = [fits.open(eventlist[i])[1].header['TSTART'] for i in range(len(eventlist))]
    time_ordered = np.argsort(start_times)
    ordered_eventlist = np.array(eventlist)[time_ordered]

    if mjd1 != 0 and mjd2 != 0:
        truncated_eventlist = []
        for i in range(len(ordered_eventlist)):
            obs_time = fits.open(ordered_eventlist[i])[1].header['MJD-OBS']
            if obs_time >= mjd1 and obs_time <= mjd2:
                truncated_eventlist.append(ordered_eventlist[i])

        return truncated_eventlist

    else:
        return ordered_eventlist

def att_file_use(eventfile):
    """
    For a given event file, determines what attitude file to use

    eventfile - path to the event file
    """
    obsid = str(pathlib.Path(eventfile).name)[:13]

    attflag = fits.open(eventfile)[1].header['ATTFLAG']
    if attflag == '110':
        return '/Volumes/Samsung_T5/NGC300_ULX_Swift/auxil/' + obsid + 'pat.fits.gz'
    elif attflag == '100':
        return '/Volumes/Samsung_T5/NGC300_ULX_Swift/auxil/' + obsid + 'sat.fits.gz'
    elif attflag == '111':
        return '/Volumes/Samsung_T5/NGC300_ULX_Swift/auxil/' + obsid + 'uat.fits.gz'

def barycorr(eventfile,outfile,refframe,orbit_file,output_folder):
    """
    General function to perform the barycenter corrections for a Swift event file

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    outfile - path to the output event file with barycenter corrections applied
    refframe - reference frame for barycenter corrections (usually ICRS)
    orbit_file - path to the orbit file of the observation
    output_folder - path to the folder where the outfile will be
    """
    obsid = eventfile[2:13]

    logfile = output_folder + 'barycorr_notes.txt'
    ra,dec = get_ra_dec(eventfile)

    with open(logfile,'w') as logtextfile:
        output = subprocess.run(['barycorr',eventfile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(ra),'dec='+str(dec),'refframe='+str(refframe),'clobber=YES'],capture_output=True,text=True)
        logtextfile.write(output.stdout)
        logtextfile.write('*------------------------------* \n')
        logtextfile.write(output.stderr)
        logtextfile.close()

    return

def xselect_script(eventlist,regfile,binsize,mjd1,mjd2):
    """
    Writes a script so that XSELECT can take in the events, bins them, and
    outputs them into a light curve (.lc) form

    eventlist - a list of event files
    binsize - desired bin size for the light curve
    """
    parent_folder = str(pathlib.Path(eventlist[0]).parent)
    script_name = parent_folder + '/xselect_earlier_ulx1_instructions.txt'

    writing = open(script_name,'w')
    writing.write('set mission swift' + '\n')
    writing.write('set inst xrt' + '\n')
    for i in range(len(eventlist)):
        if (fits.open(eventlist[i])[1].header['MJD-OBS'] >= mjd1) and (fits.open(eventlist[i])[1].header['MJD-OBS'] <= mjd2):
            eventfile = str(pathlib.Path(eventlist[i]).name)
            obsid = eventfile[:13]
            writing.write('set datadir ' + parent_folder + '\n')
            writing.write('read events ' + eventfile + '\n')
            writing.write('filter region ' + regfile + '\n')
            writing.write('extract curve binsize=' + str(binsize) + '\n')
            writing.write('save curve ' + obsid + '_ulx1.lc' + '\n')
            writing.write('clear data' + '\n')
            writing.write('clear region' + '\n')
            writing.write(' ' + '\n')
    writing.close()

def xrtlccorr(event_files,lc_files,att_files,hk_files,mjd1,mjd2):
    """
    Runs xrtlccorr to produce the corrected light curves

    event_files - list of the barycenter-corrected event files
    lc_files - list of the un-corrected .lc files
    att_files - list of barycenter-corrected att_files
    hk_files - list of barycenter-corrected hk_files
    """
    parent_folder = str(pathlib.Path(lc_files[0]).parent)
    instruct_file = parent_folder + '/xrtlccorr_earlier_ulx1_instruc.txt'

    ordered_events = time_order(event_files,0,0)
    ordered_lc = time_order(lc_files,0,0)
    ordered_att = time_order(att_files,0,0)
    ordered_hk = time_order(hk_files,0,0)

    writing = open(instruct_file,'w')
    for i in range(len(ordered_lc)):
        if (fits.open(ordered_lc[i])[1].header['MJD-OBS']>=mjd1) and (fits.open(ordered_lc[i])[1].header['MJD-OBS']<=mjd2):
            outfile = ordered_lc[i][:-3] + '_corr.lc'
            output = 'xrtlccorr'+ ' lcfile='+ordered_lc[i] + ' regionfile=NONE' + ' outfile='+outfile + ' corrfile=DEFAULT' + ' attfile='+ordered_att[i] + ' hdfile='+ordered_hk[i] + ' outinstrfile=DEFAULT' + ' pcnframe=0' + ' psfflag=yes' + ' energy=1.0' + ' createinstrmap=yes' + ' infile='+ordered_events[i] + ' clobber=yes'
            writing.write(output + '\n')
    writing.close()

def lcmath(corr_lc_files,corr_bg_files,bg_scale):
    """
    Running lcmath to do background subtraction on the xrtlccorr-corrected light curves

    corr_lc_files - list of corrected sw*_corr.lc files
    corr_bg_files - list of corrected sw*_bg_corr.lc files
    bg_scale - scaling factor for background
    """
    parent_folder = str(pathlib.Path(corr_lc_files[0]).parent)
    lcmath = open(parent_folder + '/lcmath_instruct.txt','w')

    for i in range(len(corr_lc_files)):
        inputfile = corr_lc_files[i]
        bgfile = corr_bg_files[i]
        outputfile = corr_lc_files[i][:-7] + 'bgsub_corr.lc'
        lcmath.write('lcmath ' + 'infile='+str(inputfile) + ' bgfile=' + bgfile + ' outfile=' + outputfile + ' multi=1' + ' multb=' + str(bg_scale) + ' addsubr=no' + ' err_mode=2' + '\n')

    lcmath.close()

def get_bgsub(corr_lc_files,corr_bg_files,bg_scale):
    """
    Getting background-subtracted light curves in the case of low counts (i.e.,
    cannot use lcmath!)

    corr_lc_files - list of corrected sw*_corr.lc files
    corr_bg_files - list of corrected sw*_bg_corr.lc files
    bg_scale - scaling factor for background
    """
    lc_start_times = [fits.open(corr_lc_files[i])[1].header['TIMEZERO'] for i in range(len(corr_lc_files))]
    lc_time_ordered = np.argsort(lc_start_times)
    lc_ordered_eventlist = np.array(corr_lc_files)[lc_time_ordered]

    lc_times = []
    lc_rates = []
    lc_errors = []
    lc_fracexp = []
    for i in tqdm(range(len(lc_ordered_eventlist))): #getting every corrected light curve FITS file
        data_arrays = fits.open(lc_ordered_eventlist[i])[1].data
        tstart = fits.open(lc_ordered_eventlist[i])[1].header['TIMEZERO'] #get TIMEZERO
        lc_times += list(tstart+data_arrays['TIME']) #add to an overall list, the 'TIME' array from the corrected light curve FITS file, plus timezero
        lc_rates += list(data_arrays['RATE']) #get the corresponding column for the rates
        lc_errors += list(data_arrays['ERROR']) #get the corresponding column for the errors
        lc_fracexp += list(data_arrays['FRACEXP']) #get the corresponding column for the fractional exposure

    lc_times = np.array(lc_times)
    lc_trunc_times = lc_times-lc_times[0]
    lc_rates = np.array(lc_rates)
    lc_errors = np.array(lc_errors)
    lc_fracexp = np.array(lc_fracexp)

    ##########

    bg_start_times = [fits.open(corr_bg_files[i])[1].header['TIMEZERO'] for i in range(len(corr_bg_files))]
    bg_time_ordered = np.argsort(bg_start_times)
    bg_ordered_eventlist = np.array(corr_bg_files)[bg_time_ordered]

    bg_times = []
    bg_rates = []
    bg_errors = []
    bg_fracexp = []
    for i in tqdm(range(len(bg_ordered_eventlist))): #getting every corrected light curve FITS file
        data_arrays = fits.open(bg_ordered_eventlist[i])[1].data
        tstart = fits.open(bg_ordered_eventlist[i])[1].header['TIMEZERO'] #get TIMEZERO
        bg_times += list(tstart+data_arrays['TIME']) #add to an overall list, the 'TIME' array from the corrected light curve FITS file, plus timezero
        bg_rates += list(data_arrays['RATE']) #get the corresponding column for the rates
        bg_errors += list(data_arrays['ERROR']) #get the corresponding column for the errors
        bg_fracexp += list(data_arrays['FRACEXP']) #get the corresponding column for the fractional exposure

    bg_times = np.array(bg_times)
    bg_trunc_times = bg_times-bg_times[0]
    bg_rates = np.array(bg_rates)
    bg_errors = np.array(bg_errors)
    bg_fracexp = np.array(bg_fracexp)

    ##########
    if False in lc_times==bg_times:
        return 'lc_times is not equal to bg_times!'

    if False in lc_fracexp==bg_fracexp:
        return 'lc_fracexp is not aequal to bg_fracexp!'

    new_rates = lc_rates - bg_rates*bg_scale
    new_errors = np.sqrt(lc_errors**2 + (bg_scale*bg_errors)**2)

    #for i in range(100):
    #    print(i,lc_rates[i],bg_rates[i],new_rates[i],lc_errors[i],bg_errors[i],new_errors[i])

    return lc_times,new_rates,new_errors,lc_fracexp


if __name__ == "__main__":
    eventlist = glob.glob('/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/sw*pc*po_cl.evt')
    #eventlist = ['/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/sw00049834027xpcw3po_cl.evt','/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/sw00049834087xpcw3po_cl.evt']
    bary_outputfolder = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/lightcurve/'
    trunc_eventlist = time_order(eventlist,57739.00,58238.99) #originally was 58239 to 58605

    do_barycorr = False

    if do_barycorr == True:
        print('Doing barycenter corrections...')
        for i in tqdm(range(len(trunc_eventlist))):
            inputfile = trunc_eventlist[i]
            obsid = str(pathlib.Path(inputfile).name)[:13]

            outputfile = bary_outputfolder + str(pathlib.Path(inputfile).name)[:20] + '_bary_cl.evt'
            orbitfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/auxil/' + obsid + 'sao.fits'

            attfile = att_file_use(trunc_eventlist[i])
            output_attfile = bary_outputfolder + str(pathlib.Path(attfile).name)[:16] + '_bary.fits'

            hkfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/hk/' + obsid + 'xhd.hk.gz'
            output_hkfile = bary_outputfolder + obsid + '_bary_xhd.hk.gz'

            ### applying barycenter corrections to the hd hk file
            barycorr(hkfile,output_hkfile,'ICRS',orbitfile,bary_outputfolder)

            ### applying barycenter corrections to the attitude file
            barycorr(attfile,output_attfile,'ICRS',orbitfile,bary_outputfolder)

            ### applying barycenter corrections to the event file
            barycorr(inputfile,outputfile,'ICRS',orbitfile,bary_outputfolder)

    bary_eventlist = glob.glob(bary_outputfolder + 'sw*bary_cl.evt')
    #xselect_script(bary_eventlist,'/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300ulx1.reg',10,57739.00,58238.99)

    #xrtlccorr(glob.glob(bary_outputfolder+'sw*bary_cl.evt'),glob.glob(bary_outputfolder+'sw*_ulx1.lc'),glob.glob(bary_outputfolder+'sw*at_bary.fits'),glob.glob(bary_outputfolder+'sw*_bary_*hd*hk*'),57739.00,58238.99)

    #obsids = [str(i) for i in range(49834027,49834042)] + [str(i) for i in range(49834043,49834062)] + [str(i) for i in range(49834063,49834066)] + ['88810002'] + [str(i) for i in range(49834066,49834069)] + [str(i) for i in range(49834070,49834079)] + [str(i) for i in range(49834080,49834088)]
    #corr_lc_files = [bary_outputfolder + 'sw000' + obsids[i] + '_corr.lc' for i in range(len(obsids))]
    #corr_bg_files = [bary_outputfolder + 'sw000' + obsids[i] + '_bg_corr.lc' for i in range(len(obsids))]
    #bg_scale = (30/120)**2

    #lcmath(corr_lc_files,corr_bg_files,bg_scale)
    #new_times,new_rates,new_errs = get_bgsub(corr_lc_files,corr_bg_files,bg_scale)
