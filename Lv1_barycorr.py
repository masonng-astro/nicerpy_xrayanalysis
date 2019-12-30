#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 1:37pm 2019

Program for barycorr - doing barycenter corrections to the data

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
from tqdm import tqdm
import Lv0_dirs,Lv0_call_eventcl,Lv0_call_nicersoft_eventcl
import subprocess

Lv0_dirs.global_par()

def get_ra_dec(obsid):
    """
    Obtain the RA_OBJ and DEC_OBJ corresponding to the observation!

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    event = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/cleanfilt.evt'
    event = fits.open(event)
    event_header = event[1].header

    return event_header['RA_OBJ'], event_header['DEC_OBJ']

def nicerdata_barycorr(obsid,refframe):
    """
    Applying the barycenter corrections to the X-ray timing data (for NICER in this case)

    obsid - Observation ID of the object of interest (10-digit str)
    refframe - reference frame for barycenter corrections (usually ICRS)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if refframe != 'ICRS' and refframe != 'FK5':
        raise ValueError("refframe should either be ICRS or FK5! Otherwise, update Lv1_barycorr.py if there are options I was unaware of.")

    event = Lv0_call_eventcl.open_fits(obsid,False)
    event_header = event[1].header
    ra,dec = event_header['RA_OBJ'], event_header['DEC_OBJ']

    output_folder = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_cl/'
    infile = output_folder + 'ni' + obsid + '_0mpu7_cl.evt'
    outfile = output_folder + 'ni' + obsid + '_0mpu7_cl_bary.evt'
    orbit_file = Lv0_dirs.NICER_DATADIR + obsid + '/auxil/' + 'ni' + obsid + '.orb'
    logfile = output_folder + 'barycorr_notes.txt'

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['barycorr',infile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(ra),'dec='+str(dec),'refframe='+str(refframe),'clobber=YES']))
        logtextfile.close()

    return

def nicersoft_barycorr(obsid,refframe):
    """
    Applying the barycenter corrections to the X-ray timing data (for NICER in this case)

    obsid - Observation ID of the object of interest (10-digit str)
    refframe - reference frame for barycenter corrections (usually ICRS)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if refframe != 'ICRS' and refframe != 'FK5':
        raise ValueError("refframe should either be ICRS or FK5! Otherwise, update Lv1_barycorr.py if there are options I was unaware of.")

    ra,dec = get_ra_dec(obsid)
    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    infile = nicersoft_output_folder + 'cleanfilt.evt'
    outfile = nicersoft_output_folder + 'ni' + obsid + '_nicersoft_bary.evt'
    orbit_file = nicersoft_output_folder + 'ni' + obsid + '.orb'
    logfile = nicersoft_output_folder + 'barycorr_notes.txt'

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['barycorr',infile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(ra),'dec='+str(dec),'refframe='+str(refframe),'clobber=YES']))
        logtextfile.close()

    return

def B1957_20_barycorr(obsid):
    """
    Applying barycenter corrections to the B1957+20 event files! Want to do barycenter
    corrections ObsID by ObsID, i.e., take into account PMRA, PMDEC when using the RA
    and DEC for barycenter corrections.

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    POSEPOCH = 58000 #reference epoch for position
    old_ra = 19*15 + 59/60 * 15 + 36.7404469/3600 * 15 # 299.90317048875
    old_dec = 20 + 48/60 + 14.42810/3600 # 20.804141347222224
    PMRA = -15.428307750281743369
    PMDEC = -24.933123414972509821

    gti_card = Lv0_call_eventcl.open_fits(obsid,False)

    MJDREFI = gti_card[2].header['MJDREFI']
    MJDREFF = gti_card[2].header['MJDREFF']

    gtis = gti_card[2].data #for GTIs
    centroid_time = (gtis[-1][1]-gtis[0][0])/2

    time_elapsed = centroid_time-POSEPOCH

    new_ra = old_ra + (PMRA/np.cos(old_dec*np.pi/180) * 1E-3/3600) * time_elapsed/365.2425 #365.2425 days in a year
    new_dec = old_dec + (PMDEC*1E-3/3600) * time_elapsed/365.2425

    refframe = 'ICRS'
    output_folder = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_cl/'
    infile = output_folder + 'ni' + obsid + '_0mpu7_cl.evt'
    outfile = output_folder + 'ni' + obsid + '_0mpu7_cl_bary.evt'
    orbit_file = Lv0_dirs.NICER_DATADIR + obsid + '/auxil/' + 'ni' + obsid + '.orb'
    logfile = output_folder + 'barycorr_notes.txt'

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['barycorr',infile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(new_ra),'dec='+str(new_dec),'refframe='+str(refframe),'clobber=YES']))
        logtextfile.close()

    return

def B1957_20_nicersoft_barycorr(obsid):
    """
    Applying barycenter corrections to the B1957+20 event files! Want to do barycenter
    corrections ObsID by ObsID, i.e., take into account PMRA, PMDEC when using the RA
    and DEC for barycenter corrections.

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    POSEPOCH = 58000 #reference epoch for position
    old_ra = 19*15 + 59/60 * 15 + 36.7404469/3600 * 15 # 299.90317048875
    old_dec = 20 + 48/60 + 14.42810/3600 # 20.804141347222224
    PMRA = -15.428307750281743369
    PMDEC = -24.933123414972509821

    gti_card = Lv0_call_eventcl.open_fits(obsid,False)

    MJDREFI = gti_card[2].header['MJDREFI']
    MJDREFF = gti_card[2].header['MJDREFF']

    gtis = gti_card[2].data #for GTIs
    centroid_time = (gtis[-1][1]-gtis[0][0])/2

    time_elapsed = centroid_time-POSEPOCH

    new_ra = old_ra + (PMRA/np.cos(old_dec*np.pi/180) * 1E-3/3600) * time_elapsed/365.2425 #365.2425 days in a year
    new_dec = old_dec + (PMDEC*1E-3/3600) * time_elapsed/365.2425

    refframe = 'ICRS'
    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    infile = nicersoft_output_folder + 'cleanfilt.evt'
    outfile = nicersoft_output_folder + 'ni' + obsid + '_nicersoft_bary.evt'
    orbit_file = nicersoft_output_folder + 'ni' + obsid + '.orb'
    logfile = nicersoft_output_folder + 'barycorr_notes.txt'

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['barycorr',infile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(new_ra),'dec='+str(new_dec),'refframe='+str(refframe),'clobber=YES']))
        logtextfile.close()

    return

def J1231_nicersoft_barycorr(obsid):
    """
    Applying barycenter corrections to the J1231-1411 event files! Want to do barycenter
    corrections ObsID by ObsID, i.e., take into account PMRA, PMDEC when using the RA
    and DEC for barycenter corrections.

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    POSEPOCH = 55000 #reference epoch for position
    old_ra = 12*15 + 31/60 * 15 + 11.31331900/3600 * 15
    old_dec = -14 + 11/60 + 43.63908325/3600
    PMRA = -62.6015823614
    PMDEC = 6.75015701922

    gti_card = Lv0_call_eventcl.open_fits(obsid,False)

    MJDREFI = gti_card[2].header['MJDREFI']
    MJDREFF = gti_card[2].header['MJDREFF']

    gtis = gti_card[2].data #for GTIs
    centroid_time = (gtis[-1][1]-gtis[0][0])/2

    time_elapsed = centroid_time-POSEPOCH

    new_ra = old_ra + (PMRA/np.cos(old_dec*np.pi/180) * 1E-3/3600) * time_elapsed/365.2425 #365.2425 days in a year
    new_dec = old_dec + (PMDEC*1E-3/3600) * time_elapsed/365.2425

    refframe = 'ICRS'
    nicersoft_output_folder = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
    infile = nicersoft_output_folder + 'cleanfilt.evt'
    outfile = nicersoft_output_folder + 'ni' + obsid + '_nicersoft_bary.evt'
    orbit_file = nicersoft_output_folder + 'ni' + obsid + '.orb'
    logfile = nicersoft_output_folder + 'barycorr_notes.txt'

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['barycorr',infile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(new_ra),'dec='+str(new_dec),'refframe='+str(refframe),'clobber=YES']))
        logtextfile.close()

    return

if __name__ == "__main__":
    #obsids = ['10301801'+str(i).zfill(2) for i in range(1,32)] + ['10301801'+str(i) for i in range(49,88)]
    #for i in range(len(obsids)):
    #    B1957_20_barycorr(obsids[i])
    #    B1957_20_nicersoft_barycorr(obsids[i])
    #barycorr('1060060127','ICRS')
    obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113'] + [str(i) for i in range(1060060101,1060060200)] + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
    for i in tqdm(range(len(obsids))):
        #nicerdata_barycorr(obsids[i],'ICRS')
        #nicersoft_barycorr(obsids[i],'ICRS')
        J1231_nicersoft_barycorr(obsids[i])
