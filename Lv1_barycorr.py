#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 1:37pm 2019

Program for barycorr - doing barycenter corrections to the data

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv0_call_eventcl
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

def

if __name__ == "__main__":
    barycorr('1060060127','ICRS')
