#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 1:37pm 2019

Program for barycorr - doing barycenter corrections to the data

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs
import subprocess

Lv0_dirs.global_par()

def get_ra_dec(eventfile):
    """
    Obtain the RA_OBJ and DEC_OBJ corresponding to the observation!

    obsid - Observation ID of the object of interest (10-digit str)
    """
    event = fits.open(eventfile)
    event_header = event[1].header

    return event_header['RA_OBJ'], event_header['DEC_OBJ']

def read_par(parfile):
    """
    Function that reads a par file. In particular, for the purposes of barycorr,
    it will return POSEPOCH, RAJ, DECJ, PMRA, and PMDEC.

    Step 1: Read par file line by line, where each line is stored as a string in the 'contents' array
    Step 2a: For PSRJ, RAJ, DECJ, PMRA, and PMDEC, those lines are teased out
    Step 2b: The corresponding strings are split up without whitespace
    Step 3: Extract the values accordingly
    """
    if parfile[-4:] != '.par':
        raise ValueError("parfile is neither an empty string nor a .par file. Is this right?")

    contents = open(parfile,'r').read().split('\n')
    posepoch = [contents[i] for i in range(len(contents)) if 'POSEPOCH' in contents[i]][0].split()
    raj = [contents[i] for i in range(len(contents)) if 'RAJ' in contents[i]][0].split()[1].split(':')
    decj = [contents[i] for i in range(len(contents)) if 'DECJ' in contents[i]][0].split()[1].split(':')
    pmra = [contents[i] for i in range(len(contents)) if 'PMRA' in contents[i]][0].split()[1] #milliarcsec per year
    pmdec = [contents[i] for i in range(len(contents)) if 'PMDEC' in contents[i]][0].split()[1] #milliarcsec per year

    ra_deg = np.float(raj[0])*15 + np.float(raj[1])/60 * 15 + np.float(raj[2])/3600 * 15 #RA in HH:MM:SS to deg
    dec_deg = np.float(decj[0]) + np.float(decj[1])/60 + np.float(decj[2])/3600 #DEC in HH:MM:SS to deg

    return posepoch[1], ra_deg, dec_deg, pmra, pmdec #returns object name, RA, DEC, PMRA, PMDEC

def barycorr(eventfile,outfile,refframe,orbit_file,parfile,output_folder):
    """
    General function to perform the barycenter corrections for an event file

    obsid - Observation ID of the object of interest (10-digit str)
    refframe - reference frame for barycenter corrections (usually ICRS)
    parfile - name of the .par file
    """
    if refframe != 'ICRS' and refframe != 'FK5':
        raise ValueError("refframe should either be ICRS or FK5! Otherwise, update Lv1_barycorr.py if there are options I was unaware of.")

    TIMEZERO = fits.open(eventfile)[1].header['TIMEZERO']
    logfile = output_folder + 'barycorr_notes.txt'

    if parfile[-4:] == '.par': #i.e., if we have a par file:
        posepoch,old_ra,old_dec,PMRA,PMDEC = read_par(parfile)
        #get epoch of position, old RA, old DEC, and PMRA/PMDEC proper motion corrections

        gti_card = fits.open(eventfile)[2]
        gtis = gti_card.data #aim is to get start/stop times of observation to get time of middle of observation
        MJDREFI = gti_card.header['MJDREFI']
        MJDREFF = gti_card.header['MJDREFF']
        centroid_time = (gtis[-1][1]-gtis[0][0])/2 #not going to take TIMEZERO into account because difference is small? Plus it's a -1 - (-1) thing

        centroid_MJD = (MJDREFI+MJDREF) + (TIMEZERO+centroid_time)/86400 #convert centroid MET to MJD
        time_elapsed = centroid_MJD-POSEPOCH #centroid_MJD is later than POSEPOCH!

        new_ra = old_ra + (PMRA/np.cos(old_dec*np.pi/180) * 1E-3/3600) * time_elapsed/365.2425 #365.2425 days in a year
        new_dec = old_dec + (PMDEC*1E-3/3600) * time_elapsed/365.2425

        with open(logfile,'w') as logtextfile:
            output = subprocess.run(['barycorr',eventfile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(new_ra),'dec='+str(new_dec),'refframe='+str(refframe),'clobber=YES'],capture_output=True,text=True)
            logtextfile.write(output.stdout)
            logtextfile.write('*------------------------------* \n')
            logtextfile.write(output.stderr)
            logtextfile.close()

    elif parfile == '':
        ra,dec = get_ra_dec(eventfile)
        with open(logfile,'w') as logtextfile:
            output = subprocess.run(['barycorr',eventfile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(ra),'dec='+str(dec),'refframe='+str(refframe),'clobber=YES'],capture_output=True,text=True)
            logtextfile.write(output.stdout)
            logtextfile.write('*------------------------------* \n')
            logtextfile.write(output.stderr)
            logtextfile.close()

    else:
        raise ValueError("parfile is neither an empty string nor a .par file. Is this right?")

if __name__ == "__main__":
    #print(read_par('/Users/masonng/Downloads/test.par'))
    obsid = '1911212213'
    eventfile = Lv0_dirs.NICER_DATADIR + 'rxj0209/rxj0209kgfilt.evt'
    outfile = Lv0_dirs.NICER_DATADIR + 'rxj0209/rxj0209kgfilt_bary.evt'
    orbitfile = Lv0_dirs.NICER_DATADIR + 'rxj0209/rxj0209.orb'
    parfile = ''
    output_folder = Lv0_dirs.NICER_DATADIR + 'rxj0209/'
    refframe = 'ICRS'

    barycorr(eventfile,outfile,refframe,orbitfile,parfile,output_folder)
