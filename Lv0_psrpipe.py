#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:27am 2019

Program for psrpipe

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs
from astropy.io import fits
from tqdm import tqdm
import os
import subprocess

Lv0_dirs.global_par()

def psrpipe(eventfile,flags):
    """
    Running psrpipe on the observation, to make more cuts! I decided not to
    put in pre-determined options for fully flexibility. Though standard flags
    would be ['--emin','0.3','--emax','12.0','--shrinkelvcut'], though there are
    others. Check out "psrpipe.py -h"! Also made sure that I moved $OBSID_pipe from
    the working directory to where NICERSOFT_DATADIR is, though I need to temporarily
    store the output folder in the working directory.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    flags - a LIST of input flags for psrpipe
    """
    if type(eventfile) != str:
        raise TypeError("eventfile should be a string!")
    if type(flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    event = fits.open(eventfile)
    obsid = event[0].header['OBS_ID']
    logfile = obsid + '_psrpipe.log'

    if os.path.isdir(Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe'): #to prevent duplicate files ; not likely to be the case, but just in case...
        subprocess.run(['rm','-r',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe'])

    with open(logfile,'w') as logtextfile:
        command = ['psrpipe.py',Lv0_dirs.NICER_DATADIR+obsid] + flags
        output = subprocess.run(command,capture_output=True,text=True)
        logtextfile.write(output.stdout)
        logtextfile.write('*------------------------------* \n')
        logtextfile.write(output.stderr)
        logtextfile.close()

    subprocess.run(['mv',obsid+'_pipe/',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/'])
    subprocess.run(['mv',obsid+'_psrpipe.log',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/'])
    #done later because the log file was created BEFORE $OBSID_pipe is created

    return

if __name__ == "__main__":
    #obsid = '1034070101'
    #eventfile = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_cl/ni' + obsid+'_0mpu7_cl.evt'
    #psrpipe(eventfile,['--emin','0.3','--emax','12.0'])

    eventfiles = [Lv0_dirs.NICER_DATADIR + str(i) + '/xti/event_cl/ni' + str(i) + '_0mpu7_cl.evt' for i in range(1030180101,1030180188)]
    for i in tqdm(range(len(eventfiles))):
        psrpipe(eventfiles[i],['--emin','0.3','--emax','12.0','--nounderfilt'])
