#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:27am 2019

Program for psrpipe

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs
import os
import subprocess

Lv0_dirs.global_par()

def psrpipe(obsid,flags):
    """
    Running psrpipe on the observation, to make more cuts! I decided not to
    put in pre-determined options for fully flexibility. Though standard flags
    would be ['--emin','0.3','--emax','12.0','--shrinkelvcut'], though there are
    others. Check out "psrpipe.py -h"! Also made sure that I moved $OBSID_pipe from
    the working directory to where NICERSOFT_DATADIR is, though I need to temporarily
    store the output folder in the working directory.

    obsid - Observation ID of the object of interest (10-digit str)
    flags - a LIST of input flags for psrpipe
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    ## Convert the .mkf file from the niprefilter to niprefilter2 version
    indir_nip2 = 'indir='+Lv0_dirs.NICER_DATADIR+obsid
    infile_nip2 = 'infile='+Lv0_dirs.NICER_DATADIR+obsid+'/auxil/ni'+obsid+'.mkf'
    outfile_nip2 = 'outfile='+Lv0_dirs.NICER_DATADIR+obsid+'/auxil/ni'+obsid+'.mkf'
    logfile = obsid + '_psrpipe.log'

    subprocess.check_call(['cp',Lv0_dirs.NICER_DATADIR+obsid+'/auxil/ni'+obsid+'.mkf',Lv0_dirs.NICER_DATADIR+obsid+'/auxil/oldmkf_'+obsid])
    subprocess.check_call(['niprefilter2',indir_nip2,infile_nip2,outfile_nip2,'clobber=YES'])

    if os.path.isdir(Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe'): #to prevent duplicate files ; not likely to be the case, but just in case...
        subprocess.check_call(['rm','-r',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe'])

    with open(logfile,'w') as logtextfile:
        command = ['psrpipe.py',Lv0_dirs.NICER_DATADIR+obsid] + flags
        logtextfile.write(subprocess.check_output(command))
        logtextfile.close()

    if os.path.isdir(Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe'): #to prevent duplicate files ; not likely to be the case, but just in case...
        subprocess.check_call(['rm','-r',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe'])

    tools_dir = '/Users/masonng/Documents/MIT/Research/nicerpy_xrayanalysis/'
    subprocess.check_call(['mv',tools_dir+obsid+'_pipe',Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe'])

    return

if __name__ == "__main__":
    #psrpipe('1060060127',['--emin','0.3','--emax','12.0','--shrinkelvcut'])
    for i in range(1030180133,1030180137):
        psrpipe(str(i),['--emin','0.3','--emax','12.0','--mask','14','34','54'])
    #psrpipe('1013010105',['--emin','0.3','--emax','12.0'])
