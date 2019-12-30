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
from tqdm import tqdm
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

    for i in range(1030180101,1030180188):
        psrpipe(str(i),['--emin','0.3','--emax','12.0','--mask','14','34','54'])
    #obsids = [str(i) for i in range(1060060101,1060060200)]# + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
    #obsids = [str(i) for i in range(1060060251,1060060300)]
    #[str(i) for i in range(1060060301,1060060313)] +
    #for i in range(len(obsids)):
    #obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113']
    #obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113'] + [str(i) for i in range(1060060101,1060060200)] + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
    #logfile = 'psrpipe_merge.log'
    #with open(logfile,'w') as logtextfile:
    #    logtextfile.write(subprocess.check_output(['psrpipe.py',Lv0_dirs.NICER_DATADIR+obsids[0],'--emin','0.3','--emax','8','--mask','14','34','54','--cormin','1.5','--merge','--crcut']))
    #    logtextfile.close()

    #for i in tqdm(range(len(obsids))):
    #    psrpipe(obsids[i],['--emin','0.3','--emax','12.0','--mask','14','34','54','--cormin','1.5'])

    #psrpipe('1013010105',['--emin','0.3','--emax','12.0'])
