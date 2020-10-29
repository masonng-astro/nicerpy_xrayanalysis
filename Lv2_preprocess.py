#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 2:22pm 2019

Program for doing Lv0_scp, Lv0_gunzip, Lv0_psrpipe, Lv1_barycorr, that is,
this is pre-processing data before running them through to PRESTO!

(update documentation later - meant for archival data, really! Use Lv3_quicklook.py or equivalent for
the more urgent analyses)

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv0_gunzip,Lv0_nicerl2,Lv0_psrpipe,Lv1_barycorr
import os
from tqdm import tqdm
import subprocess
import glob

Lv0_dirs.global_par()

def preprocess(obsdir,nicerl2_flags,psrpipe_flags,refframe,orbitfile,parfile,nicer_datafile,nicer_output,nicersoft_datafile,nicersoft_output,nicersoft_folder,custom_coords):
    """
    Preprocessing the NICER data for use in PRESTO, so running gunzip, psrpipe, and barycorr.

    obsdir - NICER data directory containing all the data files (e.g., path_to_NICER_dir/1034070101)
    nicerl2_flags - a LIST of input flags for nicerl2
    psrpipe_flags - a LIST of input flags for psrpipe
    refframe - reference frame for barycenter corrections (usually ICRS)
    orbitfile - orbit file (from the NICER data directory) for the barycenter corrections
    parfile - .par file if any
    nicer_datafile - event file from the original NICER data directory
    nicer_output - output/barycenter-corrected event file from the original NICER data directory
    nicersoft_datafile - event file (usually cleanfilt.evt) from the NICERsoft data directory
    nicersoft_output - output/barycenter-corrected event file from the NICERsoft data directory
    nicersoft_folder - output folder for the NICERsoft data
    custom_coords - either an empty list/array or a list/array with two elements (RA/DEC)
    """
    if type(psrpipe_flags) != list:
        raise TypeError("flags should be a list! Not even an array.")
    if type(nicerl2_flags) != list:
        raise TypeError("flags should be a list! Not even an array.")
    if refframe != 'ICRS' and refframe != 'FK5':
        raise ValueError("refframe should either be ICRS or FK5! Otherwise, update Lv1_barycorr.py if there are options I was unaware of.")

    print('Now unzipping all the files!')
    Lv0_gunzip.unzip_all(obsdir) #unzipping the contents within the observation
    print('Now running nicerl2!')
    Lv0_nicerl2.nicerl2(obsdir,nicerl2_flags)
    print('Now running psrpipe.py from NICERsoft!')
    Lv0_psrpipe.psrpipe(nicer_datafile,psrpipe_flags) #applying custom cuts (though no need --shrinkelv after HEASOFT 6.26)

    print('Now running barycorr from HEASOFT!')
    ##### For the NICER data
    Lv1_barycorr.barycorr(nicer_datafile,nicer_output,refframe,orbitfile,parfile,obsdir,custom_coords)
    ##### For the NICERsoft file (cleanfilt.evt, usually)
    Lv1_barycorr.barycorr(nicersoft_datafile,nicersoft_output,refframe,orbitfile,parfile,nicersoft_folder,custom_coords)

    return

if __name__ == "__main__":
    #obsdirs = [Lv0_dirs.NICER_DATADIR + str(i) + '/' for i in range(1034200124,1034200200)] + [Lv0_dirs.NICER_DATADIR + str(i) + '/' for i in range(1034200201,1034200241)] + [Lv0_dirs.NICER_DATADIR + str(i) + '/' for i in range(2034200201,2034200206)]
    obsdirs = [Lv0_dirs.NICER_DATADIR + '/2584010501/']
    nicerl2_flags = ['clobber=YES']
    psrpipe_flags = ['--emin','0.3','--emax','12.0'] #for psrpipe in Lv0_psrpipe
    refframe = 'ICRS' #for barycorr in Lv1_barycorr
    orbitfile = [obsdirs[i] + 'auxil/ni' + obsdirs[i][-11:-1] + '.orb' for i in range(len(obsdirs))]
    parfile = ''
    custom_coords = np.array([])

    nicer_datafile = [obsdirs[i] + 'xti/event_cl/ni' + obsdirs[i][-11:-1] + '_0mpu7_cl.evt' for i in range(len(obsdirs))]
    nicer_output = [obsdirs[i] + 'xti/event_cl/ni' + obsdirs[i][-11:-1] + '_0mpu7_cl_bary.evt' for i in range(len(obsdirs))]
    nicersoft_datafile = [Lv0_dirs.NICERSOFT_DATADIR + obsdirs[i][-11:-1] + '_pipe/cleanfilt.evt' for i in range(len(obsdirs))]
    nicersoft_output = [Lv0_dirs.NICERSOFT_DATADIR + obsdirs[i][-11:-1] + '_pipe/ni' + obsdirs[i][-11:-1] + '_nicersoft_bary.evt' for i in range(len(obsdirs))]
    nicersoft_folder = [Lv0_dirs.NICERSOFT_DATADIR + obsdirs[i][-11:-1] + '_pipe/' for i in range(len(obsdirs))]

    for i in tqdm(range(len(obsdirs))):
        preprocess(obsdirs[i],nicerl2_flags,psrpipe_flags,refframe,orbitfile[i],parfile,nicer_datafile[i],nicer_output[i],nicersoft_datafile[i],nicersoft_output[i],nicersoft_folder[i],custom_coords)
