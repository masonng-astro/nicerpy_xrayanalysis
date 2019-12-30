#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:30am 2019

Program for gunzip - to unzip the files that I got from ciri!

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs
import subprocess
import glob
from tqdm import tqdm
import os
import sh

Lv0_dirs.global_par()

def auxil(obsid):
    """
    Unzipping the auxil files

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    auxil_dir = Lv0_dirs.NICER_DATADIR + obsid + '/auxil/'
    if os.path.exists(auxil_dir+'ni'+obsid+'.cat'):
        subprocess.check_call(['rm','-r',auxil_dir+'ni'+obsid+'.cat']) #Jun 4 - seems like .cat.gz is no longer there?
    auxil_files = glob.glob(auxil_dir + '*.gz')
    for i in range(len(auxil_files)):
        sh.gunzip(auxil_files[i])

    return

def cl(obsid):
    """
    Unzipping the cl and ufa files

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    cl_dir = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_cl/'
    cl_files = glob.glob(cl_dir + '*.gz')
    for i in range(len(cl_files)):
        sh.gunzip(cl_files[i])

    return

def uf(obsid):
    """
    Unzipping the uf files

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    uf_dir = Lv0_dirs.NICER_DATADIR + obsid + '/xti/event_uf/'
    uf_files = glob.glob(uf_dir + '*.gz')
    for i in range(len(uf_files)):
        sh.gunzip(uf_files[i])

    return

def hk(obsid):
    """
    Unzipping the hk files

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    hk_dir = Lv0_dirs.NICER_DATADIR + obsid + '/xti/hk/'
    hk_files = glob.glob(hk_dir + '*.gz')
    for i in range(len(hk_files)):
        sh.gunzip(hk_files[i])

    return

def unzip_all(obsid):
    """
    Unzipping all the auxil/cl/uf/hk files

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if not os.path.isdir(Lv0_dirs.NICER_DATADIR + obsid):
        raise ValueError("This folder with ObsID " + obsid + " does not exist! Have you imported it with Lv0_scp.py?")

    auxil(obsid)
    cl(obsid)
    uf(obsid)
    hk(obsid)

    return

if __name__ == "__main__":
    #for i in range(1030180101,1030180188):
    #    unzip_all(str(i))
    #obsids = [str(i) for i in range(1060060101,1060060200)] + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
    #obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113']
    obsids = [str(i) for i in range(1030180101,1030180188)]
    for i in tqdm(range(len(obsids))):
        unzip_all(obsids[i])
