#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:17am 2019

Program for scp (secure copy protocol) - to copy files from ciri onto the hard drive

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs
import os
import subprocess

Lv0_dirs.global_par()

def scp(obsid):
    """
    To securely copy the files from ciri onto /Volumes/Samsung_T5/NICER-data/

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")

    obsid_dir = 'masonng@ciri:/nfs/ciri/nicer/decrypted/' + obsid
    nicer_dir = Lv0_dirs.NICER_DATADIR

    if os.path.isdir(nicer_dir+obsid): #to prevent duplicate files ; not likely to be the case, but just in case...
        subprocess.check_call(['rm','-r',nicer_dir+obsid])

    subprocess.check_call(['scp','-r',obsid_dir,nicer_dir])
    return

##### for i in range(11,25):
##    scp('10600601' + str(i)) WORKED.

#scp('1060060127')
