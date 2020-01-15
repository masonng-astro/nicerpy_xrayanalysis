#!/usr/bin/env python
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

    obsid_dir = 'masonng@ciri.mit.edu:/nfs/ciri/nicer/decrypted/' + obsid
    nicer_dir = Lv0_dirs.NICER_DATADIR

    if os.path.isdir(nicer_dir+obsid): #to prevent duplicate files ; not likely to be the case, but just in case...
        subprocess.check_call(['rm','-r',nicer_dir+obsid])

    subprocess.check_call(['scp','-r',obsid_dir,nicer_dir])
    return

if __name__ == "__main__":
    obsids = ['203420020' + str(i) for i in range(1,6)] + ['103419010' + str(i) for i in range(1,5)] + ['1034200' + str(i) for i in range(201,241)]
    for i in range(len(obsids)):
        scp(obsids[i])

##### for i in range(11,25):
##    scp('10600601' + str(i)) WORKED.
