#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs May 30 9:31am 2019

Program for doing nicerl2 (mix of nicercal, niprefilter2, nimaketime, nicermergeclean)

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs
from tqdm import tqdm
import os
import subprocess
import glob

Lv0_dirs.global_par()

def nicerl2(obsid,nicerl2_flags):
    """
    Running nicerl2 to do initial filtering of the ufa file!

    obsid - Observation ID of the object of interest (10-digit str)
    flags - a LIST of input flags for nicerl2
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(nicerl2_flags) != list:
        raise TypeError("flags should be a list! Not even an array.")

    logfile = obsid + '_nicerl2.log'
    nicerdata_dir = Lv0_dirs.NICER_DATADIR + obsid

    with open(logfile,'w') as logtextfile:
        logtextfile.write(subprocess.check_output(['nicerl2',nicerdata_dir]+nicerl2_flags))
        logtextfile.close()

    subprocess.check_call(['mv',obsid+'_nicerl2.log',nicerdata_dir])

    return

if __name__ == "__main__":
    #for i in tqdm(range(1030180101,1030180188)):
    #    print('Doing ObsID ' + str(i))
    #    nicerl2(str(i),['clobber=YES'])
    #nicerl2('1013010105',['ang_dist=0.035'])
    #obsids = [str(i) for i in range(1060060101,1060060200)] + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
    #obsids = [str(i) for i in range(1060060224,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
    #obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113']
    #obsids = ['0060060101','0060060102','0060060103','0060060104','0060060105','0060060106','0060060107','0060060108','0060060109','0060060110','0060060111','0060060112','0060060113'] + [str(i) for i in range(1060060101,1060060200)] + [str(i) for i in range(1060060201,1060060300)] + [str(i) for i in range(1060060301,1060060313)]
    obsids = [str(i) for i in range(1030180101,1030180188)]
    bad_ids = [str(i) for i in range(1030180131,1030180148)]
    for i in tqdm(range(len(obsids))):
        if obsids[i] not in bad_ids:
            nicerl2(obsids[i],['clobber=YES','underonly_range=0-1000'])
        if obsids[i] in bad_ids:
            nicerl2(obsids[i],['clobber=YES','underonly_range=0-1000','mpulist=0,1,2,4,5,6'])
