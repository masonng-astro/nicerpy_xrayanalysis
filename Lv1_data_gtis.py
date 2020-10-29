#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 4:55pm 2020

Script to manipulate GTIs where necessary. In the Feb 22 2020 edition, I will
be combining GTIs, such that consecutive GTIs that are a short time apart will be combined.
This is so that I can run acceleration searches on individual GTIs!

https://stackoverflow.com/questions/3704918/python-way-to-restart-a-for-loop-similar-to-continue-for-while-loops
For GTI_bunching
"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import pathlib
import subprocess
import Lv0_dirs,Lv0_fits2dict

Lv0_dirs.global_par() #obtaining the global parameters

def GTI_bunching(eventfile,gap,gtifile):
    """
    To bunch up GTIs which are separated by less than a few seconds apart.

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    gap - maximum separation between end time of 1st GTI and start time of 2nd GTI allowed
    gtifile - name of GTI file
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    gtifile_path = parent_folder + '/' + gtifile
    gtis = list(fits.open(eventfile)[2].data)
    N = len(gtis)

    should_restart = True
    while should_restart: #to restart the for loop if I reached a pair of GTIs that can be combined
        should_restart = False
        for i in range(N-1): #variable N is there as we need to UPDATE the length of the GTIs each time we delete one!
            if (gtis[i+1][0]-gtis[i][1]) <= gap: #if (start time of next GTI - end time of previous GTI) is less than 5s, say
                new_start = gtis[i][0] #defined here, because "i" will change after deleting the 2 previous GTIs!
                new_end = gtis[i+1][1]
                del gtis[i]
                del gtis[i] #use the SAME index since the list after the first deletion will have N-1 items!
                gtis.insert(i,(new_start,new_end))
                N = len(gtis)
                should_restart = True
                break
            else:
                N = len(gtis)
                continue

    gtitxt = open(parent_folder + '/bunchedgtis.txt','w')
    for i in range(len(gtis)):
        gtitxt.write(str(gtis[i][0]) + ' ' + str(gtis[i][1]) + '\n')
    gtitxt.close()

    gti_col_file = Lv0_dirs.NICER_DATADIR + 'gti_columns.txt'
    gti_header_file = Lv0_dirs.NICER_DATADIR + 'gti_header.txt'
    ftcreate_cmd = ['ftcreate',gti_col_file,parent_folder+'/bunchedgtis.txt',gtifile_path,'headfile='+gti_header_file,'extname="GTI"','clobber=yes']
    subprocess.run(ftcreate_cmd)

if __name__ == "__main__":
    eventfile = Lv0_dirs.NICER_DATADIR + 'xtej1739/2002131540_bary.evt'
    gtifile = Lv0_dirs.NICER_DATADIR + 'xtej1739/bunched.gti'
    GTI_bunching(eventfile,5,gtifile)
