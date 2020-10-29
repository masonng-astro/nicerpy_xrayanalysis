#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 3:35pm 2019

Getting averaged power spectra from M segments to MERGED OBSERVATIONS, where the data
was pre-processed using NICERsoft!

By MERGED OBSERVATIONS, we mean event files that have been created from merging
event times from multiple ObsIDs!

"""
from __future__ import division, print_function
import numpy as np
from scipy import stats, signal
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
#from presto import binary_psr
import subprocess
import glob
import os
import Lv0_dirs,Lv2_mkdir

Lv0_dirs.global_par()

def merging(obsids,merged_id):
    """
    Given a list of ObsIDs, create a file which merges all the rows in the EVENTS
    extension of the FITS files!
    Updated 8/14/2020: Use a list of event files given in a text file!

    obsids - list (or array) of ObsIDs
    merged_id - 6-digit ID for the merged event file
    """
    if type(obsids) != list and type(obsids) != np.ndarray:
        raise TypeError("obsids should either be a list or an array!")
    if type(merged_id) != str or len(merged_id) != 6:
        raise TypeError("merged_id should be a string and have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/'
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    Lv2_mkdir.makedir(merged_dir)

    filelist = np.array([Lv0_dirs.NICERSOFT_DATADIR + str(obsids[i]) + '_pipe/ni' + str(obsids[i]) + '_nicersoft_bary.evt' for i in range(len(obsids))])
    int_text = open(merged_dir+'filelist.txt','w')
    for i in range(len(filelist)): #writing the paths to each file in filelist
        int_text.write(filelist[i] + '\n')
    int_text.close()

    subprocess.run(['ftmerge','@'+merged_dir+'filelist.txt',merged_dir+'intermediate.evt'])

    gti_filelist = np.array([Lv0_dirs.NICERSOFT_DATADIR + str(obsids[i]) + '_pipe/ni' + str(obsids[i]) + '_nicersoft_bary.evt[GTI]' for i in range(1,len(obsids))])
    final_text = open(merged_dir+'filelist_gti.txt','w') #list of files with the GTI extension referenced
    final_text.write(merged_dir+'intermediate.evt[GTI]'+'\n')
    for i in range(len(gti_filelist)):
        final_text.write(gti_filelist[i] + '\n')
    final_text.close()

    merged_file = merged_dir+'merged'+merged_id+'_nicersoft_bary.evt'
    subprocess.run(['ftmerge','@'+merged_dir+'filelist_gti.txt',merged_file,'clobber=YES'])
    subprocess.run(['rm','-r',merged_dir+'intermediate.evt'])

    merged_gtis = fits.open(merged_file)[2].data
    gti_sum = sum(np.array([merged_gtis[i][1]-merged_gtis[i][0] for i in range(len(merged_gtis))]))
    subprocess.run(['fparkey',str(gti_sum),merged_file,'EXPOSURE'])

    return

################################## DEPRECATED ##################################

"""
def merging(obsids):

    Given a list of ObsIDs, create a file which merges all the rows in the EVENTS
    extension of the FITS files!

    obsids - list (or array) of ObsIDs

    if type(obsids) != list and type(obsids) != np.ndarray:
        raise TypeError("obsids should either be a list or an array!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/'
    no_existing_merged = len(glob.glob(all_merged_dir+'merged*')) #number of existing merged event files/directories
    merged_id = str(no_existing_merged+1).zfill(6) #e.g., if there were 8 merged files, the ID for the next file is 00009
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'
    Lv2_mkdir.makedir(merged_dir)

    inputfile_string = ''
    for i in range(len(obsids)):
        abs_filename = Lv0_dirs.NICERSOFT_DATADIR + str(obsids[i]) + '_pipe/ni' + str(obsids[i]) + '_nicersoft_bary.evt'
        if i != len(obsids)-1:
            inputfile_string += abs_filename +','
        else:
            inputfile_string += abs_filename

    ##### Writing the ObsIDs into the text file
    merged_text_filename = merged_dir + 'merged' + merged_id + '.txt'
    merged_text = open(merged_text_filename,'w')
    merged_text.write('merged'+merged_id+': '+ ' '.join(obsids))
    merged_text.close()

    subprocess.run(['ftmerge',inputfile_string,merged_dir+'merged'+merged_id + '_nicersoft_bary.evt','clobber=YES'])
    #print(inputfile_string)

    return

def merging_GTIs(obsids,merged_id):

    Given a list of ObsIDs and the merged_id, create the final event file, which
    already has all the rows in the EVENTS extension of the FITS files, BUT also
    including the GTI extension of the FITS files!

    obsids - list (or array) of ObsIDs
    merged_id - 6-digit ID for the merged event file

    if type(obsids) != list and type(obsids) != np.ndarray:
        raise TypeError("obsids should either be a list or an array!")
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/'
    merged_dir = all_merged_dir + 'merged' + merged_id + '/'

    inputfile_string = merged_dir+'merged' + merged_id + '_nicersoft_bary.evt[gti],'
    for i in range(1,len(obsids)):
        abs_filename = Lv0_dirs.NICERSOFT_DATADIR + str(obsids[i]) + '_pipe/ni' + str(obsids[i]) + '_nicersoft_bary.evt[gti]'
        if i != len(obsids)-1:
            inputfile_string += abs_filename +','
        else:
            inputfile_string += abs_filename

    subprocess.run(['ftmerge',inputfile_string,merged_dir+'merged'+merged_id + '_nicersoft_bary.evt','clobber=YES'])
    #print(inputfile_string)

    return
"""

if __name__ == "__main__":
    #obsids = ['2060060363','2060060364','2060060365']
    #merged_id = '000013'

    obsids = ['2050280' +str(i) for i in range(119,125)]
    merged_id = '000020'

    merging(obsids)
    merging_GTIs(obsids,merged_id)
