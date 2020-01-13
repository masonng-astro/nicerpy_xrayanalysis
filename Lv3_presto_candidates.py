#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 5 1:28pm 2019

Script to automate getting pulsation candidates of a certain frequency range,
and reporting other germane information?

"""
from __future__ import division, print_function
import numpy as np
from scipy import stats, signal
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import Lv0_dirs

Lv0_dirs.global_par()

def get_candidates(obsid,name_par_list,zmax,f1,f2):
    """
    Getting pulsation candidates within some frequency range. If I want the full
    frequency range, just do f1 = 0, f2 = some large number.

    obsid - Observation ID of the object of interest (10-digit str)
    name_par_list - list of parameters specifying parameters like GTI number and/or energy range
    zmax - maximum acceleration
    f1 - lower limit of frequency range
    f2 - upper limit of frequency range

    name_par_list should be [GTI_true,E_true,GTIno,segment_length,PI1,PI2]
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if type(name_par_list) != list and type(name_par_list) != np.ndarray:
        raise TypeError("name_par_list should either be a list or an array!")
    if len(name_par_list) != 6:
        raise ValueError("There seems to be fewer or more values in the list/array than there should be! You should have [GTI_true, E_true, GTIno, segment length, PI1, PI2]")

    header1 = "             Summed  Coherent  Num        Period          Frequency         FFT 'r'        Freq Deriv       FFT 'z'         Accel                           "
    header2 = "                        Power /          Raw           FFT 'r'          Pred 'r'       FFT 'z'     Pred 'z'      Phase       Centroid     Purity                        "

    obsid_file = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt'
    header_card = fits.open(obsid_file)[0].header
    date_obs = str(header_card['DATE-OBS'])
    date_end = str(header_card['DATE-END'])
    tstart = str(header_card['TSTART'])
    tstop = str(header_card['TSTOP'])

    if name_par_list[0] == True and name_par_list[1] == False: #if we're looking at just time segments!
        working_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(name_par_list[3]) + 's/'
        ACCEL_files = sorted(glob.glob(working_dir+'*_' + str(name_par_list[3]) + 's_ACCEL_' + str(zmax)))
        cands_txt = open(working_dir+'candidates_'+str(name_par_list[3])+'s_raw.txt','w')
        cands_txt.write('ObsID      Start-date/time-of-obs  End-date/time-of-obs    MET_start        MET_end     seg_no  MET_centroid  Cand_No  Sigma  Frequency  Freq_Deriv       z       Acceleration' + '\n')
        """
        JULY 8: Got to edit the below as appropriate! Mainly got to think about how to replace seg_no!
        ACTUALLY, BREAK THIS UP INTO 3 SEPARATE FUNCTIONS! Integrate into Lv3_presto_main as well...
    elif name_par_list[0] == False and name_par_list[1] == True: #if we're looking at just energy segments!
        working_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/'
        ACCEL_files = sorted(glob.glob(working_dir+'*E'+str(name_par_list[4]) + '-' + str(name_par_list[5])))
        cands_txt = open(working_dir+'candidates_raw.txt','w')
        cands_txt.write('ObsID      Start-date/time-of-obs  End-date/time-of-obs    MET_start        MET_end     seg_no  MET_centroid  Cand_No  Sigma  Frequency  Freq_Deriv       z       Acceleration' + '\n')
    else: #if we're looking at BOTH time AND energy segments!
        working_dir = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/accelsearch_' + str(name_par_list[3]) + 's/'
        ACCEL_files = sorted(glob.glob(working_dir+'*_' + str(name_par_list[3]) + 's_ACCEL_' + str(zmax)))
        cands_txt = open(working_dir+'candidates_raw.txt','w')
        cands_txt.write('ObsID      Start-date/time-of-obs  End-date/time-of-obs    MET_start        MET_end     seg_no  MET_centroid  Cand_No  Sigma  Frequency  Freq_Deriv       z       Acceleration' + '\n')
        """
    for i in range(len(ACCEL_files)):
        accel_textfile = np.array(open(ACCEL_files[i],'r').read().split('\n')) #read the data from the ACCEL_$zmax files
        index_header1 = np.where(accel_textfile==header1)[0][0] #index for "Summed, Coherent, Num, Period etc
        index_header2 = np.where(accel_textfile==header2)[0][0] #index for "Power / Raw  FFT  'r'  etc
        no_cands = index_header2 - index_header1 - 5 #to obtain number of candidates

        segment_no = '0004'
        MET_centroid = '141080121.942' #test
        candidates = np.genfromtxt(ACCEL_files[i],dtype='str',skip_header=3,usecols=(0,1,6,8,9,10),unpack=True,max_rows=no_cands)
        if len(candidates) == candidates.size: #meaning if there's ONE pulsation candidate in the *ACCEL_$zmax file
            if (np.float(candidates[2][:-3]) > f1) and (np.float(candidates[2][:-3]) < f2):
                cands_txt.write(obsid + '   ' + date_obs + '   ' + date_end + '   ' + tstart + '   ' + tstop + '   ' + segment_no + '   ' + MET_centroid + '   ' + candidates[0].zfill(4) + '   ' + candidates[1] + '   ' + candidates[2] + '   ' + candidates[3] + '   ' + candidates[4] + '   ' + candidates[5] + '\n')
        else: #if there are multiple pulsation candidates in the *ACCEL_$zmax file
            for j in range(candidates.shape[1]): #for EACH pulsation candidate
                if (np.float(candidates[2][j][:-3]) > f1) and (np.float(candidates[2][j][:-3]) < f2):
                    cands_txt.write(obsid + '   ' + date_obs + '   ' + date_end + '   ' + tstart + '   ' + tstop + '   ' + segment_no + '   ' + MET_centroid + '   ' + candidates[0][j].zfill(4) + '   ' + candidates[1][j] + '   ' + candidates[2][j] + '   ' + candidates[3][j] + '   ' + candidates[4][j] + '   ' + candidates[5][j] + '\n')

if __name__ == "__main__":
    get_candidates('1200250101',[True,False,0,64,0,0],100,0,100)
