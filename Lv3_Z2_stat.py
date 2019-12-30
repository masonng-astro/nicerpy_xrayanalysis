#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs Oct 24th 9:54am 2019

Calculates the Z^2 statistic given a list of ObsIDs, where the event files
must already have 'PULSE_PHASE' in them! The input will be a par file. Will also
have a way to iterate through a grid of orbital parameters.

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from pint.eventstats import sf_z2m,z2m,sig2sigma
import subprocess
from tqdm import tqdm
import glob
import os,shutil
import Lv0_dirs,Lv2_mkdir

Lv0_dirs.global_par()

def get_Z2(obsid,model,par_file,m):
    """
    Calculates the Z^2 statistic given the ObsID, binary model name (if there's one),
    and the parameter file.

    obsid - Observation ID of the object of interest (10-digit str)
    model - binary model being used (ELL1, BT, DD, or DDK for now)
    par_file - orbital parameter file for input into PINT's photonphase
    m - number of harmonics
    """
    filename = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/cleanfilt_phase'
    if model == '':
        filename = filename + '.evt'
    if model != 'ELL1' and model != 'BT' and model != 'DD' and model != 'DDK':
        raise ValueError("The model should be one of ELL1, BT, DD, or DDK! Update if a new model has been added to PINT.")
    else:
        filename = filename + '_' + model + '.evt'

    event = fits.open(filename)
    if 'PULSE_PHASE' not in event[1].columns.names:
        raise ValueError('The PULSE_PHASE column is not in the PINT-processed event file! Check your steps again.')

    phases = event[1].data['PULSE_PHASE']
    z_vals = z2m(phases,m=m)
    probs = sf_z2m(z_vals)
    significances = sig2sigma(probs)

    return significances

def get_merged(merged_id,obsids,model,par_file,m,ratecut):
    """
    Given the list of ObsIDs, make the merged file!
    Will have all the outputs of merge.py in the merged_id folder.

    merged_id - 6-digit ID for the merged event file
    obsids - list (or array) of ObsIDs
    model - binary model being used (ELL1, BT, DD, or DDK for now)
    par_file - orbital parameter file for input into PINT's photonphase
    m - number of harmonics
    ratecut - count rate cut in cts/s for cr_cut.py!
    """
    if type(obsids) != list and type(obsids) != np.ndarray:
        raise TypeError("obsids should either be a list or an array!")
    if type(merged_id) != str:
        raise TypeError("merged_id should be a string!")
    if len(merged_id) != 6:
        raise ValueError("merged_id should have 6 digits in the string!")

    merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/merged' + merged_id + '/'
    Lv2_mkdir.makedir(merged_dir)### making directory for the merged_id

    ### create a text file which lists all the cleanfilt event files with PULSE_PHASE!
    cleanfilt_files = open(merged_dir + 'merged' + merged_id + '.txt','w')
    orb_files = open(merged_dir + 'merged' + merged_id + '_orb.txt','w')
    for i in range(len(obsids)):
        filename = Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/cleanfilt_phase'
        if model == '':
            filename = filename + '.evt'
        if model != 'ELL1' and model != 'BT' and model != 'DD' and model != 'DDK':
            raise ValueError("The model should be one of ELL1, BT, DD, or DDK! Update if a new model has been added to PINT.")
        else:
            filename = filename + '_' + model + '.evt'
        cleanfilt_files.write(filename + '\n')

        orbs_filename = Lv0_dirs.NICERSOFT_DATADIR + obsids[i] + '_pipe/ni' + obsids[i] + '.orb'
        orb_files.write(orbs_filename + '\n')

    cleanfilt_files.close()
    orb_files.close()

    ##### merging the files with merge.py
    infile = '@' + merged_dir + 'merged' + merged_id + '.txt'
    outroot = 'merged' + merged_id
    outdir = merged_dir
    orb = merged_dir + 'merged' + merged_id + '_orb.txt'

    subprocess.check_call(['merge.py',infile,outroot,outdir,'--par',par_file,'--orb',orb,'--cut','--cut_rate',str(ratecut)])
    #####

    return

def get_Z2(eventfile,m):
    """
    Calculate the Z^2 significances given the event file and harmonic number m

    eventfile - name of event file
    m - number of harmonics
    """
    phases = fits.open(eventfile)[1].data['PULSE_PHASE']

    z_vals = z2m(phases,m=m)
    probs = sf_z2m(z_vals)
    significances = sig2sigma(probs)

    return significances

def niextract(eventfile,E1,E2):
    """
    Doing energy cuts only for the rate cut event file!

    eventfile - event file name
    E1 - lower energy bound (should be in a 4-digit PI string)
    E2 - upper energy bound (should be in a 4-digit PI string)
    """
    new_event = eventfile[:-4] + '_' + E1 + '-' + E2 + '.evt'
    subprocess.check_call(['niextract-events',eventfile+'[PI='+str(int(E1))+':'+str(int(E2))+']',new_event])

    return

def edit_par(par_file,var_dict):
    """
    Editing the input par file with updated parameter values in the form of a
    dictionary. Each key will have 1 value though.

    par_file - orbital parameter file for input into PINT's photonphase
    var_dict - dictionary, where each key corresponds to a variable to change in the par file, and has a 1-entry list!
    """
    dict_keys = var_dict.keys() #the variables that we want to iterate over
    new_par = par_file[:-4] + '_iter.par'

    line_no_dict = {}
    contents = open(par_file,'r').read().split('\n')
    for i in range(len(dict_keys)):
        line_no = [j for j in range(len(contents)) if dict_keys[i] in contents[j]][0] #finding indices for the variables we want to iterate over
        line_no_dict[dict_keys[i]] = [line_no] #dictionary key corresponds to a 1-entry list

    output_par = open(new_par,'w')
    for i in range(len(contents)):
        if contents[i].split(' ')[0] in dict_keys: #so if we reach a line corresponding to a variable we're iterating over
            output_par.write(contents[i].split(' ')[0] + ' '*10 + str(var_dict[contents[i].split(' ')[0]][0]) + '\n')
        else:
            output_par.write(contents[i] + '\n')
    output_par.close()

    return

def call_photonphase(obsid,model,par_file):
    """
    Calls photonphase from PINT to calculate the phase value for each event.

    obsid - Observation ID of the object of interest (10-digit str)
    model - binary model being used (ELL1, BT, DD, or DDK for now)
    par_file - orbital parameter file for input into PINT's photonphase
    """
    filename = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/cleanfilt.evt'
    orbfile = Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '.orb'

    if model == '':
        outfile_name = filename[:-4] + '_phase.evt'
    if model != 'ELL1' and model != 'BT' and model != 'DD' and model != 'DDK':
        raise ValueError("The model should be one of ELL1, BT, DD, or DDK! Update if a new model has been added to PINT.")
    else:
        outfile_name = filename[:-4] + '_phase_' + model + '.evt'

    subprocess.check_call(['photonphase',filename,par_file,'--orbfile',orbfile,'--addphase','--barytime','--outfile',outfile_name])

    return

if __name__ == "__main__":
    #print(get_Z2('1030180105','BT','',2))
    E_cut = [True,'0035','0150']
    obsids = ['2060060' + str(i) for i in range(351,372)]
    par_file = Lv0_dirs.NICERSOFT_DATADIR + 'J1231-1411_paul.par'
    #print(get_Z2_merged('000009',obsids,'ELL1',parfile,2))
    merged_id = '000012'
    source = 'J1231-1411'
    var = 'F0'
    merged_dir = Lv0_dirs.NICERSOFT_DATADIR+'/merged_events/merged' + merged_id
    m = 2 #harmonic number
    model = 'ELL1'

    #A1 = [round(i,2) for i in np.arange(2.03,2.05,0.01)]
    #PB = []
    F0 = [round(i,5) for i in np.arange(271.45295,271.45305,0.00001)]
    #TASC = []

    for i in tqdm(range(len(F0))): #for each TRIAL frequency
        if os.path.exists(merged_dir):
            shutil.rmtree(merged_dir) #so that I don't get millions of these folders!
            #If I want to get the optimized one, I can re-do them once I get the optimized set of parameters.
        edit_par(par_file,{'F0':[F0[i]]}) #change the parameter file
        for j in range(len(obsids)): #for each ObsID in the intended merge set
            call_photonphase(obsids[j],model,par_file[:-4]+'_iter.par') #call the iterated par file
        get_merged(merged_id,obsids,model,par_file[:-4]+'_iter.par',m,2.5)

        if E_cut[0] == False:
            merged_file = merged_dir + '/merged' + merged_id + '.evt'
            cut_merged_file = merged_dir + '/merged' + merged_id + '_cut.evt'
            new_merged_file = Lv0_dirs.NICERSOFT_DATADIR + source + '_' + var + '/merged' + merged_id + '_' + var + '_' + str(F0[i]) +'.evt'
            new_cut_file = Lv0_dirs.NICERSOFT_DATADIR + source + '_' + var + '/merged' + merged_id + '_cut_' + var + '_' + str(F0[i]) + '.evt'

            subprocess.check_call(['cp',merged_file,new_merged_file])
            subprocess.check_call(['cp',cut_merged_file,new_cut_file])

        else:
            original_cut = merged_dir + '/merged' + merged_id + '_cut.evt'
            cut_merged_file = original_cut[:-4] + '_' + E_cut[1] + '-' + E_cut[2] + '.evt'
            new_cut_file = Lv0_dirs.NICERSOFT_DATADIR + source + '_' + var + '/merged' + merged_id + '_cut_' + E_cut[1] + '-' + E_cut[2] + '_' + var + '_' + str(F0[i]) + '.evt'
            niextract(original_cut,E_cut[1],E_cut[2])

            subprocess.check_call(['cp',cut_merged_file,new_cut_file])

    Z2_table = open(Lv0_dirs.NICERSOFT_DATADIR + source + '_' + var + '.txt','w')
    for i in tqdm(range(len(F0))):
        if E_cut[0] == False:
            eventfile = Lv0_dirs.NICERSOFT_DATADIR + source + '_' + var + '/merged' + merged_id + '_cut_' + var + '_' + str(F0[i]) + '.evt'
        else:
            eventfile = Lv0_dirs.NICERSOFT_DATADIR + source + '_' + var + '/merged' + merged_id + '_cut_' + E_cut[1] + '-' + E_cut[2] + '_' + var + '_' + str(F0[i]) + '.evt'
        Z2_table.write('F0: ' + str(F0[i]) + ', Z^2 = ' + str(get_Z2(eventfile,m)) + '\n')
    Z2_table.close()

    #throw event files into /NICERSOFT_DATADIR/J1231-1411_F0/ or something!
