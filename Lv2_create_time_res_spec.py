#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun April 5 5:32pm 2020

Program that splits up an event file to do time-resolved spectroscopy

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv1_data_gtis,Lv2_presto_subroutines,Lv2_mkdir
import os
from tqdm import tqdm
import subprocess
import pathlib
import glob

Lv0_dirs.global_par()

def niextract_gti(eventfile,gap,gtifile,min_exp):
    """
    Using niextract-events to get segmented data based on the individual GTIs created with
    GTI_bunching in Lv1_data_gtis. (Very similar to that in Lv2_presto_subroutines,
    except I create a list of paths to the event files.)

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    gap - maximum separation between end time of 1st GTI and start time of 2nd GTI allowed
    gtifile - name of GTI file
    min_exp - minimum exposure to be used for the spectra
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    Lv1_data_gtis.GTI_bunching(eventfile,gap,gtifile)

    gtis = list(fits.open(parent_folder+'/'+gtifile)[1].data)
    niextract_folder = parent_folder + '/accelsearch_GTIs/'
    Lv2_mkdir.makedir(niextract_folder)
    for i in tqdm(range(len(gtis))):
        gti_start = gtis[i][0]
        gti_end = gtis[i][1]
        if gti_end - gti_start >= min_exp:
            subprocess.run(['mkgti.py','--gtiname',niextract_folder+'GTI'+str(i+1)+'.gti',str(gti_start),str(gti_end)])
            subprocess.run(['niextract-events',eventfile,niextract_folder+str(pathlib.Path(eventfile).name)[:-4]+'_GTI'+str(i+1).zfill(4)+'.evt','timefile='+niextract_folder+'GTI'+str(i+1)+'.gti'])

    return

def instructions(eventfile):
    """
    Writing a set of instructions to use in XSELECT, in order to extract the spectra and all!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    niextract_folder = parent_folder + '/accelsearch_GTIs/spectra/'

    instruct_file = niextract_folder + '/instructions.txt'

    eventfiles = sorted(glob.glob(niextract_folder+'/*E0050-1200.evt'))
    instruct = open(instruct_file,'w')
    for i in range(len(eventfiles)):
        eventfile = str(pathlib.Path(eventfiles[i]).name)

        instruct.write('set mission nicer' + '\n')
        instruct.write('set inst xti' + '\n')
        instruct.write('set datadir ' + niextract_folder + '\n')
        instruct.write('read event ' + eventfile + '\n')
        instruct.write('extract spectrum' + '\n')
        instruct.write('save spectrum ' + eventfile[:-3] + 'pha' + '\n')
        instruct.write('clear data' + '\n')
        instruct.write(' ' + '\n')

    instruct.close()

    return

def set_rmf_arf(eventfile):
    """
    Setting paths to response and arf files in the individual spectra files.
    To check, can do "fkeyprint FILENAME ANCRFILE/RESPFILE"

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    niextract_folder = parent_folder + '/accelsearch_GTIs/spectra'

    spectrafiles = sorted(glob.glob(niextract_folder+'/*.pha'))

    resp = '/Volumes/Samsung_T5/nicer-rmf6s-teamonly-array52.rmf'
    arf = '/Volumes/Samsung_T5/nicer-consim135p-teamonly-array52.arf'

    for i in range(len(spectrafiles)):
        subprocess.run(['fparkey',resp,spectrafiles[i],'RESPFILE'])
        subprocess.run(['fparkey',arf,spectrafiles[i],'ANCRFILE'])

    return

def grppha(eventfile):
    """
    Function that does GRPPHA on a set of pha files.
    The function will output pha files of the format 'grp_$file'!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    niextract_folder = parent_folder + '/accelsearch_GTIs/spectra/'
    binned_phas = sorted(glob.glob(niextract_folder+'*pha'))

    command_file = niextract_folder + 'grppha_commands.go'
    writing = open(command_file,'w')
    #### now to build the grppha command
    grppha = 'grppha'
    chatter = 'chatter=0'
    for i in range(len(binned_phas)):
        #comm = 'comm="group nicer_channels_to_group.dat & systematics 30-1200 0.02 & chkey ANCRFILE /Volumes/Samsung_T5/nicer-consim135p-teamonly-array50.arf & chkey RESPFILE /Volumes/Samsung_T5/nicer-rmf6s-teamonly-array50.rmf & chkey BACKFILE ' + backfile + ' & exit"'
        comm = 'comm="group nicer_channels_to_group.dat & systematics 30-1200 0.02 & exit"'
        infile = 'infile="' + str(pathlib.Path(binned_phas[i]).name) + '"'
        outfile = 'outfile="grp_' + str(pathlib.Path(binned_phas[i]).name) + '"'
        writing.write(grppha+' '+infile+' '+outfile+' '+chatter+' '+comm+'\n')
    writing.close()

def xspec_read_all(eventfile):
    """
    To read all the spectral files (with rmf,arf already set!) and ignore 0.0-0.3, 12-higher keV

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    niextract_folder = parent_folder + '/accelsearch_GTIs/spectra/'

    spectrafiles = sorted(glob.glob(niextract_folder+'/*.pha'))

    readspectra = open(niextract_folder + 'readspectra.xcm','w')
    data_string = 'data ' + ' '.join([str(i+1)+':'+str(i+1) + ' ' + spectrafiles[i] for i in range(len(spectrafiles))])
    readspectra.write(data_string + '\n')
    readspectra.write('ignore **:0.0-0.485,12.01-**' + '\n')
    readspectra.write('abund wilm' + '\n')
    readspectra.close()

    return

if __name__ == "__main__":
    eventfile = Lv0_dirs.NICER_DATADIR + 'at2019wey/2008041401_filt_bary.evt'

    gap = 5
    gtifile = 'bunched.gti'
    min_exp = 50

    #niextract_gti(eventfile,gap,gtifile,min_exp)
    #instructions(eventfile)

    #set_rmf_arf(eventfile)

    #import pathlib
    #specs = sorted(glob.glob('/Volumes/Samsung_T5/NICER-data/at2019wey/accelsearch_GTIs/spectra/*evt'))
    #for i in range(len(specs)):
    #    print(str(pathlib.Path(specs[i]).name) + ' ' + str(fits.open(specs[i])[1].header['NAXIS2']))
    #grppha(eventfile)
    xspec_read_all(eventfile)
