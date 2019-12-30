#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 5:18pm 2019

Script that takes in an input (barycentered) merged time series from multiple
data_ids, and outputs a pulse profile! Methods script.

"""
from __future__ import division, print_function
import numpy as np
import subprocess
from astropy.io import fits

import binary_psr
import Lv0_dirs,Lv2_phase
import matplotlib.pyplot as plt

def niextract_gti_energy(merging,data_id,PI1,PI2):
    """
    Using niextract-events to get segmented data based on the energy range

    merging - True/False - whether to use merged data
    data_id - 10-digit ObsID or 6-digit ID for the merged event file
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    if type(data_id) != str:
        raise TypeError("data_id should be a string!")
    if len(data_id) != 6 and len(data_id) != 10:
        raise ValueError("data_id should have 6 or 10 digits in the string!")

    if merging == True:
        print('Extracting energy segments, merging is True!')
        all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
        merged_dir = all_merged_dir + 'merged' + data_id + '/'
        merged_filename = merged_dir + 'merged' + data_id + '_nicersoft_bary.evt'
        output_file = merged_dir + 'merged' + data_id + '_nicersoft_bary_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '.evt'

        subprocess.check_call(['niextract-events',merged_filename + '[PI='+str(PI1)+':'+str(PI2)+']',output_file])

    else:
        print('Extracting energy segments, merging is False!')
        event = Lv0_dirs.NICERSOFT_DATADIR + data_id + '_pipe/ni' + data_id + '_nicersoft_bary.evt'

        working_dir = Lv0_dirs.NICERSOFT_DATADIR+data_id+'_pipe/'
        inputfile = working_dir + 'ni' + data_id + '_nicersoft_bary.evt'
        outputfile = working_dir + 'ni' + data_id + '_nicersoft_bary_E'+str(PI1).zfill(4)+'-'+str(PI2).zfill(4)+'.evt'

        subprocess.check_call(['niextract-events',inputfile+'[PI='+str(PI1)+':'+str(PI2)+']',outputfile])


    return

def do_demodulate(merging,data_id,par_file,E_trunc,PI1,PI2):
    """
    Using do_demodulate in binary_psr.py in Scott Ransom's PRESTO Python library to
    demodulate the time series for the merged event file!

    merging - True/False - whether to use merged data
    data_id - 10-digit ObsID or 6-digit ID for the merged event file
    par_file - orbital parameter file for input into binary_psr
    E_trunc - True/False - whether to do energy truncation
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    """
    if type(data_id) != str:
        raise TypeError("data_id should be a string!")
    if len(data_id) != 6 and len(data_id) != 10:
        raise ValueError("data_id should have 6 or 10 digits in the string!")
    if E_trunc != True and E_trunc != False:
        raise ValueError("E_trunc (energy truncation) should either be True or False!")

    if merging == True:
        print("Doing demodulation, merging is True!")
        all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
        merged_dir = all_merged_dir + 'merged' + data_id + '/'
        if E_trunc == True:
            old_file = merged_dir + 'merged' + data_id + '_nicersoft_bary_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '.evt'
        else:
            old_file = merged_dir + 'merged' + data_id + '_nicersoft_bary.evt'

    else:
        print("Doing demodulation, merging is False!")
        if E_trunc == True:
            old_file = Lv0_dirs.NICERSOFT_DATADIR + data_id + '_pipe/ni' + data_id + '_nicersoft_bary_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '.evt'
        else:
            old_file = Lv0_dirs.NICERSOFT_DATADIR + data_id + '_pipe/ni' + data_id + '_nicersoft_bary.evt'

    new_file = old_file[:-4] + '_demod.evt'
    subprocess.check_call(['cp',old_file,new_file])
    with fits.open(new_file,mode='update') as fitsfile_demod:
        MJDREFI = fitsfile_demod[1].header['MJDREFI'] #integer for MJD reference
        MJDREFF = fitsfile_demod[1].header['MJDREFF'] #float decimal for MJD reference

        times = fitsfile_demod[1].data['TIME'] #original time series
        gtis_start = fitsfile_demod[2].data['START'] #original GTI start times
        gtis_stop = fitsfile_demod[2].data['STOP'] #original GTI end times

        times_MJD = MJDREFI + MJDREFF + times/86400 #converting METs to MJD
        gtis_start_MJD = MJDREFI + MJDREFF + gtis_start/86400 #converting GTIs in METs to MJD
        gtis_stop_MJD = MJDREFI + MJDREFF + gtis_stop/86400 #converting GTIs in METs to MJD

        try:
            times_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(times_MJD) #demodulated event times
            gtis_start_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(gtis_start_MJD) #demodulated GTI start times
            gtis_stop_demod = binary_psr.binary_psr(par_file).demodulate_TOAs(gtis_stop_MJD) #demodulated GTI end times

            fitsfile_demod[1].data['TIME'] = (times_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs
            fitsfile_demod[2].data['START'] = (gtis_start_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs
            fitsfile_demod[2].data['STOP'] = (gtis_stop_demod - MJDREFI - MJDREFF) * 86400 #convert back to METs

            fitsfile_demod.flush()

        except ValueError:
            pass
    return

def pulse_profile(merging,data_id,E_trunc,PI1,PI2,pulse_pars,no_phase_bins):
    """
    Extracts the time series from the demodulated merged event file, and creates
    the pulse profile from Lv2_phase.py!

    merging - True/False - whether to use merged data
    data_id - 10-digit ObsID or 6-digit ID for the merged event file
    E_trunc - True/False - whether to do energy truncation
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    pulse_pars - parameters corresponding to the pulse
    no_phase_bins - number of phase bins desired
    """
    if type(data_id) != str:
        raise TypeError("data_id should be a string!")
    if len(data_id) != 6 and len(data_id) != 10:
        raise ValueError("data_id should have 6 or 10 digits in the string!")
    if E_trunc != True and E_trunc != False:
        raise ValueError("E_trunc (energy truncation) should either be True or False!")
    if type(pulse_pars) != list and type(pulse_pars) != np.ndarray:
        raise TypeError("pulse_pars should either be a list or an array!")

    if merging == True:
        print("Making pulse profile, merging is True!")
        all_merged_dir = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/' #directory for the merged events
        merged_dir = all_merged_dir + 'merged' + data_id + '/'
        if E_trunc == True:
            demod_merged_file = merged_dir + 'merged' + data_id + '_nicersoft_bary_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '_demod.evt'
        else:
            demod_merged_file = merged_dir + 'merged' + data_id + '_nicersoft_bary_demod.evt'

    else:
        print("Making pulse profile, merging is False!")
        if E_trunc == True:
            demod_merged_file = Lv0_dirs.NICERSOFT_DATADIR + data_id + '_pipe/ni' + data_id + '_nicersoft_bary_E' + str(PI1).zfill(4) + '-' + str(PI2).zfill(4) + '_demod.evt'
        else:
            demod_merged_file = Lv0_dirs.NICERSOFT_DATADIR + data_id + '_pipe/ni' + data_id + '_nicersoft_bary_demod.evt'

    #demod_merged_file = demod_merged_file[:-10] + '.evt'
    print('Will be using ' + str(demod_merged_file))
    event = fits.open(demod_merged_file)
    times = event[1].data['TIME']

    gtis = event[2].data
    T = sum([ (gtis[i][1]-gtis[i][0]) for i in range(len(gtis)) ])
    print('Will be using ' + str(T/1000) + 'ks worth of data!')

    phase_bins,summed_profile = Lv2_phase.pulse_folding(times,T,times[0],pulse_pars[0],pulse_pars[1],pulse_pars[2],no_phase_bins)

    return phase_bins,summed_profile

def plot_pf(merging,data_id,E_trunc,PI1,PI2,pulse_pars,no_phase_bins):
    """
    Extracts the time series from the demodulated merged event file, and creates
    the pulse profile from Lv2_phase.py!

    merging - True/False - whether to use merged data
    data_id - 10-digit ObsID or 6-digit ID for the merged event file
    E_trunc - True/False - whether to do energy truncation
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    pulse_pars - parameters corresponding to the pulse
    no_phase_bins - number of phase bins desired
    """
    if type(data_id) != str:
        raise TypeError("data_id should be a string!")
    if len(data_id) != 6 and len(data_id) != 10:
        raise ValueError("data_id should have 6 or 10 digits in the string!")
    if E_trunc != True and E_trunc != False:
        raise ValueError("E_trunc (energy truncation) should either be True or False!")
    if type(pulse_pars) != list and type(pulse_pars) != np.ndarray:
        raise TypeError("pulse_pars should either be a list or an array!")

    phase_bins,summed_profile = pulse_profile(merging,data_id,E_trunc,PI1,PI2,pulse_pars,no_phase_bins)

    plt.figure()
    plt.step(phase_bins[:-1],summed_profile)
    plt.xlabel('Phase', fontsize=12)
    plt.ylabel('Count/s',fontsize=12)
    plt.show()

    return

if __name__ == "__main__":
    merging = False
    data_id = '1013010105'
    E_trunc = False
    PI1 = 30
    PI2 = 1200
    pulse_pars = [29.639575,-3.77535E-10,1.1147E-20]
    no_phase_bins = 50
    plot_pf(merging,data_id,E_trunc,PI1,PI2,pulse_pars,no_phase_bins)
