#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Saturday Jan 18 6:04pm 2020

Performing pulsation search techniques on archival data coming from NICER (or others!)

Don't need to search for burst candidates here - already have it in Lv3_TBOs.py!

Pretty similar to Lv3_quicklook.py, but takes as input a list of NICER observation directories.
Search techniques to add which are not part of the quick look scheme:
- Sideband searches, but with a way to filter out known modulations (Lv2_spurious that I was thinking of before? e.g., ISS orbit)
- Power spectral stacking
- semi-coherent searches

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import pathlib
import Lv0_dirs,Lv0_fits2dict,Lv0_scp
import Lv1_barycorr
import Lv2_preprocess,Lv2_lc,Lv2_ps,Lv2_color,Lv2_phase,Lv2_efsearch,Lv2_TBOs_method,Lv2_merging_events,Lv2_average_ps_methods,Lv2_presto_subroutines
import Lv3_E_boundary,Lv3_detection_level
from matplotlib.backends.backend_pdf import PdfPages
import time
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import subprocess
import glob
import time

Lv0_dirs.global_par() #obtaining the global parameters

start_time = time.time()

if __name__ == "__main__":
    obsdirs = [Lv0_dirs.NICER_DATADIR + '00' + str(i) +'/' for i in range(34070101,34070105)]
    nicer_obsids = [obsdirs[i][-11:-1] for i in range(len(obsdirs))]

    eventfiles = [Lv0_dirs.NICERSOFT_DATADIR + nicer_obsids[i] + '_pipe/ni' + nicer_obsids[i] + '_nicersoft_bary.evt' for i in range(len(obsdirs))]

############################# PREPROCESSING PARAMETERS ############################
    do_preprocess = False #scp, unzips, runs nicerl2, runs psrpipe, and then barycorr
    if do_preprocess == True:
        for i in range(len(nicer_obsids)):
            Lv0_scp.scp(nicer_obsids[i])

    nicerl2_flags = ['clobber=YES','ang_dist=0.035']
    psrpipe_flags = ['--emin','0.3','--emax','12.0','--angdist','0.035'] #for psrpipe in Lv0_psrpipe
    refframe = 'ICRS'
    parfile = ''

    custom_coords = np.array([])
################################################################################

############################## MERGING PARAMETERS ##############################
    do_merge = True
    merged_yet = True
    merged_id = '000014'
    if do_merge == True: #Lv2_merging_events is put higher up because it's quite important!
        if merged_yet == False: #if the ObsIDs haven't been merged yet
            Lv2_merging_events.merging(nicer_obsids)
            Lv2_merging_events.merging_GTIs(nicer_obsids,merged_id)
            #### will add in an interface with cr_cut.py at some point...!

        eventfiles = [Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/merged' + merged_id + '/merged' + merged_id + '_nicersoft_bary.evt']
################################################################################

######################### LC_PS_PHASE_COLOR PARAMETERS #########################
    do_lc_ps_phase_color = False
##### fill out the parameters below ; there are quite a number
################################################################################

############################### PRESTO PARAMETERS ##############################
    do_presto = False
    if do_presto == True:
        conv_fits2presto = False #only for process_all
        process_all = False #whether to process the ENTIRE observation or not

        preprocess_segments = True #whether to get GTI files, run niextract, edit_inf, edit_binary, nicerfits2presto
        process_segments = True #whether to make truncations to the data

        time_segments = True #truncate by time segments
        energy_segments = True #truncate by energy range
        time_energy_segments = True #truncate both by time segments and energy range
        accelsearch = True
        prepfold = True

        tbin = '1' #time bin for PRESTO in seconds

        ##### From Lv2_presto_all
        ##### For when we analyze the entire data set (no truncations in energy or time)
        accelsearch_flags = ['-numharm','4','-zmax','10','-photon','-flo','0.1','-fhi','1']

        ##### Parameters for prepfold
        zmax = 10

        ##### From Lv2_presto_subroutines
        segment_lengths = [10000] #desired length of segments in seconds
        PI1 = [30,300]
        PI2 = [300,1200]
################################################################################

############################ AVERAGE_PS PARAMETERS #############################
    do_average_ps = False
    if do_average_ps == True:
        demod = False
        preprocessing = False #getting GTIs, running niextract, nicerfits2presto, edit_inf, edit_bin, and realfft!
        time_segments = True
        time_energy_segments = True
        mode = 't'
        segment_lengths = [500] #desired length of segments in seconds
        PI1 = ['']
        PI2 = ['']
        par_file = ''
        tbin = 1 #bin size in s
        threshold = 5 #threshold for counts in each segment (in 1s bins; in terms of a percentage)
        W = 1 #number of consecutive Fourier bins to average over
        starting_freq = 0.1 #for the noise histogram
        xlims = np.array([])
        plot_mode = "save"

################################################################################

############################# SEARCH_BIN PARAMETERS ############################
##### SIDEBAND SEARCHES
    do_searchbin = False
    fft_files = [eventfiles[i][:-3] + 'fft' for i in range(len(eventfiles))]
    ncand = 1000 #no. of candidates to try to return ; default 100 (1 to 1000)
    flo = 1 #the lowest frequency to check
    fhi = 1000 #the highest frequency to check
    overlap = 0.05
    extra_args = []
################################################################################

############################# EF_SEARCH PARAMETERS #############################
##### EPOCH FOLDING

    do_efsearch = True
    if do_efsearch == True:
        n_segments = 100 #number of segments to break the epoch folding search into
        dper = 4.80 #Value for the period used in the folding. In 'efsearch' the
        #input period represents the centre of the range of the trial periods.
        nphase = 32 #Number of phases in the folded light curve(s). Typing 'INDEF'
        #forces the task to use the default value (see parameter "nbdf").
        nper = 64 #The number of periods over which the search is carried out
        dres = 1E-4 # The period resolution is the spacing between two contiguous periods in the search.
        #'INDEF' uses the default value of: half the Fourier resolution in the interval (e.g., P^2/T(i)/2 ; T(i) is interval duration)
        plot_efsearch = 'no' #to plot the results from efsearch ; do "exit" to see the next plot!

################################################################################

    if do_preprocess == True: #unzips, runs nicerl2, runs psrpipe, and then barycorr
        for i in range(len(obsdirs)):
            obsid = obsdirs[i][-11:-1]
            orbitfiles = [obsdirs[i] + '/auxil/ni' + obsid + '.orb' for i in range(len(obsdirs))]
            nicer_datafile = [obsdirs[i] + 'xti/event_cl/ni' + obsid + '_0mpu7_cl.evt' for i in range(len(obsdirs))]
            nicer_output = [obsdirs[i] + 'xti/event_cl/ni' + obsid + '_0mpu7_cl_bary.evt' for i in range(len(obsdirs))]
            nicersoft_datafile = [Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/cleanfilt.evt' for i in range(len(obsdirs))]
            nicersoft_output = [Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/ni' + obsid + '_nicersoft_bary.evt' for i in range(len(obsdirs))]
            nicersoft_folder = [Lv0_dirs.NICERSOFT_DATADIR + obsid + '_pipe/' for i in range(len(obsdirs))]

            Lv2_preprocess.preprocess(obsdirs[i],nicerl2_flags,psrpipe_flags,refframe,orbitfiles[i],parfile,nicer_datafile[i],nicer_output[i],nicersoft_datafile[i],nicersoft_output[i],nicersoft_folder[i],custom_coords)

    if do_lc_ps_phase_color == True: #would only be useful for looking at the whole time series? Might be of limited usage depending on what I want to use it for.
        par_list = ['PI','PI_FAST','TIME'] #parameter list from event_cl
        tbin_size = 1 #how you want to bin the light curve data
        Ebin_size = 0.05 #in keV
        mode = 'show' #probably best to 'save' if using a LARGE set of ObsIDs!
        truncations = 'all' #'all', 't', 'E', or 'tE', depending on whether we want to look at entire time series (all), or truncation by time interval (t), or time truncation by energy range (E), or truncation by both (tE)

        lc = True
        ps = True
        phase = True
        color = True
        ###############################################################################

        #### DEFINE DESIRED TIME INTERVALS AND ENERGY RANGES HERE FOR:
        # Lv2_ps - partial_t, partial_E, partial_tE
        # Lv2_phase - partial_t, partial_E, partial_tE
        # Lv2_color - plotting_t

        t1 = 2629580
        t2 = 2630180
        E1 = 0.3
        E2 = 2.5

        #for Lv2_ps
        ps_type = 'manual' # 'period' (for periodogram) or 'manual' (for FFT) or 'both'
        oversampling = [False,5] # [False to NOT oversample, oversampling factor - 5 to oversample by factor of 5. (factor-1) sets of 0s are padded.]
        xlims = [False,0,5] # [False to NOT impose xlimit on plots; 2nd/3rd entries are the desired x-limits if needed.]
        vlines = [True,0.2081] # [False to NOT draw a vertical line on the plot; 2nd entry is the equation for the vertical line, e.g. x=2]

        #for Lv2_phase
        ### For an unknown observation, one should run JUST Lv2_lc and Lv2_ps first to get
        ### the pulsation frequencies. Pulse profiles come LATER.
        ### If I have pulse_pars[1] and pulse_pars[2] != 0, then time binning DOES NOT MATTER, i.e., it'll be counts/s!
        pulse_pars = [0.2081,0,0]
        shift = 0.4 # how much to shift the pulse by in the phase axis. It only affects how the pulse profile is 'displaced'.
        no_phase_bins = 20 # number of phase bins desired

        #for Lv2_color
        E1_data = 0.3 #data is reliable between 0.3 and 12 keV
        E2_data = 12 # in keV
        cut_type = 'manual' # 'manual' cut for boundary energy, or 'median' - for half number of counts
        bound = 2.7 # boundary energy for when cut_type = 'manual'!

        for i in range(len(eventfiles)):
            E_bound = Lv3_E_boundary.E_bound(eventfiles[i],par_list,E1_data,E2_data,cut_type,bound) #use Lv3_E_boundary to get boundary energy
            ############################ FOR WHOLE OBSERVATION ############################
            if truncations == 'all':
                if lc == True:
                    Lv2_lc.whole(eventfiles[i],par_list,tbin_size,mode) #light curve
                    time.sleep(1)
                if ps == True:
                    Lv2_ps.whole(eventfiles[i],par_list,tbin_size,mode,ps_type,oversampling,xlims,vlines) #power spectra
                    time.sleep(1)
                if phase == True:
                    Lv2_phase.whole(eventfiles[i],par_list,tbin_size,pulse_pars,shift,no_phase_bins,mode)
                    time.sleep(1)
                if color == True:
                    Lv2_color.plotting(eventfiles[i],par_list,E_bound,tbin_size,mode)

            ########################## FOR DESIRED TIME INTERVAL ##########################
            if truncations == 't':
                if lc == True:
                    Lv2_lc.partial_t(eventfiles[i],par_list,tbin_size,t1,t2,mode) #light curve
                    time.sleep(1)
                if ps == True:
                    Lv2_ps.partial_t(eventfiles[i],par_list,tbin_size,t1,t2,mode,ps_type,oversampling,xlims,vlines) #power spectra
                    time.sleep(1)
                if phase == True:
                    Lv2_phase.partial_t(eventfiles[i],par_list,tbin_size,pulse_pars,shift,no_phase_bins,t1,t2,mode)
                    time.sleep(1)
                if color == True:
                    Lv2_color.plotting_t(eventfiles[i],par_list,E_bound,tbin_size,t1,t2,mode)

            ########################### FOR DESIRED ENERGY RANGE ##########################
            # anticipate that this will be used much?
            if truncations == 'E':
                if lc == True:
                    Lv2_lc.partial_E(eventfiles[i],par_list,tbin_size,Ebin_size,E1,E2,mode)
                    time.sleep(1)
                if ps == True:
                    Lv2_ps.partial_E(eventfiles[i],par_list,tbin_size,Ebin_size,E1,E2,mode,ps_type,oversampling,xlims,vlines)
                    time.sleep(1)
                if phase == True:
                    Lv2_phase.partial_E(eventfiles[i],par_list,tbin_size,Ebin_size,pulse_pars,shift,no_phase_bins,E1,E2,mode)

            ################# FOR DESIRED TIME INTERVAL AND ENERGY RANGE #################
            if truncations == 'tE':
                if lc == True:
                    Lv2_lc.partial_tE(eventfiles[i],par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode)
                    time.sleep(1)
                if ps == True:
                    Lv2_ps.partial_tE(eventfiles[i],par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode,ps_type,oversampling,xlims,vlines)
                    time.sleep(1)
                if phase == True:
                    Lv2_phase.partial_tE(eventfiles[i],par_list,tbin_size,Ebin_size,pulse_pars,shift,no_phase_bins,t1,t2,E1,E2,mode)
                    time.sleep(1)
                if color == True:
                    Lv2_color.plotting_t(eventfiles[i],par_list,E_bound,tbin_size,t1,t2,mode)

    if do_average_ps == True:
        for k in range(0,len(PI1),2):
            for j in range(len(segment_lengths)):
                for i in range(len(eventfiles)):
                    N = Lv3_detection_level.N_trials(tbin,segment_lengths[j])

                    if preprocessing == True:
                        if time_segments == True or time_energy_segments == True:
                            Lv2_presto_subroutines.get_gti_file(eventfiles[i],segment_lengths[j])
                        if time_segments == True:
                            Lv2_presto_subroutines.niextract_gti_time(eventfiles[i],segment_lengths[j])
                        if time_energy_segments == True:
                            Lv2_presto_subroutines.niextract_gti_time_energy(eventfiles[i],segment_lengths[j],PI1[k],PI2[k])

                        if demod == True:
                            Lv2_average_ps_methods.do_demodulate(eventfiles[i],segment_lengths[j],mode,par_file)

                        Lv2_average_ps_methods.do_nicerfits2presto(eventfiles[i],tbin,segment_lengths[j])
                        Lv2_average_ps_methods.edit_inf(eventfiles[i],tbin,segment_lengths[j])
                        Lv2_average_ps_methods.edit_binary(eventfiles[i],tbin,segment_lengths[j])
                        Lv2_average_ps_methods.realfft(eventfiles[i],segment_lengths[j])

                    Lv2_average_ps_methods.plotting(eventfiles[i],segment_lengths[j],demod,tbin,threshold,PI1[k],PI2[k],starting_freq,W,N,xlims,plot_mode)

    if do_presto == True:
        if process_all == True:
            if conv_fits2presto == True:
                for i in range(len(eventfiles)):
                    Lv2_presto_subroutines.do_nicerfits2presto(eventfiles[i],tbin,0)

            if accelsearch == True:
                print('Doing realfft/accelsearch now!')
                ### Running realfft and accelsearch from PRESTO
                for i in range(len(eventfiles)):
                    Lv2_presto_subroutines.realfft(eventfiles[i],0,'all')
                    Lv2_presto_subroutines.accelsearch(eventfiles[i],0,'all',accelsearch_flags)

            ## no_cand might be a list, if I'm processing multiple eventfiles at once...
            if prepfold == True:
                for i in range(len(eventfiles)):
                    Lv2_presto_subroutines.prepfold(eventfiles[i],0,'all',zmax)
                    Lv2_presto_subroutines.ps2pdf(eventfiles[i],0,'all')

################################################################################
############################### PRESTO_SEGMENTS ################################
################################################################################

        if process_segments == True:
            if time_segments == True and preprocess_segments == True:
                for j in range(len(segment_lengths)):
                    for i in range(len(eventfiles)):
                        Lv2_presto_subroutines.get_gti_file(eventfiles[i],segment_lengths[j]) #make GTI files for each segment
                        Lv2_presto_subroutines.niextract_gti_time(eventfiles[i],segment_lengths[j]) #performing niextract-events

                        Lv2_presto_subroutines.do_nicerfits2presto(eventfiles[i],tbin,segment_lengths[j])
                        Lv2_presto_subroutines.edit_inf(eventfiles[i],tbin,segment_lengths[j])
                        Lv2_presto_subroutines.edit_binary(eventfiles[i],tbin,segment_lengths[j])

        #                Lv3_duty_cycle.duty_cycle(eventfiles[i],tbin,segment_lengths[j],duty_cycle_bin,threshold)
        #                Lv3_duty_cycle.duty_cycle_dist(eventfiles[i],tbin,segment_lengths[j],duty_cycle_bin,threshold)

            if energy_segments == True and preprocess_segments == True:
                if len(PI1) != len(PI2):
                    raise ValueError("Make sure that the length of PI1 and PI2 are the same! Need pairs of PI values.")
                for j in range(len(PI1)):
                    for i in range(len(eventfiles)):
                        Lv2_presto_subroutines.niextract_gti_energy(eventfiles[i],PI1[j],PI2[j])
                        Lv2_presto_subroutines.do_nicerfits2presto(eventfiles[i],tbin,0) #segment_length makes no sense, so 0 is a placeholder

            if time_energy_segments == True and preprocess_segments == True:
                for j in range(len(segment_lengths)):
                    for i in range(len(eventfiles)):
                        Lv2_presto_subroutines.get_gti_file(eventfiles[i],segment_lengths[j]) #make GTI files for each segment
                        for k in range(len(PI1)):
                            Lv2_presto_subroutines.niextract_gti_time_energy(eventfiles[i],segment_lengths[j],PI1[k],PI2[k])

                        Lv2_presto_subroutines.do_nicerfits2presto(eventfiles[i],tbin,segment_lengths[j])
                        Lv2_presto_subroutines.edit_inf(eventfiles[i],tbin,segment_lengths[j])
                        Lv2_presto_subroutines.edit_binary(eventfiles[i],tbin,segment_lengths[j])

                    #for k in range(len(PI1)):
                    #    Lv3_duty_cycle.duty_cycle_tE(eventfiles[i],tbin,segment_lengths[j],PI1[k],PI2[k],duty_cycle_bin,threshold)
                    #    Lv3_duty_cycle.duty_cycle_tE_dist(eventfiles[i],tbin,segment_lengths[j],PI1[k],PI2[k],duty_cycle_bin,threshold)

            ### Running realfft and accelsearch from PRESTO
            if accelsearch == True:
                if time_segments == True or time_energy_segments == True:
                    for j in range(len(segment_lengths)):
                        for i in range(len(eventfiles)):
                            Lv2_presto_subroutines.realfft(eventfiles[i],segment_lengths[j],'t')
                            Lv2_presto_subroutines.accelsearch(eventfiles[i],segment_lengths[j],'t',accelsearch_flags)
                if energy_segments == True:
                    for i in range(len(eventfiles)):
                        Lv2_presto_subroutines.realfft(eventfiles[i],0,'E')
                        Lv2_presto_subroutines.accelsearch(eventfiles[i],0,'E',accelsearch_flags)
                else:
                    "None of time_segments, time_energy_segments, or energy_segments are True!"

            ## no_cand might be a list, if I'm processing multiple eventfiles at once...
            if prepfold == True:
                if time_segments == True or time_energy_segments == True:
                    for j in range(len(segment_lengths)):
                        for i in range(len(eventfiles)):
                            Lv2_presto_subroutines.prepfold(eventfiles[i],segment_lengths[j],'t',zmax)
                if energy_segments == True:
                    for i in range(len(eventfiles)):
                        Lv2_presto_subroutines.prepfold(eventfiles[i],0,'E',zmax)

                ### doing ps2pdf
                if time_segments == True or time_energy_segments == True:
                    for j in range(len(segment_lengths)):
                        for i in range(len(eventfiles)):
                            Lv2_presto_subroutines.ps2pdf(eventfiles[i],segment_lengths[j],'t')
                if energy_segments == True:
                    for i in range(len(eventfiles)):
                        Lv2_presto_subroutines.ps2pdf(eventfiles[i],0,'E')

                else:
                    "None of time_segments, time_energy_segments, or energy_segments are True!"

    if do_searchbin == True:
        for i in range(len(fft_files)):
            search_bin_cmd = ['search_bin','-ncand',str(ncand),'-flo',str(flo),'-fhi',str(fhi),'-overlap',str(overlap)] + extra_args + [fft_files[i]]
            subprocess.run(search_bin_cmd)

    if do_efsearch == True:
        print('Doing epoch-folding searches!')
        for i in tqdm(range(len(eventfiles))):
            eventfile_header = fits.open(eventfiles[i])[1].header
            T = eventfile_header['TSTOP'] - eventfile_header['TSTART']
            nbint = int((T/(dper/nphase))/n_segments) # The number of newbins per interval used in the analysis. The
            #"nbint" together with the NEWBIN duration determines the length in time of an interval
            #and therefore the total number of intervals within the start and stop time over which the
            #analysis will be carried out. Typing 'INDEF' forces the task to use the default value
            #(see parameter "nbdf"). NOTE: By pressing return "nbint" is set to the value found in the
            #parameter file used in a previous run."
            outfile_root = nicer_obsids[i] + '_' + str(n_segments) + 'segs_' + str(nper)
            Lv2_efsearch.efsearch(eventfiles[i],n_segments,dper,nphase,nbint,nper,dres,outfile_root,plot_efsearch)

end_time = time.time()

print('Time elapsed: ' + str(end_time-start_time) + ' seconds.')
