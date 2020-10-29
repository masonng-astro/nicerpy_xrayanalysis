#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Friday Jan 10 11:01am 2020

Performing quick look routines on fresh data coming from NICER (or others!)

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import pathlib
import Lv0_dirs,Lv0_fits2dict
import Lv1_barycorr,Lv1_data_bin
import Lv2_lc,Lv2_ps,Lv2_color,Lv2_phase,Lv2_efsearch,Lv2_TBOs_method,Lv2_presto_subroutines
import Lv3_E_boundary
from matplotlib.backends.backend_pdf import PdfPages
import time
import matplotlib.pyplot as plt
import os
import subprocess
import glob
import time

Lv0_dirs.global_par() #obtaining the global parameters

start_time = time.time()

def nicerql(eventfile,extra_nicerql_args):
    """
    Probably the second step in the process, but this is to generate the psrpipe
    diagnostic plots, to see if there are any obvious red flags in the data.

    Will just really need eventfile ; orb file and mkf file is assumed to be in the SAME folder

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    extra_nicerql_args - extra arguments to the nicerql script
    """
    parent_folder = str(pathlib.Path(eventfile).parent)
    orbfiles = glob.glob(parent_folder+'/*.orb')
    mkffiles = glob.glob(parent_folder+'/*.mkf*')
    if len(orbfiles) != 1 and len(mkffiles) != 1:
        raise ValueError("Either there's no orb/mkf file (in which case, make sure they're in the same folder as the event file!) or there are multiple files - do check!")

    logfile = parent_folder + '/nicerql.log'
    print('Running nicerql now for ' + eventfile + ':')
    with open(logfile,'w') as logtextfile:
        output = subprocess.run(['nicerql.py',eventfile,'--orb',orbfiles[0],'--mkf',mkffiles[0]] + extra_nicerql_args,capture_output=True,text=True)
        logtextfile.write(output.stdout)
        logtextfile.write('*------------------------------* \n')
        logtextfile.write(output.stderr)
        logtextfile.close()

    png_files = glob.glob('*evt*png')
    for i in range(len(png_files)):
        subprocess.run(['mv',png_files[i],parent_folder+'/'])

    return

def filtering(eventfile,outfile,maskdet,eventflags):
    """
    Function that will filter out bad detectors and impose eventflag restrictions.
    Will expect to add to this function over time...

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    outfile - output file from the filtering
    maskdet - list of DET_IDs to mask
    eventflags - NICER event flags to filter out (e.g., in the format '(EVENT_FLAGS==bx1x000)')
    """
    evfilt_expr = eventflags
    for i in range(len(maskdet)):
        evfilt_expr += '.and.(DET_ID!='+str(maskdet[i]) + ')'

    subprocess.run(['ftcopy',eventfile+'['+evfilt_expr+']',outfile,'clobber=yes','history=yes'])

if __name__ == "__main__":
    #eventfile = '/Volumes/Samsung_T5/NICER-data/swiftj0044/2004110232.evt'
    eventfile = '/Volumes/Samsung_T5/NICER-data/at2020wey/2008041401.evt'
    #eventfile = '/Users/masonng/Documents/MIT/Research/xtej1739_burst8/2002251704.evt'

############################# DIAGNOSTIC PARAMETERS ############################
    do_diagnostic = False
    extra_nicerql_args = ['--save']
################################################################################

############################# FILTERING PARAMETERS #############################
    do_filter = False
    filtfile = '/Volumes/Samsung_T5/NICER-data/xtej1739_burst19+/2002131540_filt.evt'
    maskdets = [34]
    eventflags = '(EVENT_FLAGS==bx1x000)'
################################################################################

############################# BARYCORR PARAMETERS ##############################
    do_bary = True
    out_baryfile = Lv0_dirs.NICER_DATADIR + 'at2020wey/2008041401.evt'
    #out_baryfile = '/Users/masonng/Documents/MIT/Research/xtej1739_burst8/2002251704_bary.evt'
    refframe = 'ICRS'
    orbitfile = Lv0_dirs.NICER_DATADIR + 'at2020wey/2008041401.orb'
    #orbitfile = '/Users/masonng/Documents/MIT/Research/xtej1739_burst8/2002251704.orb'
    parfile = ''
    output_folder = Lv0_dirs.NICER_DATADIR + 'at2020wey/'
    #output_folder = '/Users/masonng/Documents/MIT/Research/xtej1739_burst8/'
    custom_coords = np.array([])
################################################################################

########################### BURST_CANDS PARAMETERS #############################
    do_burst_cands = False
################################################################################

######################### LC_PS_PHASE_COLOR PARAMETERS #########################
    do_lc_ps_phase_color = False
##### fill out the parameters below ; there are quite a number
################################################################################

############################### PRESTO PARAMETERS ##############################
    do_presto = False
    conv_fits2presto = True #only for process_all
    process_all = True #whether to process the ENTIRE observation or not

    process_segments = True #whether to make truncations to the data
    time_segments = False #truncate by time segments
    gti_segments = False #truncate by GTI segments
    energy_segments = True #truncate by energy range
    time_energy_segments = False #truncate both by time segments and energy range

    accelsearch = True
    prepfold = True

    ### for gti_segments
    gap = 5 #seconds
    gtifile = 'bunched.gti'

    tbin = '0.01' #time bin for PRESTO in seconds

    ##### From Lv2_presto_all
    ##### For when we analyze the entire data set (no truncations in energy or time)
    accelsearch_flags = ['-numharm','8','-zmax','100','-photon','-flo','1','-fhi','2000']

    ##### Parameters for prepfold
    zmax = 100

    ##### From Lv2_presto_subroutines
    segment_lengths = [10000] #desired length of segments (and 200)
    ## NOTE: nicerfits2presto tries to find a 'nice number' of bins to bin the data, and
    ## sometimes that will be longer than the segment. I'm not sure what the best remedy is,
    ## but currently, the best thing is just to track the segment size on the terminal, then increase
    ## accordingly. For example, I tried 100s segments for 0034070101 (Cen X-3), but nicerfits2presto.py
    ## is binning at 100.8s, so I went up to the next nice number - 102..
    PI1 = [30,200]
    PI2 = [200,1200]
################################################################################

############################# SEARCH_BIN PARAMETERS ############################
##### SIDEBAND SEARCHES
    do_searchbin = False
    fft_file = out_baryfile[:-3] + 'fft'
    ncand = 1000 #no. of candidates to try to return ; default 100 (1 to 1000)
    flo = 1 #the lowest frequency to check
    fhi = 1000 #the highest frequency to check
    overlap = 0.05
    extra_args = []
################################################################################

############################# EF_SEARCH PARAMETERS #############################
##### EPOCH FOLDING

    do_efsearch = False
    if do_efsearch == True:
        eventfile_header = fits.open(out_baryfile)[1].header
        T = eventfile_header['TSTOP'] - eventfile_header['TSTART']
        n_segments = 1 #number of segments to break the epoch folding search into
        dper = 1.36 #Value for the period  used  in the folding. In 'efsearch' the
        #input period represents the centre of the range of the trial periods.
        nphase = 16 #Number of phases in the folded light curve(s). Typing 'INDEF'
        #forces the task to use the default value (see parameter "nbdf").
        nbint = int((T/(dper/nphase))/n_segments) # The number of newbins per interval used in the analysis. The
        #"nbint" together with the NEWBIN duration determines the length in time of an interval
        #and therefore the total number of intervals within the start and stop time over which the
        #analysis will be carried out. Typing 'INDEF' forces the task to use the default value
        #(see parameter "nbdf"). NOTE: By pressing return "nbint" is set to the value found in the
        #parameter file used in a previous run."
        nper = 128 #The number of periods over which the search is carried out
        dres = 1E-4 # The period resolution is the spacing between two contiguous periods in the search.
        #'INDEF' uses the default value of: half the Fourier resolution in the interval (e.g., P^2/T(i)/2 ; T(i) is interval duration)
        outfile_root = 'swift1818_' + str(n_segments) + 'segs_' + str(nper)
        plot_efsearch = 'no' #to plot the results from efsearch ; do "exit" to see the next plot!

################################################################################

    ##### For Lv3_duty_cycle:
    duty_cycle_bin = 1 #for 1s bins to do duty cycle calculations
    threshold = 20 #for 10%

    if do_diagnostic == True:
        ##### Quick look at diagnostic plots
        nicerql(eventfile,extra_nicerql_args)

    if do_filter == True:
        ##### Filtering if needed
        filtering(eventfile,filtfile,maskdets,eventflags)

    if do_bary == True:
        if do_filter == True:
            Lv1_barycorr.barycorr(filtfile,out_baryfile,refframe,orbitfile,parfile,output_folder,custom_coords)
        else:
            Lv1_barycorr.barycorr(eventfile,out_baryfile,refframe,orbitfile,parfile,output_folder,custom_coords)

    if do_burst_cands == True:
        Lv2_TBOs_method.burst_cands(out_baryfile)

    if do_lc_ps_phase_color == True:
        par_list = ['PI','PI_FAST','TIME'] #parameter list from event_cl
        tbin_size = 0.001 #how you want to bin the light curve data
        Ebin_size = 0.05 #in keV
        mode = 'show'
        truncations = 't' #'all', 't', 'E', or 'tE', depending on whether we want to look at entire time series (all), or truncation by time interval (t), or time truncation by energy range (E), or truncation by both (tE)

        lc = True
        ps = True
        phase = False
        color = False
        ###############################################################################

        #### DEFINE DESIRED TIME INTERVALS AND ENERGY RANGES HERE FOR:
        # Lv2_ps - partial_t, partial_E, partial_tE
        # Lv2_phase - partial_t, partial_E, partial_tE
        # Lv2_color - plotting_t

        t1 = 5450
        t2 = 6300
        E1 = 2
        E2 = 12

        #for Lv2_ps
        ps_type = 'manual' # 'period' (for periodogram) or 'manual' (for FFT) or 'both'
        oversampling = [False,5] # [False to NOT oversample, oversampling factor - 5 to oversample by factor of 5. (factor-1) sets of 0s are padded.]
        xlims = [False,0,1] # [False to NOT impose xlimit on plots; 2nd/3rd entries are the desired x-limits if needed.]
        vlines = [False,1/1.3632] # [False to NOT draw a vertical line on the plot; 2nd entry is the equation for the vertical line, e.g. x=2]

        #for Lv2_phase
        ### For an unknown observation, one should run JUST Lv2_lc and Lv2_ps first to get
        ### the pulsation frequencies. Pulse profiles come LATER.
        ### If I have pulse_pars[1] and pulse_pars[2] != 0, then time binning DOES NOT MATTER, i.e., it'll be counts/s!
        pulse_pars = [1/1.3632,0,0] #J0209.6-7427
        shift = 0.4 # how much to shift the pulse by in the phase axis. It only affects how the pulse profile is 'displaced'.
        no_phase_bins = 20 # number of phase bins desired

        #for Lv2_color
        E1_data = 0.3 #data is reliable between 0.3 and 12 keV
        E2_data = 12 # in keV
        cut_type = 'manual' # 'manual' cut for boundary energy, or 'median' - for half number of counts
        bound = 2.7 # boundary energy for when cut_type = 'manual'!

        E_bound = Lv3_E_boundary.E_bound(out_baryfile,par_list,E1_data,E2_data,cut_type,bound) #use Lv3_E_boundary to get boundary energy
        ############################ FOR WHOLE OBSERVATION ############################
        if truncations == 'all':
            if lc == True:
                Lv2_lc.whole(out_baryfile,par_list,tbin_size,mode) #light curve
                time.sleep(1)
            if ps == True:
                Lv2_ps.whole(out_baryfile,par_list,tbin_size,mode,ps_type,oversampling,xlims,vlines) #power spectra
                time.sleep(1)
            if phase == True:
                Lv2_phase.whole(out_baryfile,par_list,tbin_size,pulse_pars,shift,no_phase_bins,mode)
                time.sleep(1)
            if color == True:
                Lv2_color.plotting(out_baryfile,par_list,E_bound,tbin_size,mode)

        ########################## FOR DESIRED TIME INTERVAL ##########################
        if truncations == 't':
            if lc == True:
                Lv2_lc.partial_t(out_baryfile,par_list,tbin_size,t1,t2,mode) #light curve
                time.sleep(1)
            if ps == True:
                Lv2_ps.partial_t(out_baryfile,par_list,tbin_size,t1,t2,mode,ps_type,oversampling,xlims,vlines) #power spectra
                time.sleep(1)
            if phase == True:
                Lv2_phase.partial_t(out_baryfile,par_list,tbin_size,pulse_pars,shift,no_phase_bins,t1,t2,mode)
                time.sleep(1)
            if color == True:
                Lv2_color.plotting_t(out_baryfile,par_list,E_bound,tbin_size,t1,t2,mode)

        ########################### FOR DESIRED ENERGY RANGE ##########################
        # anticipate that this will be used much?
        if truncations == 'E':
            if lc == True:
                Lv2_lc.partial_E(out_baryfile,par_list,tbin_size,Ebin_size,E1,E2,mode)
                time.sleep(1)
            if ps == True:
                Lv2_ps.partial_E(out_baryfile,par_list,tbin_size,Ebin_size,E1,E2,mode,ps_type,oversampling,xlims,vlines)
                time.sleep(1)
            if phase == True:
                Lv2_phase.partial_E(out_baryfile,par_list,tbin_size,Ebin_size,pulse_pars,shift,no_phase_bins,E1,E2,mode)

        ################# FOR DESIRED TIME INTERVAL AND ENERGY RANGE #################
        if truncations == 'tE':
            if lc == True:
                Lv2_lc.partial_tE(out_baryfile,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode)
                time.sleep(1)
            if ps == True:
                Lv2_ps.partial_tE(out_baryfile,par_list,tbin_size,Ebin_size,t1,t2,E1,E2,mode,ps_type,oversampling,xlims,vlines)
                time.sleep(1)
            if phase == True:
                Lv2_phase.partial_tE(out_baryfile,par_list,tbin_size,Ebin_size,pulse_pars,shift,no_phase_bins,t1,t2,E1,E2,mode)
                time.sleep(1)
            if color == True:
                Lv2_color.plotting_t(out_baryfile,par_list,E_bound,tbin_size,t1,t2,mode)

    if do_presto == True:
        if process_all == True:
            if conv_fits2presto == True:
                Lv2_presto_subroutines.do_nicerfits2presto(out_baryfile,tbin,0,'all')

            if accelsearch == True:
                print('Doing realfft/accelsearch now!')
                ### Running realfft and accelsearch from PRESTO
                Lv2_presto_subroutines.realfft(out_baryfile,0,'all')
                Lv2_presto_subroutines.accelsearch(out_baryfile,0,'all',accelsearch_flags)

            ## no_cand might be a list, if I'm processing multiple out_baryfiles at once...
            if prepfold == True:
                print("Running prepfold now!")
                Lv2_presto_subroutines.prepfold(out_baryfile,0,'all',zmax)
                Lv2_presto_subroutines.ps2pdf(out_baryfile,0,'all')

################################################################################
############################### PRESTO_SEGMENTS ################################
################################################################################

        if process_segments == True:
            if time_segments == True:
                for j in range(len(segment_lengths)):
                    Lv2_presto_subroutines.get_gti_file(out_baryfile,segment_lengths[j]) #make GTI files for each segment
                    Lv2_presto_subroutines.niextract_gti_time(out_baryfile,segment_lengths[j]) #performing niextract-events

                    Lv2_presto_subroutines.do_nicerfits2presto(out_baryfile,tbin,segment_lengths[j],'t')
                    Lv2_presto_subroutines.edit_inf(out_baryfile,tbin,segment_lengths[j])
                    Lv2_presto_subroutines.edit_binary(out_baryfile,tbin,segment_lengths[j])
    #
    #                Lv3_duty_cycle.duty_cycle(out_baryfile,tbin,segment_lengths[j],duty_cycle_bin,threshold)
    #                Lv3_duty_cycle.duty_cycle_dist(out_baryfile,tbin,segment_lengths[j],duty_cycle_bin,threshold)

            if gti_segments == True:
                Lv2_presto_subroutines.niextract_gti(out_baryfile,gap,gtifile)
                Lv2_presto_subroutines.do_nicerfits2presto(out_baryfile,tbin,0,'gtis')

            if energy_segments == True:
                print("Doing acceleration searches with energy cuts!")
                if len(PI1) != len(PI2):
                    raise ValueError("Make sure that the length of PI1 and PI2 are the same! Need pairs of PI values.")
                parent_folder = str(pathlib.Path(out_baryfile).parent)
                accelE_folder = parent_folder + '/accelsearch_E/'
                file_suffix = str(pathlib.Path(out_baryfile).name)

                for j in range(len(PI1)):
                    Lv2_presto_subroutines.niextract_gti_energy(out_baryfile,PI1[j],PI2[j])
                    Lv2_presto_subroutines.do_nicerfits2presto(accelE_folder+file_suffix[:-4]+'_E'+str(PI1[j])+'-'+str(PI2[j])+'.evt',tbin,0,'E') #segment_length makes no sense, so 0 is a placeholder

            if time_energy_segments == True:
                for j in range(len(segment_lengths)):
                    Lv2_presto_subroutines.get_gti_file(out_baryfile,segment_lengths[j]) #make GTI files for each segment
                    for k in range(len(PI1)):
                        Lv2_presto_subroutines.niextract_gti_time_energy(out_baryfile,segment_lengths[j],PI1[k],PI2[k])

                    Lv2_presto_subroutines.do_nicerfits2presto(out_baryfile,tbin,segment_lengths[j],'t')
                    Lv2_presto_subroutines.edit_inf(out_baryfile,tbin,segment_lengths[j])
                    Lv2_presto_subroutines.edit_binary(out_baryfile,tbin,segment_lengths[j])

                    #for k in range(len(PI1)):
                    #    Lv3_duty_cycle.duty_cycle_tE(out_baryfile,tbin,segment_lengths[j],PI1[k],PI2[k],duty_cycle_bin,threshold)
                    #    Lv3_duty_cycle.duty_cycle_tE_dist(out_baryfile,tbin,segment_lengths[j],PI1[k],PI2[k],duty_cycle_bin,threshold)

            ### Running realfft and accelsearch from PRESTO
            if accelsearch == True:
                print("Running acceleration searches now!")
                if time_segments == True or time_energy_segments == True:
                    for j in range(len(segment_lengths)):
                        Lv2_presto_subroutines.realfft(out_baryfile,segment_lengths[j],'t')
                        Lv2_presto_subroutines.accelsearch(out_baryfile,segment_lengths[j],'t',accelsearch_flags)
                if gti_segments == True:
                    Lv2_presto_subroutines.realfft(out_baryfile,0,'gtis')
                    Lv2_presto_subroutines.accelsearch(out_baryfile,0,'gtis',accelsearch_flags)
                if energy_segments == True:
                    Lv2_presto_subroutines.realfft(out_baryfile,0,'E')
                    Lv2_presto_subroutines.accelsearch(out_baryfile,0,'E',accelsearch_flags)
                else:
                    "None of time_segments, gti_segments, time_energy_segments, or energy_segments are True!"

            ## no_cand might be a list, if I'm processing multiple eventfiles at once...
            if prepfold == True:
                print("Running prepfold now!")
                if time_segments == True or time_energy_segments == True:
                    for j in range(len(segment_lengths)):
                        Lv2_presto_subroutines.prepfold(out_baryfile,segment_lengths[j],'t',zmax)
                if gti_segments == True:
                    Lv2_presto_subroutines.prepfold(out_baryfile,0,'gtis',zmax)
                if energy_segments == True:
                    Lv2_presto_subroutines.prepfold(out_baryfile,0,'E',zmax)

                ### doing ps2pdf
                if time_segments == True or time_energy_segments == True:
                    for j in range(len(segment_lengths)):
                        Lv2_presto_subroutines.ps2pdf(out_baryfile,segment_lengths[j],'t')
                if gti_segments == True:
                    Lv2_presto_subroutines.ps2pdf(out_baryfile,0,'gtis')
                if energy_segments == True:
                    Lv2_presto_subroutines.ps2pdf(out_baryfile,0,'E')

                else:
                    "None of time_segments, time_energy_segments, or energy_segments are True!"

    if do_searchbin == True:
        search_bin_cmd = ['search_bin','-ncand',str(ncand),'-flo',str(flo),'-fhi',str(fhi),'-overlap',str(overlap)] + extra_args + [fft_file]
        subprocess.run(search_bin_cmd)

    if do_efsearch == True:
        Lv2_efsearch.efsearch(out_baryfile,n_segments,dper,nphase,nbint,nper,dres,outfile_root,plot_efsearch)

end_time = time.time()

print('Time elapsed: ' + str(end_time-start_time) + ' seconds.')
