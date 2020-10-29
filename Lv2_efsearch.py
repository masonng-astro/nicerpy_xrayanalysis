#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Saturday Jan 11 7:22pm 2020

Running epoch folding searches through FTOOLS' efsearch!

The algorithm works in the following way (thinking of RXJ0209 as an example):
- It will look at the first segment (defined by nbint) of the data and do the
epoch folding there. It will then run through the time series until there are photons,
and defines a new segment starting from there, with a segment length defined before.
It repeats until it ends. So unless you have 100% duty cycle, you won't always get
the number of segments you're after! This also does mean that not all segments
are of equal length, which means no. of events vary.

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import pathlib
import Lv0_dirs,Lv0_fits2dict
import Lv1_barycorr,Lv1_data_bin
import Lv2_lc,Lv2_ps,Lv2_color,Lv2_phase
import Lv3_E_boundary
from matplotlib.backends.backend_pdf import PdfPages
import time
from tqdm import tqdm
import peakutils
import matplotlib.pyplot as plt
import os
import subprocess
import glob

def efsearch(eventfile,n_segments,dper,nphase,nbint,nper,dres,outfile_root,plot_efsearch):
    """
    Performing FTOOLS' efsearch!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    n_segments - no. of segments to break the epoch folding search into
    dper - value for period used in the folding; input represents center of range of trial periods
    nphase - no. of phases in the folded light curve(s) ; "INDEF" uses the default value
    nbint - the number of newbins per interval used in the analysis
    nper - no. of periods over which the search is carried out
    dres - period resolution is the spacing between two contiguous periods in the search
           "INDEF" uses default value of half the Fourier resolution in the interval
    outfile_root - prefix for the end file name
    plot_efsearch - 'yes' or 'no' to plot the results from efsearch; do "exit" for the next plot!
    """

    efsearch_cmd = ['efsearch',eventfile,'window="-"','sepoch=INDEF','dper='+str(dper),'nphase='+str(nphase),'nbint='+str(nbint),'nper='+str(nper),'dres='+str(dres),'outfile='+outfile_root,'outfiletype=2','plot='+plot_efsearch]

    parent_folder = str(pathlib.Path(eventfile).parent)
    efsearch_logfile = parent_folder + '/efsearch.log'

    with open(efsearch_logfile,'w') as logtextfile:
        output = subprocess.run(efsearch_cmd,capture_output=True,text=True)
        logtextfile.write(output.stdout)
        logtextfile.write('*------------------------------* \n')
        logtextfile.write(output.stderr)
    logtextfile.close()

    subprocess.run(['mv',outfile_root+'.fes',parent_folder + '/' + outfile_root + '.fes'])

    efsearch_results = fits.open(parent_folder + '/' + outfile_root + '.fes')
    efsearch_plot = parent_folder + '/' + outfile_root + '.pdf'

    with PdfPages(efsearch_plot) as pdf:
        for i in tqdm(range(1,len(efsearch_results))):
            T_start = efsearch_results[i].header['TSTARTI'] + efsearch_results[i].header['TSTARTF']
            T_stop = efsearch_results[i].header['TSTOPI'] + efsearch_results[i].header['TSTOPF']
            T_zeroref = efsearch_results[i].header['TIMEZERI'] + efsearch_results[i].header['TIMEZERF']

            ### extracting period, chi-squared (+ error) values from the efsearch results
            period = efsearch_results[i].data['PERIOD']
            chisqrd1 = efsearch_results[i].data['CHISQRD1']
            error1 = efsearch_results[i].data['ERROR1']

            ### First find the peaks in the plot
            peak_indices = peakutils.indexes(chisqrd1,thres=0.2,min_dist=0.0001)

            ### plotting the chi-squared vs period plots and the corresponding Gaussian fits
            plt.errorbar(x=period,y=chisqrd1,yerr=error1,fmt='x')
            plt.title('TSTART: ' + str((T_start-T_zeroref)*86400) + ' ; TSTOP: ' + str((T_stop-T_zeroref)*86400))
            plt.xlabel('Period (s)',fontsize=12)
            plt.ylabel('chi-squared',fontsize=12)
            pdf.savefig()
            plt.close()

if __name__ == "__main__":
    eventfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/ngc300x1/ngc300x1_merge_niceroverlap_all.evt'
    n_segments = 1
    dper = 1/8.461949e-6
    nphase = 20
    nper = 128
    dres = 1
    plot_efsearch = "yes"

    T = fits.open(eventfile)[2].data['STOP'][-1] - fits.open(eventfile)[2].data['START'][0]

    outfile_root = 'ngc300x1' + '_' + str(n_segments) + 'segs_' + str(nper)
    nbint = int((T/(dper/nphase))/n_segments)

    efsearch(eventfile,n_segments,dper,nphase,nbint,nper,dres,outfile_root,plot_efsearch)
