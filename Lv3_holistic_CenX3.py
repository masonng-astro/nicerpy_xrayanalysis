#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 11:26am 2019

Creating a holistic analysis for Cen X-3, which includes light curves, color-intensity diagrams,
pulse profiles, and spectra. The automation will go over the first 3 data products.

It will be similar in form to Lv3_main.py, but more tailored and more loopy (i.e., for loops).

"""
from __future__ import division, print_function
import numpy as np
import time

import Lv0_dirs,Lv0_obsid_info
import Lv0_call_eventcl,Lv0_call_att,Lv0_call_hk,Lv0_call_mkf,Lv0_call_orb,Lv0_call_uf,Lv0_call_ufa
import Lv1_data_bin,Lv1_data_gtis,Lv1_data_spectra
import Lv2_lc,Lv2_ps,Lv2_color,Lv2_phase,Lv2_sources
import Lv3_E_boundary,Lv3_diagnostics,Lv3_diagnostics_display

from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt

Lv0_dirs.global_par()

obsids = ['0034070101','0034070102','0034070103','0034070104','1034070101',
          '1034070102','1034070103','1034070104','1034070105','1034070106']
obsids_freq = ['0034070101','0034070102','0034070103','0034070104','1034070101',
               '1034070106']

#pulse_f =
#pulse_f_error =
tbin_size = 1 #1-second bins for the light curves
par_list = ['TIME','PI','PI_FAST']

bary = True

## Getting the boundary energy for the color-intensity diagrams
E1 = 0.3
E2 = 12
cut_type = 'manual'
bound = 2.7

Ebin_size = 0.05
shift = 0.2 #for pulse profile
no_phase_bins = 51

for i in range(len(obsids)):
    dir = Lv0_dirs.BASE_DIR+'outputs/holistic_analysis_CenX3/'
    filename = dir + obsids[i] + '_bary_bin' + str(tbin_size) + 's.pdf'
    start_end = Lv0_obsid_info.obstime(obsids[i])
    radec = Lv0_obsid_info.ra_dec(obsids[i])
    peakf = Lv0_obsid_info.peak_freq(obsids[i])
    gtis_shifted,gtis_unshifted = Lv1_data_gtis.get_gtis(obsids[i],bary,50)
    Ebound = Lv3_E_boundary.E_bound(obsids[i],bary,par_list,E1,E2,cut_type,bound)
    obj_name = Lv2_sources.obsid_to_obj(obsids[i])

    with PdfPages(filename) as pdf:
        ### COVER PAGE
        cover = plt.figure(figsize=(8.27,11.69))
        cover.clf()
        txt = obsids[i]
        cover.text(0.5,0.95,txt,transform=cover.transFigure,size=24,ha='center')
        cover.text(0.05,0.9,'Start time (UTC): ' + start_end[0] + '; End time (UTC): ' + start_end[1])
        cover.text(0.05,0.875,'Target RA: ' + radec[0] + '; Dec: ' + radec[1])
        cover.text(0.05,0.85,'Peak frequency: ' + peakf[0] + ' +- ' + peakf[1])
        cover.text(0.05,0.80,'Shifted GTIs > 50s: ' + str(gtis_shifted))
        cover.text(0.05,0.775,'Unshifted GTIs > 50s: \n')
        cover.text(0.05,0.75,str(gtis_unshifted))

        pdf.savefig()
        plt.close()

        ### SECOND PAGE
        plt.figure(figsize=(10,8))
        Lv2_lc.whole(obsids[i],bary,par_list,tbin_size,'save')
        pdf.savefig()
        plt.close()

        ### THIRD PAGE - PULSE PROFILE
        if obsids[i] in obsids_freq:
            plt.figure(figsize=(10,8))
            Lv2_phase.partial_E(obsids[i],bary,par_list,tbin_size,Ebin_size,np.float(peakf[0]),shift,no_phase_bins,0.3,Ebound,'overlap')
            Lv2_phase.partial_E(obsids[i],bary,par_list,tbin_size,Ebin_size,np.float(peakf[0]),shift,no_phase_bins,Ebound,12,'overlap')
            Lv2_phase.partial_E(obsids[i],bary,par_list,tbin_size,Ebin_size,np.float(peakf[0]),shift,no_phase_bins,0.3,12,'overlap')
            plt.title('Pulse profile for ' + obj_name + ', ObsID ' + str(obsids[i])+ '\n Energy range: ' + str(E1) + 'keV - ' + str(E2) + 'keV',fontsize=12)
            plt.legend((str(E1) + '-'+str(Ebound) + ' keV',str(Ebound)+'-'+str(E2)+' keV', str(E1)+'-'+str(E2)+' keV'),loc='best')
            pdf.savefig()
            plt.close()

        for j in range(0,len(gtis_shifted),2):
            plt.figure(figsize=(8.27,11.69))
            Lv2_lc.partial_t(obsids[i],bary,par_list,tbin_size,gtis_shifted[j],gtis_shifted[j+1],'save')
            pdf.savefig()
            plt.close()

            plt.figure(figsize=(8.27,11.69))
            Lv2_color.plotting_t(obsids[i],bary,par_list,Ebound,tbin_size,gtis_shifted[j],gtis_shifted[j+1],'save')
            pdf.savefig()
            plt.close()
