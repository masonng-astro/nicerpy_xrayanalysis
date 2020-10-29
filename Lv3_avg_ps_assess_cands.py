#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 4:08pm 2020

Script that creates a heat map (histogram with 3rd dimension) showing on the y axis,
log10(f), and has either W (fixed threshold) or threshold (fixed W) on the x axis.
The color bar will correspond to significances.

"""
from __future__ import division, print_function
import numpy as np
import subprocess
import pathlib
from tqdm import tqdm

import Lv0_dirs,Lv2_TBOs_method
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import mplcursors

Lv0_dirs.global_par()

def assess_cands(eventfile,segment_length,PI1,PI2,hist_min_f,W,threshold):
    """
    Creating the heat map

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    segment_length - desired length of segments in seconds
    PI1 - lower bound of PI (not energy in keV!) desired for the energy range
    PI2 - upper bound of PI (not energy in keV!) desired for the energy range
    hist_min_f - minimum frequency to calculate the histogram from
    W - list of W values (number of consecutive Fourier bins to average over)
    threshold - list of thresholds for counts in each segment (in 1s bins; in terms of a percentage)
    """
    parent_folder = str(pathlib.Path(eventfile).parent)

    powers = []
    for k in range(len(W)):
        for l in range(len(threshold)):
            if PI1 != '':
                f,ps = np.genfromtxt(parent_folder+'/S'+str(segment_length)+'_W'+str(W[k])+'_T'+str(threshold[l])+'_E'+str(PI1)+'-'+str(PI2)+'.txt',usecols=(0,1),unpack=True)
            else:
                f,ps = np.genfromtxt(parent_folder+'/S'+str(segment_length)+'_W'+str(W[k])+'_T'+str(threshold[l]) + '.txt',usecols=(0,1),unpack=True)
            cand_f = f[f>hist_min_f]
            cand_ps = ps[f>hist_min_f]
            powers += list(cand_ps)
    powers = np.array(powers)
    powers_min = np.min(powers[powers!=np.inf])
    powers_max = np.max(powers[powers!=np.inf])

    ### keeping W constant
    for k in range(len(W)):
        if PI1 != '':
            filename = parent_folder + '/S'+str(segment_length)+'_W'+str(W[k])+'_E'+str(PI1)+'-'+str(PI2)+'.pdf'
        else:
            filename = parent_folder + '/S'+str(segment_length)+'_W'+str(W[k])+'.pdf'
        with PdfPages(filename) as pdf:
            for l in range(len(threshold)):
                if PI1 != '':
                    f,ps = np.genfromtxt(parent_folder+'/S'+str(segment_length)+'_W'+str(W[k])+'_T'+str(threshold[l])+'_E'+str(PI1)+'-'+str(PI2)+'.txt',usecols=(0,1),unpack=True)
                else:
                    f,ps = np.genfromtxt(parent_folder+'/S'+str(segment_length)+'_W'+str(W[k])+'_T'+str(threshold[l]) + '.txt',usecols=(0,1),unpack=True)
                cand_f = f[f>hist_min_f]
                cand_ps = ps[f>hist_min_f]

                plt.scatter(x=np.ones(len(cand_f))*threshold[l],y=cand_f,c=cand_ps,marker='o',cmap=cm.gist_heat,vmin=powers_min,vmax=powers_max,edgecolor='k')
                if PI1 != '':
                    title = 'Segment Length: ' + str(segment_length) + 's; E = ' + str(PI1) + '-' + str(PI2) + '; W = '+str(W[k])
                else:
                    title = 'Segment Length: ' + str(segment_length) + 's; W = '+str(W[k])
                plt.title(title,fontsize=12)
                plt.xlabel('Threshold (percentage)',fontsize=12)
                plt.ylabel('log10(f)',fontsize=12)
                plt.yscale('log')
            plt.colorbar()
            pdf.savefig()
            plt.close()

    ### keeping threshold constant
    for l in range(len(threshold)):
        if PI1 != '':
            filename = parent_folder + '/S'+str(segment_length)+'_T'+str(threshold[l])+'_E'+str(PI1)+'-'+str(PI2)+'.pdf'
        else:
            filename = parent_folder + '/S'+str(segment_length)+'_T'+str(threshold[l])+'.pdf'
        with PdfPages(filename) as pdf:
            for k in range(len(W)):
                if PI1 != '':
                    f,ps = np.genfromtxt(parent_folder+'/S'+str(segment_length)+'_W'+str(W[k])+'_T'+str(threshold[l])+'_E'+str(PI1)+'-'+str(PI2)+'.txt',usecols=(0,1),unpack=True)
                else:
                    f,ps = np.genfromtxt(parent_folder+'/S'+str(segment_length)+'_W'+str(W[k])+'_T'+str(threshold[l]) + '.txt',usecols=(0,1),unpack=True)
                cand_f = f[f>hist_min_f]
                cand_ps = ps[f>hist_min_f]

                plt.scatter(x=np.ones(len(cand_f))*W[k],y=cand_f,c=cand_ps,marker='o',cmap=cm.gist_heat,vmin=powers_min,vmax=powers_max,edgecolor='k')
                if PI1 != '':
                    title = 'Segment Length: ' + str(segment_length) + 's; E = ' + str(PI1) + '-' + str(PI2) + '; Threshold = '+str(threshold[l])
                else:
                    title = 'Segment Length: ' + str(segment_length) + 's; Threshold = '+str(threshold[l])
                plt.title(title,fontsize=12)
                plt.xlabel('W',fontsize=12)
                plt.ylabel('log10(f)',fontsize=12)
                plt.yscale('log')
            plt.colorbar()
            pdf.savefig()
            plt.close()

    ### keeping the threshold constant

if __name__ == "__main__":
    merged_id = '000019'
    eventfile = Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/merged' + merged_id + '/merged' + merged_id + '_nicersoft_bary.evt'
    #eventfile = '/Volumes/Samsung_T5/NICER-data/xtej1739_mostrec/2002131540_filt_bary.evt'

    segment_lengths = [500] #desired length of segments in seconds
    PI1 = [30,200]
    PI2 = [200,1200]
    threshold = [5,10,20,30,40,50] #threshold for counts in each segment (in 1s bins; in terms of a percentage)
    W = [1,5,10] #number of consecutive Fourier bins to average over
    hist_min_f = 1

    for i in tqdm(range(len(segment_lengths))):
        for j in tqdm(range(len(PI1))):
            assess_cands(eventfile,segment_lengths[i],PI1[j],PI2[j],hist_min_f,W,threshold)
