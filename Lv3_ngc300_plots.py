#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28th 6:35pm 2019

Script to create these PDF sets:

* soft_1 vs time AND spectrum. Also make sure I do color the corresponding value…
* soft_2 vs time AND spectrum.
* A_band vs time AND spectrum.
* B_band vs time AND spectrum.
* C_band vs time AND spectrum.
* D_band vs time AND spectrum.
* “in_band” vs time AND spectrum.

Want to do then, plots in the shape of the spectrum taking the top half, then the count rate and color for the bottom half…

EDITED on Thurs Aug 8th:

Incorporate Lv2_ngc300_color.py!
"""

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import glob
import Lv0_dirs,Lv2_ngc300_color

Lv0_dirs.global_par()

### this is for the individual GTIs though! Use later...
spectra = sorted(glob.glob(Lv0_dirs.NGC300 + 'bgsub_cl50*.txt'))

################################## FOR COLORS ##################################
bin_size = '05d' #for 1 day
band1 = 'A'
band2 = 'B'
band3 = 'C'
band4 = 'D'

mjds_soft_color,soft_color,soft_color_unc = Lv2_ngc300_color.get_color(bin_size,band1,band2)
mjds_hard_color,hard_color,hard_color_unc = Lv2_ngc300_color.get_color(bin_size,band3,band4)
################################################################################

################################# FOR INTENSITY ################################
intensity_band = 'inband'
counts_file = Lv0_dirs.NGC300 + 'n300_ulx.bgsub_cl50_RGnorm_' + bin_size + '.ffphot'
unc_file = Lv0_dirs.NGC300 + 'n300_ulx.bgsub_cl50_RGerr_' + bin_size + '.ffphot'
mjds_band = np.genfromtxt(counts_file,usecols=(0),unpack=True)
if intensity_band == 'soft1':
    counts_band = np.genfromtxt(counts_file,usecols=(1),unpack=True)
    unc_band = np.genfromtxt(unc_file,usecols=(1),unpack=True)
if intensity_band == 'soft2':
    counts_band = np.genfromtxt(counts_file,usecols=(2),unpack=True)
    unc_band = np.genfromtxt(unc_file,usecols=(2),unpack=True)
if intensity_band == 'A':
    counts_band = np.genfromtxt(counts_file,usecols=(3),unpack=True)
    unc_band = np.genfromtxt(unc_file,usecols=(3),unpack=True)
if intensity_band == 'B':
    counts_band = np.genfromtxt(counts_file,usecols=(4),unpack=True)
    unc_band = np.genfromtxt(unc_file,usecols=(4),unpack=True)
if intensity_band == 'C':
    counts_band = np.genfromtxt(counts_file,usecols=(5),unpack=True)
    unc_band = np.genfromtxt(unc_file,usecols=(5),unpack=True)
if intensity_band == 'D':
    counts_band = np.genfromtxt(counts_file,usecols=(6),unpack=True)
    unc_band = np.genfromtxt(unc_file,usecols=(6),unpack=True)
if intensity_band == 'inband':
    counts_band = np.genfromtxt(counts_file,usecols=(7),unpack=True)
    unc_band = np.genfromtxt(unc_file,usecols=(7),unpack=True)
################################################################################

boundaries = np.array([58230,58280,58290,58300,58310,58335,58340,58355,58380,58420,58460,58520])

fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
ax1.errorbar(mjds_band,counts_band,yerr=unc_band,fmt='^',mfc='none')
ax1.axhline(y=0,ls='--',lw=0.5,alpha=0.5)
ax1.set_ylabel('Counts/s for band: ' + intensity_band, fontsize=12)
for i in range(len(boundaries)):
    ax1.axvline(x=boundaries[i],lw=0.5,alpha=0.5,color='r')
#ax1.set_ylim([-2,4])

ax2.errorbar(mjds_soft_color,soft_color,yerr=soft_color_unc,fmt='^',mfc='none')
ax2.set_xlabel('Time (MJD)',fontsize=12)
ax2.set_ylabel('Soft Color: ' + band2+'/'+band1,fontsize=12)
for i in range(len(boundaries)):
    ax2.axvline(x=boundaries[i],lw=0.5,alpha=0.5,color='r')

ax3.errorbar(mjds_hard_color,hard_color,yerr=hard_color_unc,fmt='^',mfc='none')
ax3.set_xlabel('Time (MJD)',fontsize=12)
ax3.set_ylabel('Hard Color: ' + band4+'/'+band3,fontsize=12)
for i in range(len(boundaries)):
    ax3.axvline(x=boundaries[i],lw=0.5,alpha=0.5,color='r')

plt.subplots_adjust(hspace=0)

################################################################################
## doing a color-color diagram

soft_color_common = []
soft_color_unc_common = []
hard_color_common = []
hard_color_unc_common = []
mjds_common = []
for i in range(len(mjds_hard_color)): #because there are fewer of these...
    if mjds_hard_color[i] in mjds_soft_color:
        mjds_common.append(mjds_hard_color[i])

        soft_color_common.append(soft_color[mjds_soft_color==mjds_hard_color[i]][0])
        soft_color_unc_common.append(soft_color_unc[mjds_soft_color==mjds_hard_color[i]][0])

        hard_color_common.append(hard_color[i])
        hard_color_unc_common.append(hard_color_unc[i])

plt.figure()

mjds_common = np.array(mjds_common)
soft_color_common = np.array(soft_color_common)
soft_color_unc_common = np.array(soft_color_unc_common)
hard_color_common = np.array(hard_color_common)
hard_color_unc_common = np.array(hard_color_unc_common)

soft_counts1 = soft_color_common[(mjds_common>=58230)&(mjds_common<=58280)]
soft_unc_counts1 = soft_color_unc_common[(mjds_common>=58230)&(mjds_common<=58280)]
hard_counts1 = hard_color_common[(mjds_common>=58230)&(mjds_common<=58280)]
hard_unc_counts1 = hard_color_unc_common[(mjds_common>=58230)&(mjds_common<=58280)]

soft_counts2 = soft_color_common[(mjds_common>=58290)&(mjds_common<=58300)]
soft_unc_counts2 = soft_color_unc_common[(mjds_common>=58290)&(mjds_common<=58300)]
hard_counts2 = hard_color_common[(mjds_common>=58290)&(mjds_common<=58300)]
hard_unc_counts2 = hard_color_unc_common[(mjds_common>=58290)&(mjds_common<=58300)]

soft_counts3 = soft_color_common[(mjds_common>=58310)&(mjds_common<=58335)]
soft_unc_counts3 = soft_color_unc_common[(mjds_common>=58310)&(mjds_common<=58335)]
hard_counts3 = hard_color_common[(mjds_common>=58310)&(mjds_common<=58335)]
hard_unc_counts3 = hard_color_unc_common[(mjds_common>=58310)&(mjds_common<=58335)]

soft_counts4 = soft_color_common[(mjds_common>=58340)&(mjds_common<=58355)]
soft_unc_counts4 = soft_color_unc_common[(mjds_common>=58340)&(mjds_common<=58355)]
hard_counts4 = hard_color_common[(mjds_common>=58340)&(mjds_common<=58355)]
hard_unc_counts4 = hard_color_unc_common[(mjds_common>=58340)&(mjds_common<=58355)]

soft_counts5 = soft_color_common[(mjds_common>=58380)&(mjds_common<=58420)]
soft_unc_counts5 = soft_color_unc_common[(mjds_common>=58380)&(mjds_common<=58420)]
hard_counts5 = hard_color_common[(mjds_common>=58380)&(mjds_common<=58420)]
hard_unc_counts5 = hard_color_unc_common[(mjds_common>=58380)&(mjds_common<=58420)]

soft_counts6 = soft_color_common[(mjds_common>=58460)&(mjds_common<=58520)]
soft_unc_counts6 = soft_color_unc_common[(mjds_common>=58460)&(mjds_common<=58520)]
hard_counts6 = hard_color_common[(mjds_common>=58460)&(mjds_common<=58520)]
hard_unc_counts6 = hard_color_unc_common[(mjds_common>=58460)&(mjds_common<=58520)]

#plt.errorbar(hard_color_common,soft_color_common,xerr=hard_color_unc_common,yerr=soft_color_unc_common,fmt='^',mfc='none')
plt.errorbar(hard_counts1,soft_counts1,xerr=hard_unc_counts1,yerr=soft_unc_counts1,fmt='^',mfc='none')
plt.errorbar(hard_counts2,soft_counts2,xerr=hard_unc_counts2,yerr=soft_unc_counts2,fmt='^',mfc='none')
plt.errorbar(hard_counts3,soft_counts3,xerr=hard_unc_counts3,yerr=soft_unc_counts3,fmt='^',mfc='none')
plt.errorbar(hard_counts4,soft_counts4,xerr=hard_unc_counts4,yerr=soft_unc_counts4,fmt='^',mfc='none')
plt.errorbar(hard_counts5,soft_counts5,xerr=hard_unc_counts5,yerr=soft_unc_counts5,fmt='^',mfc='none')
plt.errorbar(hard_counts6,soft_counts6,xerr=hard_unc_counts6,yerr=soft_unc_counts6,fmt='^',mfc='none')
plt.legend(('58230-58280','58290-58300','58310-58335','58340-58355','58380-58420','58460-58520'),loc='best')

plt.xlabel('Hard Color: ' + band4+'/'+band3,fontsize=12)
plt.ylabel('Soft Color: ' + band2+'/'+band1,fontsize=12)

################################################################################
## doing a color-intensity diagrams. Can check with Ron about this, but to maximize
## the number of points, don't just plot colors whereby you have BOTH soft AND hard colors.
## The main thing is that for the soft color, both A AND B bands had to have >0 count rate.
## Construct them separately

#### first intensity vs soft color
soft_color_common = []
soft_color_unc_common = []
intensity_common = []
intensity_unc_common = []
mjds_common = []
for i in range(len(mjds_soft_color)): #because there are fewer of these...
    if mjds_soft_color[i] in mjds_band:
        mjds_common.append(mjds_soft_color[i])

        soft_color_common.append(soft_color[i])
        soft_color_unc_common.append(soft_color_unc[i])

        intensity_common.append(counts_band[mjds_band==mjds_soft_color[i]][0])
        intensity_unc_common.append(unc_band[mjds_band==mjds_soft_color[i]][0])

plt.figure()

mjds_common = np.array(mjds_common)
soft_color_common = np.array(soft_color_common)
soft_color_unc_common = np.array(soft_color_unc_common)
intensity_common = np.array(intensity_common)
intensity_unc_common = np.array(intensity_unc_common)

soft_counts1 = soft_color_common[(mjds_common>=58230)&(mjds_common<=58280)]
soft_unc_counts1 = soft_color_unc_common[(mjds_common>=58230)&(mjds_common<=58280)]
intensity1 = intensity_common[(mjds_common>=58230)&(mjds_common<=58280)]
intensity_unc1 = intensity_unc_common[(mjds_common>=58230)&(mjds_common<=58280)]

soft_counts2 = soft_color_common[(mjds_common>=58290)&(mjds_common<=58300)]
soft_unc_counts2 = soft_color_unc_common[(mjds_common>=58290)&(mjds_common<=58300)]
intensity2 = intensity_common[(mjds_common>=58290)&(mjds_common<=58300)]
intensity_unc2 = intensity_unc_common[(mjds_common>=58290)&(mjds_common<=58300)]

soft_counts3 = soft_color_common[(mjds_common>=58310)&(mjds_common<=58335)]
soft_unc_counts3 = soft_color_unc_common[(mjds_common>=58310)&(mjds_common<=58335)]
intensity3 = intensity_common[(mjds_common>=58310)&(mjds_common<=58335)]
intensity_unc3 = intensity_unc_common[(mjds_common>=58310)&(mjds_common<=58335)]

soft_counts4 = soft_color_common[(mjds_common>=58340)&(mjds_common<=58355)]
soft_unc_counts4 = soft_color_unc_common[(mjds_common>=58340)&(mjds_common<=58355)]
intensity4 = intensity_common[(mjds_common>=58340)&(mjds_common<=58355)]
intensity_unc4 = intensity_unc_common[(mjds_common>=58340)&(mjds_common<=58355)]

soft_counts5 = soft_color_common[(mjds_common>=58380)&(mjds_common<=58420)]
soft_unc_counts5 = soft_color_unc_common[(mjds_common>=58380)&(mjds_common<=58420)]
intensity5 = intensity_common[(mjds_common>=58380)&(mjds_common<=58420)]
intensity_unc5 = intensity_unc_common[(mjds_common>=58380)&(mjds_common<=58420)]

soft_counts6 = soft_color_common[(mjds_common>=58460)&(mjds_common<=58520)]
soft_unc_counts6 = soft_color_unc_common[(mjds_common>=58460)&(mjds_common<=58520)]
intensity6 = intensity_common[(mjds_common>=58460)&(mjds_common<=58520)]
intensity_unc6 = intensity_unc_common[(mjds_common>=58460)&(mjds_common<=58520)]

#plt.errorbar(hard_color_common,soft_color_common,xerr=hard_color_unc_common,yerr=soft_color_unc_common,fmt='^',mfc='none')
plt.errorbar(soft_counts1,intensity1,xerr=soft_unc_counts1,yerr=intensity_unc1,fmt='^',mfc='none')
plt.errorbar(soft_counts2,intensity2,xerr=soft_unc_counts2,yerr=intensity_unc2,fmt='^',mfc='none')
plt.errorbar(soft_counts3,intensity3,xerr=soft_unc_counts3,yerr=intensity_unc3,fmt='^',mfc='none')
plt.errorbar(soft_counts4,intensity4,xerr=soft_unc_counts4,yerr=intensity_unc4,fmt='^',mfc='none')
plt.errorbar(soft_counts5,intensity5,xerr=soft_unc_counts5,yerr=intensity_unc5,fmt='^',mfc='none')
plt.errorbar(soft_counts6,intensity6,xerr=soft_unc_counts6,yerr=intensity_unc6,fmt='^',mfc='none')
plt.legend(('58230-58280','58290-58300','58310-58335','58340-58355','58380-58420','58460-58520'),loc='best')

plt.xlabel('Soft Color: ' + band2+'/'+band1,fontsize=12)
plt.ylabel('Intensity (ct/s)',fontsize=12)
#plt.xlim([-0.1,2.9])
#plt.ylim([-0.9,1.6])

################################################################################

#### now intensity vs hard color
hard_color_common = []
hard_color_unc_common = []
intensity_common = []
intensity_unc_common = []
mjds_common = []
for i in range(len(mjds_hard_color)): #because there are fewer of these...
    if mjds_hard_color[i] in mjds_band:
        mjds_common.append(mjds_hard_color[i])

        hard_color_common.append(hard_color[i])
        hard_color_unc_common.append(hard_color_unc[i])

        intensity_common.append(counts_band[mjds_band==mjds_hard_color[i]][0])
        intensity_unc_common.append(unc_band[mjds_band==mjds_hard_color[i]][0])

plt.figure()

mjds_common = np.array(mjds_common)
hard_color_common = np.array(hard_color_common)
hard_color_unc_common = np.array(hard_color_unc_common)
intensity_common = np.array(intensity_common)
intensity_unc_common = np.array(intensity_unc_common)

hard_counts1 = hard_color_common[(mjds_common>=58230)&(mjds_common<=58280)]
hard_unc_counts1 = hard_color_unc_common[(mjds_common>=58230)&(mjds_common<=58280)]
intensity1 = intensity_common[(mjds_common>=58230)&(mjds_common<=58280)]
intensity_unc1 = intensity_unc_common[(mjds_common>=58230)&(mjds_common<=58280)]

hard_counts2 = hard_color_common[(mjds_common>=58290)&(mjds_common<=58300)]
hard_unc_counts2 = hard_color_unc_common[(mjds_common>=58290)&(mjds_common<=58300)]
intensity2 = intensity_common[(mjds_common>=58290)&(mjds_common<=58300)]
intensity_unc2 = intensity_unc_common[(mjds_common>=58290)&(mjds_common<=58300)]

hard_counts3 = hard_color_common[(mjds_common>=58310)&(mjds_common<=58335)]
hard_unc_counts3 = hard_color_unc_common[(mjds_common>=58310)&(mjds_common<=58335)]
intensity3 = intensity_common[(mjds_common>=58310)&(mjds_common<=58335)]
intensity_unc3 = intensity_unc_common[(mjds_common>=58310)&(mjds_common<=58335)]

hard_counts4 = hard_color_common[(mjds_common>=58340)&(mjds_common<=58355)]
hard_unc_counts4 = hard_color_unc_common[(mjds_common>=58340)&(mjds_common<=58355)]
intensity4 = intensity_common[(mjds_common>=58340)&(mjds_common<=58355)]
intensity_unc4 = intensity_unc_common[(mjds_common>=58340)&(mjds_common<=58355)]

hard_counts5 = hard_color_common[(mjds_common>=58380)&(mjds_common<=58420)]
hard_unc_counts5 = hard_color_unc_common[(mjds_common>=58380)&(mjds_common<=58420)]
intensity5 = intensity_common[(mjds_common>=58380)&(mjds_common<=58420)]
intensity_unc5 = intensity_unc_common[(mjds_common>=58380)&(mjds_common<=58420)]

hard_counts6 = hard_color_common[(mjds_common>=58460)&(mjds_common<=58520)]
hard_unc_counts6 = hard_color_unc_common[(mjds_common>=58460)&(mjds_common<=58520)]
intensity6 = intensity_common[(mjds_common>=58460)&(mjds_common<=58520)]
intensity_unc6 = intensity_unc_common[(mjds_common>=58460)&(mjds_common<=58520)]

#plt.errorbar(hard_color_common,hard_color_common,xerr=hard_color_unc_common,yerr=hard_color_unc_common,fmt='^',mfc='none')
plt.errorbar(hard_counts1,intensity1,xerr=hard_unc_counts1,yerr=intensity_unc1,fmt='^',mfc='none')
plt.errorbar(hard_counts2,intensity2,xerr=hard_unc_counts2,yerr=intensity_unc2,fmt='^',mfc='none')
plt.errorbar(hard_counts3,intensity3,xerr=hard_unc_counts3,yerr=intensity_unc3,fmt='^',mfc='none')
plt.errorbar(hard_counts4,intensity4,xerr=hard_unc_counts4,yerr=intensity_unc4,fmt='^',mfc='none')
plt.errorbar(hard_counts5,intensity5,xerr=hard_unc_counts5,yerr=intensity_unc5,fmt='^',mfc='none')
plt.errorbar(hard_counts6,intensity6,xerr=hard_unc_counts6,yerr=intensity_unc6,fmt='^',mfc='none')
plt.legend(('58230-58280','58290-58300','58310-58335','58340-58355','58380-58420','58460-58520'),loc='best')

plt.xlabel('Hard Color: ' + band4+'/'+band3,fontsize=12)
plt.ylabel('Intensity (ct/s)',fontsize=12)

plt.show()

"""

RGcms = '/Volumes/Samsung_T5/NGC300_ULX/n300_ulx.bgsub_cl50_RGcms.ffphot'
RGnorm = '/Volumes/Samsung_T5/NGC300_ULX/n300_ulx.bgsub_cl50_RGnorm.ffphot'
RGerror = '/Volumes/Samsung_T5/NGC300_ULX/n300_ulx.bgsub_cl50_RGerr_norm.ffphot'

soft1,soft2,A_band,B_band,C_band,D_band,inband,time = np.genfromtxt(RGcms,usecols=(3,4,5,6,7,8,9,11),unpack=True)
soft1_err,soft2_err,A_err,B_err,C_err,D_err,inband_err,time_error = np.genfromtxt(RGerror,usecols=(3,4,5,6,7,8,9,11),unpack=True)
# time_error does not mean error in the time value, I just mean the time value found in the error text file


def create_counts_pdf(filename,time,counts,counts_error,spectra):

    Creates a 2x1 subplot where the top plot is the spectrum, and the bottom plot
    is the light curve from array_counts

    filename - desired filename for the PDF
    time - array of time values
    array_counts - array of counts from the desired band
    spectra - list of spectra text files

    pdf_filename = '/Volumes/Samsung_T5/NGC300_ULX/plots/' + filename + '.pdf'

    with PdfPages(pdf_filename) as pdf:
        for i in tqdm(range(len(spectra))):
            if i != 397:
                E,E_err,rate,rate_error = np.genfromtxt(spectra[i],usecols=(0,1,2,3),unpack=True)

                f,(ax1,ax2) = plt.subplots(2,1)
                ax1.errorbar(x=E,y=rate,xerr=E_err,yerr=rate_error,fmt='+')
                ax1.set_xlim([0.3,12])
                ax1.set_xscale('log')
                ax2.errorbar(x=time,y=counts,fmt='x')
                #ax2.errorbar(x=time[i],y=counts[i],fmt='x')

                pdf.savefig()
                plt.close()

    return

create_counts_pdf('soft1',time,soft1,soft1_err,spectra)

"""
#do PdfPages
#then do the counts on top and spectrum below! Do plt.subplot or something
# make sure for the counts, I plot ALL and them, AND color the one where the spectrum is being shown
