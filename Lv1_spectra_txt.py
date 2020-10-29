#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28th 4:57pm 2019

Temporary script to obtain spectra, from background-subtracted spectra, that have
the response and arf matrices folded in.

The script will generate the .xcm file needed to be read into XSPEC!
Will also break up the big all_spectra.txt file into individual text files!

"""
from __future__ import division, print_function
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import glob
from mpldatacursor import datacursor
import Lv0_dirs
import pathlib
from matplotlib.backends.backend_pdf import PdfPages

Lv0_dirs.global_par()

"""
response = '/Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf'
arf = '/Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf'

spectra = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/bgsub_cl50*.pha'))

data_list = []
response_list = []
arf_list = []

for i in range(len(spectra)):
    data_list.append(str(i+1)+':'+str(i+1))
    data_list.append(spectra[i])

    response_list.append(str(i+1))
    response_list.append(response)

    arf_list.append(str(i+1))
    arf_list.append(arf)

data_string = 'data ' + ' '.join([data_list[i] for i in range(len(data_list))])
response_string = 'response ' + ' '.join(response_list[i] for i in range(len(response_list)))
arf_string = 'arf ' + ' '.join(arf_list[i] for i in range(len(arf_list)))
ignore_string = 'ignore **:0.0-0.285,12.01-**'

xcm_file = open('get_spectra.xcm','w')
xcm_file.write(data_string + '\n')
xcm_file.write(response_string + '\n')
xcm_file.write(arf_string + '\n')
xcm_file.write(ignore_string + '\n')
xcm_file.close()
"""

##### Get spectra, and save them into individual files. Remember that I got them
##### simultaneously in XSPEC!

"""
spectra_bgsub = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/bgsub_cl50*.pha'))
spectra_all = '/Volumes/Samsung_T5/NGC300_ULX/all_spectra.txt'
spectra = np.array(open(spectra_all,'r').read().split('\n'))

separator = list(np.where(spectra=='NO NO NO NO')[0])
separator.insert(0,2)
separator.insert(len(separator),len(spectra)-1)

for i in tqdm(range(len(spectra_bgsub))):
    if i != 397:
        spectrum_file = spectra_bgsub[i][:-3] + 'txt'
        spectrum = open(spectrum_file,'w')

        if i < 397:
            for j in range(separator[i]+1,separator[i+1]):
                spectrum.write(spectra[j] + '\n')
        else:
            for j in range(separator[i-1]+1,separator[i]):
                spectrum.write(spectra[j] + '\n')

        spectrum.close()
"""

"""
##### Do the same thing, but for the 5-day binned data
binspectra_bgsub = sorted(glob.glob('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_*_3C50x_*.pha'))
spectra_all = '/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_3C50x_spectra.txt'
spectra = np.array(open(spectra_all,'r').read().split('\n'))

separator = list(np.where(spectra=='NO NO NO NO NO NO')[0])
separator.insert(0,2)
separator.insert(len(separator),len(spectra)-1)

for i in tqdm(range(len(binspectra_bgsub))):
    spectrum_file = binspectra_bgsub[i][:-3] + 'txt'
    spectrum = open(spectrum_file,'w')

    for j in range(separator[i]+1,separator[i+1]):
        spectrum.write(spectra[j] + '\n')

    spectrum.close()
"""

"""
#not_to_use = ['58279','58334','58399','58464','58514','58599','58604']
not_to_use = []
#total_E = []
##### time to plot all 41 spectra and save them onto PDFs...!
#binspectra_cl = sorted(glob.glob('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_58*_cl_*.txt'))
#binspectra_cl = sorted(glob.glob('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/grp_58*_bgsub_*.txt'))
binspectra_cl = sorted(glob.glob('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_58*_3C50x_*.txt'))


with PdfPages('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_3C50x_binned_spectra.pdf') as pdf:
    for i in tqdm(range(len(binspectra_cl))):
        if binspectra_cl[i][-24:-19] not in not_to_use:
            E_cl,E_cl_unc,flux_cl,flux_cl_unc,bg,bg_unc = np.genfromtxt(binspectra_cl[i],usecols=(0,1,2,3,4,5),unpack=True,skip_footer=1)
            plt.errorbar(E_cl,flux_cl,xerr=E_cl_unc,yerr=flux_cl_unc,color='r')#,drawstyle='steps-mid')
            plt.errorbar(E_cl,bg,xerr=E_cl_unc,yerr=bg_unc,color='b')#,drawstyle='steps-mid'
            plt.xscale('log')
            plt.title('File name: ' + str(pathlib.Path(binspectra_cl[i]).name),fontsize=12)
            plt.xlabel('Energy (keV)',fontsize=10)
            plt.xlim([0.3,12])
            plt.ylabel('Normalized counts /s/keV',fontsize=10)
            plt.axhline(y=0,lw=0.5,alpha=0.5)
            #plt.legend(('cl','bg'),loc='best')
            pdf.savefig()
            plt.close()
"""


#determined on 10/3
filtered_mjds = ['58239','58244','58249','58254','58259','58264','58269','58274',
        '58289','58309','58314','58324','58329','58334','58339','58389',
        '58399','58449','58454','58459','58464','58484','58489','58504',
        '58509']

#binspectra_cl = [binspectra_cl[i] for i in range(len(binspectra_cl)) if binspectra_cl[i][-21:-16] in filtered_mjds]

"""
### to show
for i in tqdm(range(len(binspectra_cl))):
    #plt.figure(i)
    if '58389' in binspectra_cl[i]:
        plt.axhline(y=0,alpha=0.5,lw=0.5)
        E_cl,E_cl_unc,flux_cl,flux_cl_unc = np.genfromtxt(binspectra_cl[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)
        #E_bg,E_bg_unc,flux_bg,flux_bg_unc = np.genfromtxt(binspectra_bg[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)
        E_bgsub,E_bgsub_unc,flux_bgsub,flux_bgsub_unc = np.genfromtxt(binspectra_bgsub[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)
        plt.errorbar(E_cl,flux_cl,xerr=E_cl_unc,yerr=flux_cl_unc,fmt='r+-')
        #plt.errorbar(E_bg,flux_bg,xerr=E_bg_unc,yerr=flux_bg_unc,fmt='b+-')
        #plt.errorbar(E_bgsub,flux_bgsub,xerr=E_bgsub_unc,yerr=flux_bgsub_unc,fmt='b+-')
        plt.title('File name: ' + binspectra_cl[i] + '\n' + 'Background: ' + binspectra_bg[i],fontsize=12)
        plt.xlabel('Energy (keV)',fontsize=10)
        plt.ylabel('Normalized counts /s/keV',fontsize=10)
        plt.legend(('y=0','cl','bgsub'),loc='best')
        plt.xlim([0.3,12])
plt.show()
"""

"""
##### 9/29/2020 - compare 3C50 and 3C50 plot_efsearch

binspectra_3C50 = sorted(glob.glob('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_58*_3C50_*.txt'))
binspectra_3C50x = sorted(glob.glob('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_58*_3C50x_*.txt'))

with PdfPages('/Volumes/Samsung_T5/n300_ulx_2020/spectra_05d/jsgrp_3C50_3C50x_altplot.pdf') as pdf:
    for i in tqdm(range(len(binspectra_3C50))):
        E1,E1_unc,flux1,flux1_unc,bg1,bg1_unc = np.genfromtxt(binspectra_3C50[i],usecols=(0,1,2,3,4,5),unpack=True,skip_footer=1)
        E2,E2_unc,flux2,flux2_unc,bg2,bg2_unc = np.genfromtxt(binspectra_3C50x[i],usecols=(0,1,2,3,4,5),unpack=True,skip_footer=1)
        plt.errorbar(E1,flux1,xerr=E1_unc,yerr=flux1_unc,color='r',drawstyle='steps-mid')
        plt.errorbar(E2,flux2,xerr=E2_unc,yerr=flux2_unc,color='b',drawstyle='steps-mid')
        plt.xscale('log')
        plt.title('File name: ' + str(pathlib.Path(binspectra_3C50[i]).name),fontsize=12)
        plt.xlabel('Energy (keV)',fontsize=10)
        plt.xlim([0.3,12])
        plt.ylabel('Normalized counts /s/keV',fontsize=10)
        plt.legend(('3C50','3C50 + X-1'),loc='best')
        plt.axhline(y=0,lw=0.5,alpha=0.5)
        pdf.savefig()
        plt.close()
"""

"""
##### 9/29/2020 - generate photometry files like Ron's!
xbgsub_fffphot = open(Lv0_dirs.NGC300_2020 + 'n300_ulx.xsbgsub_cl50_g2020norm.fffphot','w')
bgsub_file = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.bgsub_cl50_g2020norm.fffphot',dtype='str',usecols=(0),unpack=True)
exp,telapse,mjd,norm = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.bgsub_cl50_g2020norm.fffphot',dtype='float64',usecols=(1,2,11,12),unpack=True)
for i in tqdm(range(len(bgsub_file))):
    xbgsub = Lv0_dirs.NGC300_2020 + 'bgsub_cl50/xspha/xs' + bgsub_file[i]
    rates = fits.open(xbgsub)[1].data['RATE']
    chans = fits.open(xbgsub)[1].data['CHANNEL']
    s1 = str(round(sum(rates[20-1:30-1])*norm[i],4))
    s2 = str(round(sum(rates[30-1:40-1])*norm[i],4))
    a = str(round(sum(rates[40-1:100-1])*norm[i],4))
    b = str(round(sum(rates[100-1:200-1])*norm[i],4))
    c = str(round(sum(rates[200-1:400-1])*norm[i],4))
    d = str(round(sum(rates[400-1:1200-1])*norm[i],4))
    inband = str(round(sum(rates[40-1:1200-1])*norm[i],4))
    hbg = str(round(sum(rates[1300-1:1501-1])*norm[i],4))
    xbgsub_fffphot.write(str(pathlib.Path(xbgsub).name) + '  ' + str(exp[i]) + '  ' + str(telapse[i]) + '     ' + str(s1) + ' ' + str(s2) + '      ' + str(a) + ' ' + str(b) + ' ' + str(c) + ' ' + str(d) + '      ' + str(inband) + '  ' + str(hbg) + '   ' + str(mjd[i]) + '  ' + str(norm[i]) + '\n')
xbgsub_fffphot.close()
"""


"""
a,b,c,d,inband,mjd = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.bgsub_cl50_g2020norm.fffphot',dtype='float64',usecols=(5,6,7,8,9,11),unpack=True)
a_x,b_x,c_x,d_x,inband_x,mjd_x = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.xbgsub_cl50_g2020norm.fffphot',dtype='float64',usecols=(5,6,7,8,9,11),unpack=True)
plt.figure(1)
plt.plot(mjd,a,'rx-')
plt.plot(mjd_x,a_x,'bx-')
plt.legend(('3C50 sub','3C50 + X-1 sub'),loc='best')
plt.title('A band (0.4-1 keV)',fontsize=12)
plt.figure(2)
plt.plot(mjd,b,'rx-')
plt.plot(mjd,b_x,'bx-')
plt.legend(('3C50 sub','3C50 + X-1 sub'),loc='best')
plt.title('B band (1-2 keV)',fontsize=12)
plt.figure(3)
plt.plot(mjd,c,'rx-')
plt.plot(mjd,c_x,'bx-')
plt.legend(('3C50 sub','3C50 + X-1 sub'),loc='best')
plt.title('C band (2-4 keV)',fontsize=12)
plt.figure(4)
plt.plot(mjd,inband,'rx-')
plt.plot(mjd,inband_x,'bx-')
plt.legend(('3C50 sub','3C50 + X-1 sub'),loc='best')
plt.title('In-band (0.4-12 keV)',fontsize=12)

plt.show()
"""

"""
bgsub = sorted(glob.glob(Lv0_dirs.NGC300_2020 + 'bgsub_cl50/pha/bgsub_cl50*pha'))
comparison = open('/Volumes/Samsung_T5/n300_ulx_2020/n300_ulx.bgsub_cl50_g2020norm.phot','r').read().split('\n')[:-1]
s1_normphot = np.array([np.float(comparison[i].split()[3]) for i in range(len(comparison))])
s2_normphot = np.array([np.float(comparison[i].split()[4]) for i in range(len(comparison))])
a_normphot = np.array([np.float(comparison[i].split()[5]) for i in range(len(comparison))])
b_normphot = np.array([np.float(comparison[i].split()[6]) for i in range(len(comparison))])
c_normphot = np.array([np.float(comparison[i].split()[7]) for i in range(len(comparison))])
d_normphot = np.array([np.float(comparison[i].split()[8]) for i in range(len(comparison))])
inband_normphot = np.array([np.float(comparison[i].split()[9]) for i in range(len(comparison))])
hbg_normphot = np.array([np.float(comparison[i].split()[10]) for i in range(len(comparison))])
norm = np.array([np.float(comparison[i].split()[12]) for i in range(len(comparison))])
"""

"""
##### To print the values from the data set
for i in range(len(bgsub)):
    bgsub_file = str(pathlib.Path(bgsub[i]).name)
    rates = fits.open(bgsub[i])[1].data['RATE']
    chans = fits.open(bgsub[i])[1].data['CHANNEL']
    s1 = str(round(sum(rates[20-1:30-1])*norm[i],4))
    s2 = str(round(sum(rates[30-1:40-1])*norm[i],4))
    a = str(round(sum(rates[40-1:100-1])*norm[i],4))
    b = str(round(sum(rates[100-1:200-1])*norm[i],4))
    c = str(round(sum(rates[200-1:400-1])*norm[i],4))
    d = str(round(sum(rates[400-1:1200-1])*norm[i],4))
    inband = str(round(sum(rates[40-1:1200-1])*norm[i],4))
    hbg = str(round(sum(rates[1300-1:1501-1])*norm[i],4))
    print(bgsub_file,s1,s2,a,b,c,d,inband,hbg)
#####
"""

"""
##### For the errors:
for i in range(len(bgsub)):
    bgsub_file = str(pathlib.Path(bgsub[i]).name)
    errs = fits.open(bgsub[i])[1].data['STAT_ERR']
    chans = fits.open(bgsub[i])[1].data['CHANNEL']
    s1 = str(round( np.sqrt(np.mean(errs[20-1:30-1]**2) * len(errs[20-1:30-1]))*norm[i],4))
    s2 = str(round( np.sqrt(np.mean(errs[30-1:40-1]**2) * len(errs[30-1:40-1]))*norm[i],4))
    a = str(round( np.sqrt(np.mean(errs[40-1:100-1]**2) * len(errs[40-1:100-1]))*norm[i],4))
    b = str(round( np.sqrt(np.mean(errs[100-1:200-1]**2) * len(errs[100-1:200-1]))*norm[i],4))
    c = str(round( np.sqrt(np.mean(errs[200-1:400-1]**2) * len(errs[200-1:400-1]))*norm[i],4))
    d = str(round( np.sqrt(np.mean(errs[400-1:1200-1]**2) * len(errs[400-1:1200-1]))*norm[i],4))
    inband = str(round( np.sqrt(np.mean(errs[40-1:1200-1]**2) * len(errs[40-1:1200-1]))*norm[i],4))
    hbg = str(round( np.sqrt(np.mean(errs[1300-1:1501-1]**2) * len(errs[1300-1:1501-1]))*norm[i],4))
    #print(bgsub_file,s1,s2,a,b,c,d,inband,hbg)
"""

"""
comparison_err = open('/Volumes/Samsung_T5/n300_ulx_2020/n300_ulx.bgsub_cl50_g2020err_norm.phot','r').read().split('\n')[:-1]
s1_err = np.array([np.float(comparison_err[i].split()[3]) for i in range(len(comparison_err))])
s2_err = np.array([np.float(comparison_err[i].split()[4]) for i in range(len(comparison_err))])
a_err = np.array([np.float(comparison_err[i].split()[5]) for i in range(len(comparison_err))])
b_err = np.array([np.float(comparison_err[i].split()[6]) for i in range(len(comparison_err))])
c_err = np.array([np.float(comparison_err[i].split()[7]) for i in range(len(comparison_err))])
d_err = np.array([np.float(comparison_err[i].split()[8]) for i in range(len(comparison_err))])
inband_err = np.array([np.float(comparison_err[i].split()[9]) for i in range(len(comparison_err))])
hbg_err = np.array([np.float(comparison_err[i].split()[10]) for i in range(len(comparison_err))])
norm = np.array([np.float(comparison_err[i].split()[12]) for i in range(len(comparison_err))])

##### To compare between what I see in the bgsub spectra and the *err*.phot files generated by Ron
err_file = open(Lv0_dirs.NGC300_2020 + 'compare_photerr_bgsub.txt','w')
for i in tqdm(range(len(bgsub))):
    bgsub_file = str(pathlib.Path(bgsub[i]).name)
    errs = fits.open(bgsub[i])[1].data['STAT_ERR']
    chans = fits.open(bgsub[i])[1].data['CHANNEL']
    s1 = str(round( (round( np.sqrt(np.mean(errs[20-1:30-1]**2) * len(errs[20-1:30-1]))*norm[i],4) - s1_err[i])/s1_err[i]*100,4 ))
    s2 = str(round( (round( np.sqrt(np.mean(errs[30-1:40-1]**2) * len(errs[30-1:40-1]))*norm[i],4) - s2_err[i])/s2_err[i]*100,4 ))
    a = str(round( (round( np.sqrt(np.mean(errs[40-1:100-1]**2) * len(errs[40-1:100-1]))*norm[i],4) - a_err[i])/a_err[i]*100,4 ))
    b = str(round( (round( np.sqrt(np.mean(errs[100-1:200-1]**2) * len(errs[100-1:200-1]))*norm[i],4) - b_err[i])/b_err[i]*100, 4))
    c = str(round( (round( np.sqrt(np.mean(errs[200-1:400-1]**2) * len(errs[200-1:400-1]))*norm[i],4) - c_err[i])/c_err[i]*100, 4))
    d = str(round( (round( np.sqrt(np.mean(errs[400-1:1200-1]**2) * len(errs[400-1:1200-1]))*norm[i],4) - d_err[i])/d_err[i]*100, 4))
    inband = str(round( (round( np.sqrt(np.mean(errs[40-1:1200-1]**2) * len(errs[40-1:1200-1]))*norm[i],4) - inband_err[i])/inband_err[i]*100, 4))
    hbg = str(round( (round( np.sqrt(np.mean(errs[1300-1:1501-1]**2) * len(errs[1300-1:1501-1]))*norm[i],4) - hbg_err[i])/hbg_err[i]*100, 4))
    err_file.write(bgsub_file + ' ' + s1 + ' ' + s2 + ' ' + a + ' ' + b + ' ' + c + ' ' + d + ' ' + inband + ' ' + hbg + '\n')
err_file.close()
"""

"""
##### To compare between what I see in the bgsub spectra and the .phot files generated by Ron
for i in range(len(bgsub)):
    bgsub_file = str(pathlib.Path(bgsub[i]).name)
    rates = fits.open(bgsub[i])[1].data['RATE']
    chans = fits.open(bgsub[i])[1].data['CHANNEL']
    s1 = str(round((round(sum(rates[20-1:30-1])*norm[i],4)-s1_normphot[i])/s1_normphot[i]*100,4))
    s2 = str(round((round(sum(rates[30-1:40-1])*norm[i],4)-s2_normphot[i])/s2_normphot[i]*100,4))
    a = str(round((round(sum(rates[40-1:100-1])*norm[i],4)-a_normphot[i])/a_normphot[i]*100,4))
    b = str(round((round(sum(rates[100-1:200-1])*norm[i],4)-b_normphot[i])/b_normphot[i]*100,4))
    c = str(round((round(sum(rates[200-1:400-1])*norm[i],4)-c_normphot[i])/c_normphot[i]*100,4))
    d = str(round((round(sum(rates[400-1:1200-1])*norm[i],4)-d_normphot[i])/d_normphot[i]*100,4))
    inband = str(round((round(sum(rates[40-1:1200-1])*norm[i],4)-inband_normphot[i])/inband_normphot[i]*100,4))
    hbg = str(round((round(sum(rates[1300-1:1501-1])*norm[i],4)-hbg_normphot[i])/hbg_normphot[i]*100,4))
    #print(bgsub_file,s1,s2,a,b,c,d,inband,hbg)
"""

"""
##### To insert the error calculations into xbgsub:
bgsub_fffphot = open('/Volumes/Samsung_T5/n300_ulx_2020/n300_ulx.bgsub_cl50_g2020err_norm.fffphot','r').read().split('\n')[:-1]
exp,telapse,mjd,norm = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.bgsub_cl50_g2020err_norm.fffphot',dtype='float64',usecols=(1,2,11,12),unpack=True)
xbgsub_fffphot = np.array([Lv0_dirs.NGC300_2020 + 'bgsub_cl50/xspha/xs' + str(bgsub_fffphot[i].split()[0]) for i in range(len(bgsub_fffphot))])

##### To compare between what I see in the bgsub spectra and the *err*.phot files generated by Ron
xbgsub_err = open(Lv0_dirs.NGC300_2020 + 'n300_ulx.xsbgsub_cl50_g2020err_norm.fffphot','w')
for i in tqdm(range(len(xbgsub_fffphot))):
    xbgsub_file = str(pathlib.Path(xbgsub_fffphot[i]).name)
    errs = fits.open(xbgsub_fffphot[i])[1].data['STAT_ERR']
    chans = fits.open(xbgsub_fffphot[i])[1].data['CHANNEL']
    s1 = str(round( np.sqrt(np.mean(errs[20-1:30-1]**2) * len(errs[20-1:30-1]))*norm[i],4))
    s2 = str(round( np.sqrt(np.mean(errs[30-1:40-1]**2) * len(errs[30-1:40-1]))*norm[i],4))
    a = str(round( np.sqrt(np.mean(errs[40-1:100-1]**2) * len(errs[40-1:100-1]))*norm[i],4))
    b = str(round( np.sqrt(np.mean(errs[100-1:200-1]**2) * len(errs[100-1:200-1]))*norm[i],4))
    c = str(round( np.sqrt(np.mean(errs[200-1:400-1]**2) * len(errs[200-1:400-1]))*norm[i],4))
    d = str(round( np.sqrt(np.mean(errs[400-1:1200-1]**2) * len(errs[400-1:1200-1]))*norm[i],4))
    inband = str(round( np.sqrt(np.mean(errs[40-1:1200-1]**2) * len(errs[40-1:1200-1]))*norm[i],4))
    hbg = str(round( np.sqrt(np.mean(errs[1300-1:1501-1]**2) * len(errs[1300-1:1501-1]))*norm[i],4))
    xbgsub_err.write(xbgsub_file + '  ' + str(exp[i]) + '  ' + str(telapse[i]) + '     ' + str(s1) + ' ' + str(s2) + '      ' + str(a) + ' ' + str(b) + ' ' + str(c) + ' ' + str(d) + '      ' + str(inband) + '  ' + str(hbg) + '   ' + str(mjd[i]) + '  ' + str(norm[i]) + '\n')
xbgsub_err.close()
"""

##### 10/8

def get_par(par_txt,pars):
    """
    Extract values from the parameter text file

    par_txt - path to the text file
    pars - list of parameters
    """
    par_dict = {}
    for i in range(len(pars)):
        par_dict[pars[i]] = []

    partext = open(par_txt,'r').read().split('\n')
    for i in range(len(pars)):
        for j in range(len(partext)):
            if pars[i] in partext[j]:
                parline = partext[j].split()
                par_dict[pars[i]].append((float(parline[-3]),float(parline[-1])))

    return par_dict

def get_Cstat(cstat_txt):
    """
    Extract values from the C-stat text file

    cstat_txt - path to the text file
    """
    cstattxt = open(cstat_txt,'r').read().split('\n')
    cstat = []
    for i in range(len(cstattxt)):
        if "C-Statistic" in cstattxt[i]:
            cstatline = cstattxt[i].split()
            cstat.append((float(cstatline[1]),float(cstatline[3])))
        elif "Total" in cstattxt[i]:
            total = (float(cstattxt[i].split()[-4]),float(cstattxt[i].split()[-2]))

    return cstat,total

def get_lumin(lumin_txt):
    """
    Extract values from the lumin text file

    lumin_txt - path to the text file
    """
    lumintxt = open(lumin_txt,'r').read().split('\n')
    lumin = []
    for i in range(len(lumintxt)):
        if "Luminosity" in lumintxt[i]:
            luminline = lumintxt[i].split()
            lumin.append(float(luminline[2]))

    return np.array(lumin)

def get_flux(flux_txt):
    """
    Extract values from the flux text file

    flux_txt - path to the text file
    """
    fluxtxt = open(flux_txt,'r').read().split('\n')
    flux = []
    for i in range(len(fluxtxt)):
        if "Model" in fluxtxt[i]:
            fluxline = fluxtxt[i].split()
            flux.append(float(fluxline[4][1:]))

    return np.array(flux)

def area_under_curve():
    ################ 9/30/2019 - to implement the 'area under the curve' idea...
    #mjds = ['58239','58244','58249','58254','58259','58264','58269','58274','58279',
    #        '58284','58289','58294','58309','58314','58324','58329','58334','58339',
    #        '58344','58349','58384','58389','58394','58399','58409','58449','58454',
    #        '58459','58464','58469','58474','58479','58484','58489','58494','58499',
    #        '58504','58509','58514','58599','58604']

    """
    ### 10/3 - spectra to actually use!
    mjds = ['58239','58244','58249','58254','58259','58264','58269','58274',
            '58289','58309','58314','58324','58329','58334','58339','58389',
            '58399','58449','58454','58459','58464','58484','58489','58504',
            '58509']

    binspectra_cl = [binspectra_cl[i] for i in range(len(binspectra_cl)) if binspectra_cl[i][-21:-16] in mjds]

    index1 = 0
    index2 = -1
    all_energies_cl = []
    all_cl = []

    all_energies_bgsub = []
    all_bgsub = []
    for i in tqdm(range(len(binspectra_cl))):#[index1:index2]):
        area_cl = 0
        area_cl_percent = []

        area_bgsub = 0
        area_bgsub_percent = []
        E_cl,E_cl_unc,flux_cl,flux_cl_unc = np.genfromtxt(binspectra_cl[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)
        #E_bg,E_bg_unc,flux_bg,flux_bg_unc = np.genfromtxt(binspectra_bg[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)
        E_bgsub,E_bgsub_unc,flux_bgsub,flux_bgsub_unc = np.genfromtxt(binspectra_bgsub[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)

        for j in range(len(E_cl)-1):
            if flux_cl[j] >= 0:
                area_cl += (E_cl[j+1]-E_cl[j])*flux_cl[j]

        cumulative_cl = 0
        for j in range(len(E_cl)-1):
            if flux_cl[j] >= 0:
                area_current_cl = (E_cl[j+1]-E_cl[j])*flux_cl[j]
                cumulative_cl += area_current_cl/area_cl*100
                area_cl_percent.append(cumulative_cl)
            else:
                if len(area_cl_percent) == 0:
                    area_cl_percent.append(0)
                else:
                    area_cl_percent.append(area_cl_percent[-1])

        for j in range(len(E_bgsub)-1):
            if flux_bgsub[j] >= 0:
                area_bgsub += (E_bgsub[j+1]-E_bgsub[j])*flux_bgsub[j]

        cumulative_bgsub = 0
        for j in range(len(E_bgsub)-1):
            if flux_bgsub[j] >= 0:
                area_current_bgsub = (E_bgsub[j+1]-E_bgsub[j])*flux_bgsub[j]
                cumulative_bgsub += area_current_bgsub/area_bgsub*100
                area_bgsub_percent.append(cumulative_bgsub)
            else:
                if len(area_bgsub_percent) == 0:
                    area_bgsub_percent.append(0)
                else:
                    area_bgsub_percent.append(area_bgsub_percent[-1])

        all_energies_cl.append(E_cl[:-1])
        all_cl.append(area_cl_percent)

        all_energies_bgsub.append(E_bgsub[:-1])
        all_bgsub.append(area_bgsub_percent)


    plt.figure(1)
    for i in range(len(all_energies_cl)):
        #plt.figure(i+100)
        plt.plot(all_energies_cl[i],all_cl[i],label=mjds[i+index1])
    #plt.plot(all_energies_cl[-3],all_cl[-3])
    plt.xlabel('Energy, E (keV)',fontsize=12)
    plt.ylabel('Percentage of total area under curve, %',fontsize=12)
    datacursor(formatter='{label}'.format,bbox=None)

    plt.axhline(y=80,alpha=0.5,lw=0.5)
    plt.axhline(y=90,alpha=0.5,lw=0.5)
    plt.axvline(x=4,alpha=0.5,lw=0.5)
    plt.axvline(x=5,alpha=0.5,lw=0.5)
    plt.axvline(x=6,alpha=0.5,lw=0.5)

    plt.show()
    """

    """
    plt.figure(1)
    for i in range(len(all_energies_bgsub)):
        plt.plot(all_energies_bgsub[i],all_bgsub[i],label=mjds[i])
    plt.xlabel('Energy, E (keV)',fontsize=12)
    plt.ylabel('Percentage of total area under curve, %',fontsize=12)
    datacursor(formatter='{label}'.format)
    plt.show()
    """

    """
    #bound1 = [0,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    #bound2 = [0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    bound1 = [0.0,0.5,0.6,0.7,1.0,2.0,3.0,4.0,4.3,4.6,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0]
    bound2 = [0.5,0.6,0.7,1.0,2.0,3.0,4.0,4.3,4.6,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0]
    mjds = []
    for i in range(len(binspectra_bgsub)):
        if binspectra_bgsub[i][-24:-19] not in not_to_use:
            mjd = binspectra_bgsub[i][-24:-19]
            E,E_unc,flux,flux_unc = np.genfromtxt(binspectra_bgsub[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)
            total_E += list(E)

            trunc1 = E[(E>=bound1[0])&(E<bound2[0])]
            trunc2 = E[(E>=bound1[1])&(E<bound2[1])]
            trunc3 = E[(E>=bound1[2])&(E<bound2[2])]
            trunc4 = E[(E>=bound1[3])&(E<bound2[3])]
            trunc5 = E[(E>=bound1[4])&(E<bound2[4])]
            trunc6 = E[(E>=bound1[5])&(E<bound2[5])]
            trunc7 = E[(E>=bound1[6])&(E<bound2[6])]
            trunc8 = E[(E>=bound1[7])&(E<bound2[7])]
            trunc9 = E[(E>=bound1[8])&(E<bound2[8])]
            trunc10 = E[(E>=bound1[9])&(E<bound2[9])]
            trunc11 = E[(E>=bound1[10])&(E<bound2[10])]
            trunc12 = E[(E>=bound1[11])&(E<bound2[11])]
            trunc13 = E[(E>=bound1[12])&(E<bound2[12])]
            trunc14 = E[(E>=bound1[13])&(E<bound2[13])]
            trunc15 = E[(E>=bound1[14])&(E<bound2[14])]
            trunc16 = E[(E>=bound1[15])&(E<bound2[15])]
            trunc17 = E[(E>=bound1[16])&(E<bound2[16])]
            trunc18 = E[(E>=bound1[17])&(E<bound2[17])]
            trunc19 = E[(E>=bound1[18])&(E<bound2[18])]
            trunc20 = E[(E>=bound1[19])&(E<bound2[19])]
            print(str(mjd) + ' ' + str(len(trunc1)) + ' ' + str(len(trunc2))
            + ' ' + str(len(trunc3)) + ' ' + str(len(trunc4)) + ' ' + str(len(trunc5))
            + ' ' + str(len(trunc6)) + ' ' + str(len(trunc7)) + ' ' + str(len(trunc8))
            + ' ' + str(len(trunc9)) + ' ' + str(len(trunc10)) + ' ' + str(len(trunc11))
            + ' ' + str(len(trunc12)) + ' ' + str(len(trunc13)) + ' ' + str(len(trunc14))
            + ' ' + str(len(trunc15)) + ' ' + str(len(trunc16)) + ' ' + str(len(trunc17))
            + ' ' + str(len(trunc18)) + ' ' + str(len(trunc19)) + ' ' + str(len(trunc20)))

    total_E = np.array(total_E)


    cumulative_percent = 0
    for i in range(len(bound1)):
        truncated = total_E[(total_E>=bound1[i])&(total_E<=bound2[i])]
        percentage = round(len(truncated)/len(total_E)*100,2)
        cumulative_percent += percentage
        print('The percentage of ' + str(bound1[i]) + '-' + str(bound2[i]) + ' keV out of all data is: ' + str(percentage) + '% ; cumulative: ' + str(cumulative_percent) + '%')


    plt.hist(total_E,bins=200,log=True)
    plt.xlabel('Energy, E (keV)',fontsize=12)
    plt.ylabel('Number',fontsize=12)
    plt.show()
    """

##### 10/26; data obtained from doing a more systematic set of fits using tcl scripts

def plot_fit_table(tablefile,pars):
    """
    Plot the results from the fit_$MODEL.table tables

    table - the path to the table, which includes the model name
    pars - model parameters, in the form of a list
    """
    par_dict = {}
    error_dict = {}
    sigma_dict = {}

    filenames = np.genfromtxt(tablefile,dtype='str',usecols=(0),unpack=True)
    mjds = np.array([float(filenames[i][-25:-20]) for i in range(len(filenames))])
    no_cols = 1 + 4*len(pars) + 8

    counter = 0 #to keep track of the variable index; of the same order as in "pars"!
    for i in range(1,no_cols-8,4): #for the parameter value
        data_array = np.genfromtxt(tablefile,dtype='float',usecols=(i),unpack=True)
        par_dict[pars[counter]] = data_array
        counter += 1
    counter = 0 #to keep track of the variable index
    for i in range(2,no_cols-8,4): #for the sigma value
        data_array = np.genfromtxt(tablefile,dtype='float',usecols=(i),unpack=True)
        sigma_dict[pars[counter]] = data_array
        counter += 1
    counter = 0 #to keep track of the variable index
    for i in range(3,no_cols-8,4): #for the error value
        error_low, error_high = np.genfromtxt(tablefile,dtype='float',usecols=(i,i+1),unpack=True) #two arrays for lower and upper bound!
        error_dict[pars[counter]] = np.array([(error_low[j],error_high[j]) for j in range(len(error_low)) ])
        counter += 1

    cstat = np.genfromtxt(tablefile,dtype='float',usecols=(no_cols-8),unpack=True)
    dof = np.genfromtxt(tablefile,dtype='float',usecols=(no_cols-7),unpack=True)
    flux,flux_low,flux_high = np.genfromtxt(tablefile,dtype='float',usecols=(no_cols-6,no_cols-5,no_cols-4),unpack=True)
    lumin,lumin_low,lumin_high = np.genfromtxt(tablefile,dtype='float',usecols=(no_cols-3,no_cols-2,no_cols-1),unpack=True)*1e44

    ##### Plotting
    fig,axs = plt.subplots(len(pars)+3,sharex=True,gridspec_kw={'hspace':0})
    fig.suptitle('Total C-stat: ' + str(round(sum(cstat),4)) + '/' + str(sum(dof)))

    for i in range(len(pars)):
        axs[i].errorbar(x=mjds,y=par_dict[pars[i]],yerr=sigma_dict[pars[i]],color='r',fmt='x')
        axs[i].set_ylabel(pars[i])
    axs[2].set_yscale('log')
    axs[len(pars)].plot(mjds,lumin/1e39,'rx')
    axs[len(pars)].set_yscale('log')
    axs[len(pars)+1].plot(mjds,flux/1e-12,'rx')
    axs[len(pars)+1].set_yscale('log')
    axs[len(pars)+2].plot(mjds,cstat/dof,'rx')
    axs[len(pars)+2].axhline(y=1.0,lw=0.5,alpha=0.5)

    axs[len(pars)].set_ylabel('Luminosity \n 1e39 \n ergs/s)')
    axs[len(pars)+1].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
    axs[len(pars)+2].set_ylabel('Reduced \n C-stat')

    axs[len(pars)+2].set_xlabel('Time (MJD)',fontsize=12)

    plt.show()


if __name__ == "__main__":

    tablefile = Lv0_dirs.NGC300_2020 + 'spectra_05d/fit_tbnew-diskbb.table'
    pars = ['nH','Tin','norm']
    plot_fit_table(tablefile,pars)

    ##### For POWERLAW

    def powerlaw():
        cstat_varyPL,total_cstat_varyPL = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varyPL_cstat.txt')
        lumin_varyPL = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varyPL_lumin.txt')
        flux_varyPL = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varyPL_flux.txt')

        par_dict_varynH = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varynH_par.txt',['nH','PhoIndex','norm'])
        cstat_varynH,total_cstat_varynH = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varynH_cstat.txt')
        lumin_varynH = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varynH_lumin.txt')
        flux_varynH = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varynH_flux.txt')

        par_dict_varyall = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varyall_par.txt',['nH','PhoIndex','norm'])
        cstat_varyall,total_cstat_varyall = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varyall_cstat.txt')
        lumin_varyall = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varyall_lumin.txt')
        flux_varyall = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-powerlaw_varyall_flux.txt')

        mjds = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.xbgsub_cl50_g2020norm_05d.fffphot',usecols=(0),unpack=True)

        ##### for varyPL
        def varyPL():
            fig,axs = plt.subplots(5,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep nH (TBNEW) fixed at (' + str(par_dict_varyPL['nH'][0][0]) + ' +- ' + str(par_dict_varyPL['nH'][0][1]) + ') E22 atoms cm^-2' + '\n' + 'Total C-stat: ' + str(total_cstat_varyPL[0]) + '/' + str(total_cstat_varyPL[1]))
            phoindex_varyPL = np.array([par_dict_varyPL['PhoIndex'][i][0] for i in range(len(par_dict_varyPL['PhoIndex']))])
            phoindex_err_varyPL = np.array([par_dict_varyPL['PhoIndex'][i][1] for i in range(len(par_dict_varyPL['PhoIndex']))])
            norm_varyPL = np.array([par_dict_varyPL['norm'][i][0] for i in range(len(par_dict_varyPL['norm']))])
            norm_err_varyPL = np.array([par_dict_varyPL['norm'][i][1] for i in range(len(par_dict_varyPL['norm']))])

            cstatarr_varyPL = np.array([cstat_varyPL[i][0]/cstat_varyPL[i][1] for i in range(len(cstat_varyPL)) ])

            axs[0].errorbar(x=mjds,y=phoindex_varyPL,yerr=phoindex_err_varyPL,color='r',fmt='x')
            axs[1].errorbar(x=mjds,y=norm_varyPL/1e-4,yerr=norm_err_varyPL/1e-4,color='r',fmt='x')
            axs[2].plot(mjds,lumin_varyPL/1e38,'rx')
            axs[3].plot(mjds,flux_varyPL/1e-12,'rx')
            axs[4].plot(mjds,cstatarr_varyPL,'rx')
            axs[4].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('PhoIndex')
            axs[1].set_ylabel('Norm (1e-4)')
            axs[2].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[3].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[4].set_ylabel('Reduced \n C-stat')

            axs[4].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        ##### for varynH
        def varynH():
            fig,axs = plt.subplots(4,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep PhoIndex fixed at (' + str(par_dict_varynH['PhoIndex'][0][0]) + ' +- ' + str(par_dict_varynH['PhoIndex'][0][1]) + ')' + '\n' + 'Keep norm fixed at (' + str(par_dict_varynH['norm'][0][0]) + ' +- ' + str(par_dict_varynH['norm'][0][1]) + '\n' + 'Total C-stat: ' + str(total_cstat_varynH[0]) + '/' + str(total_cstat_varynH[1]))
            nH_varynH = np.array([par_dict_varynH['nH'][i][0] for i in range(len(par_dict_varynH['nH']))])
            nH_err_varynH = np.array([par_dict_varynH['nH'][i][1] for i in range(len(par_dict_varynH['nH']))])

            cstatarr_varynH = np.array([cstat_varynH[i][0]/cstat_varynH[i][1] for i in range(len(cstat_varynH)) ])

            axs[0].errorbar(x=mjds[nH_varynH<100],y=nH_varynH[nH_varynH<100],yerr=nH_err_varynH[nH_varynH<100],color='r',fmt='x')
            axs[0].set_yscale('log')
            axs[1].plot(mjds,lumin_varynH/1e38,'rx')
            axs[2].plot(mjds,flux_varynH/1e-12,'rx')
            axs[3].plot(mjds,cstatarr_varynH,'rx')
            axs[3].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('nH (E22 atoms/cm^2)')
            axs[1].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[2].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[3].set_ylabel('Reduced \n C-stat')

            axs[3].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        ##### for varyall
        def varyall():
            fig,axs = plt.subplots(6,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('All of nH, PhoIndex, and norm are untied' + '\n' + 'Total C-stat: ' + str(total_cstat_varyall[0]) + '/' + str(total_cstat_varyall[1]))
            phoindex_varyall = np.array([par_dict_varyall['PhoIndex'][i][0] for i in range(len(par_dict_varyall['PhoIndex']))])
            phoindex_err_varyall = np.array([par_dict_varyall['PhoIndex'][i][1] for i in range(len(par_dict_varyall['PhoIndex']))])
            norm_varyall = np.array([par_dict_varyall['norm'][i][0] for i in range(len(par_dict_varyall['norm']))])
            norm_err_varyall = np.array([par_dict_varyall['norm'][i][1] for i in range(len(par_dict_varyall['norm']))])
            nH_varyall = np.array([par_dict_varyall['nH'][i][0] for i in range(len(par_dict_varyall['nH']))])
            nH_err_varyall = np.array([par_dict_varyall['nH'][i][1] for i in range(len(par_dict_varyall['nH']))])

            cstatarr_varyall = np.array([cstat_varyall[i][0]/cstat_varyall[i][1] for i in range(len(cstat_varyall)) ])

            axs[0].errorbar(x=mjds,y=phoindex_varyall,yerr=phoindex_err_varyall,color='r',fmt='x')
            axs[1].errorbar(x=mjds,y=norm_varyall/1e-4,yerr=norm_err_varyall/1e-4,color='r',fmt='x')
            axs[2].errorbar(x=mjds,y=nH_varyall,yerr=nH_err_varyall,color='r',fmt='x')
            axs[3].plot(mjds,lumin_varyall/1e38,'rx')
            axs[4].plot(mjds,flux_varyall/1e-12,'rx')
            axs[5].plot(mjds,cstatarr_varyall,'rx')
            axs[5].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('PhoIndex')
            axs[1].set_ylabel('Norm \n (1e-4)')
            axs[2].set_ylabel('nH \n 1E22 \n atoms/cm^2')
            axs[3].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[4].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[5].set_ylabel('Reduced \n C-stat')

            axs[5].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

    def bbodyrad():
        ##### vary kT, norm:
        par_dict_varyBB = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varyBB_par.txt',['nH','kT','norm'])
        cstat_varyBB,total_cstat_varyBB = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varyBB_cstat.txt')
        lumin_varyBB = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varyBB_lumin.txt')
        flux_varyBB = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varyBB_flux.txt')

        par_dict_varynH = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varynH_par.txt',['nH','kT','norm'])
        cstat_varynH,total_cstat_varynH = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varynH_cstat.txt')
        lumin_varynH = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varynH_lumin.txt')
        flux_varynH = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varynH_flux.txt')

        par_dict_varyall = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varyall_par.txt',['nH','kT','norm'])
        cstat_varyall,total_cstat_varyall = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varyall_cstat.txt')
        lumin_varyall = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varyall_lumin.txt')
        flux_varyall = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-bbodyrad_varyall_flux.txt')

        mjds = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.xbgsub_cl50_g2020norm_05d.fffphot',usecols=(0),unpack=True)

        ##### for varyBB
        def varyBB():
            fig,axs = plt.subplots(5,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep nH (TBNEW) fixed at (' + str(par_dict_varyBB['nH'][0][0]) + ' +- ' + str(par_dict_varyBB['nH'][0][1]) + ') E22 atoms cm^-2' + '\n' + 'Total C-stat: ' + str(total_cstat_varyBB[0]) + '/' + str(total_cstat_varyBB[1]))
            kT_varyBB = np.array([par_dict_varyBB['kT'][i][0] for i in range(len(par_dict_varyBB['kT']))])
            kT_err_varyBB = np.array([par_dict_varyBB['kT'][i][1] for i in range(len(par_dict_varyBB['kT']))])
            norm_varyBB = np.array([par_dict_varyBB['norm'][i][0] for i in range(len(par_dict_varyBB['norm']))])
            norm_err_varyBB = np.array([par_dict_varyBB['norm'][i][1] for i in range(len(par_dict_varyBB['norm']))])

            cstatarr_varyBB = np.array([cstat_varyBB[i][0]/cstat_varyBB[i][1] for i in range(len(cstat_varyBB)) ])

            axs[0].errorbar(x=mjds,y=kT_varyBB,yerr=kT_err_varyBB,color='r',fmt='x')
            axs[1].errorbar(x=mjds,y=norm_varyBB,yerr=norm_err_varyBB,color='r',fmt='x')
            axs[2].plot(mjds,lumin_varyBB/1e38,'rx')
            axs[3].plot(mjds,flux_varyBB/1e-12,'rx')
            axs[4].plot(mjds,cstatarr_varyBB,'rx')
            axs[4].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('kT (keV)')
            axs[1].set_ylabel('Norm')
            axs[2].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[3].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[4].set_ylabel('Reduced \n C-stat')

            axs[4].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        ##### for varynH
        def varynH():
            fig,axs = plt.subplots(4,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep kT fixed at (' + str(par_dict_varynH['kT'][0][0]) + ' +- ' + str(par_dict_varynH['kT'][0][1]) + ') keV' + '\n' + 'Keep norm fixed at (' + str(par_dict_varynH['norm'][0][0]) + ' +- ' + str(par_dict_varynH['norm'][0][1]) + '\n' + 'Total C-stat: ' + str(total_cstat_varynH[0]) + '/' + str(total_cstat_varynH[1]))
            nH_varynH = np.array([par_dict_varynH['nH'][i][0] for i in range(len(par_dict_varynH['nH']))])
            nH_err_varynH = np.array([par_dict_varynH['nH'][i][1] for i in range(len(par_dict_varynH['nH']))])

            cstatarr_varynH = np.array([cstat_varynH[i][0]/cstat_varynH[i][1] for i in range(len(cstat_varynH)) ])

            axs[0].errorbar(x=mjds,y=nH_varynH,yerr=nH_err_varynH,color='r',fmt='x')
            axs[0].set_yscale('log')
            axs[1].plot(mjds,lumin_varynH/1e38,'rx')
            axs[2].plot(mjds,flux_varynH/1e-12,'rx')
            axs[3].plot(mjds,cstatarr_varynH,'rx')
            axs[3].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('nH (E22 atoms/cm^2)')
            axs[1].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[2].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[3].set_ylabel('Reduced \n C-stat')

            axs[3].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        ##### for varyall
        def varyall():
            fig,axs = plt.subplots(6,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('All of nH, kT, and norm are untied' + '\n' + 'Total C-stat: ' + str(total_cstat_varyall[0]) + '/' + str(total_cstat_varyall[1]))
            kT_varyall = np.array([par_dict_varyall['kT'][i][0] for i in range(len(par_dict_varyall['kT']))])
            kT_err_varyall = np.array([par_dict_varyall['kT'][i][1] for i in range(len(par_dict_varyall['kT']))])
            norm_varyall = np.array([par_dict_varyall['norm'][i][0] for i in range(len(par_dict_varyall['norm']))])
            norm_err_varyall = np.array([par_dict_varyall['norm'][i][1] for i in range(len(par_dict_varyall['norm']))])
            nH_varyall = np.array([par_dict_varyall['nH'][i][0] for i in range(len(par_dict_varyall['nH']))])
            nH_err_varyall = np.array([par_dict_varyall['nH'][i][1] for i in range(len(par_dict_varyall['nH']))])

            cstatarr_varyall = np.array([cstat_varyall[i][0]/cstat_varyall[i][1] for i in range(len(cstat_varyall)) ])

            axs[0].errorbar(x=mjds,y=kT_varyall,yerr=kT_err_varyall,color='r',fmt='x')
            axs[1].errorbar(x=mjds,y=norm_varyall,yerr=norm_err_varyall,color='r',fmt='x')
            axs[2].errorbar(x=mjds,y=nH_varyall,yerr=nH_err_varyall,color='r',fmt='x')
            axs[3].plot(mjds,lumin_varyall/1e38,'rx')
            axs[4].plot(mjds,flux_varyall/1e-12,'rx')
            axs[5].plot(mjds,cstatarr_varyall,'rx')
            axs[5].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('kT')
            axs[1].set_ylabel('Norm')
            axs[2].set_ylabel('nH \n 1E22 \n atoms/cm^2')
            axs[3].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[4].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[5].set_ylabel('Reduced \n C-stat')

            axs[5].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

    def diskbb():
        ##### vary kT, norm:
        par_dict_varyBB = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varyBB_par.txt',['nH','Tin','norm'])
        cstat_varyBB,total_cstat_varyBB = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varyBB_cstat.txt')
        lumin_varyBB = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varyBB_lumin.txt')
        flux_varyBB = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varyBB_flux.txt')

        par_dict_varynH = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varynH_par.txt',['nH','Tin','norm'])
        cstat_varynH,total_cstat_varynH = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varynH_cstat.txt')
        lumin_varynH = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varynH_lumin.txt')
        flux_varynH = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varynH_flux.txt')

        par_dict_varyall = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varyall_par.txt',['nH','Tin','norm'])
        cstat_varyall,total_cstat_varyall = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varyall_cstat.txt')
        lumin_varyall = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varyall_lumin.txt')
        flux_varyall = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-diskbb_varyall_flux.txt')

        mjds = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.xbgsub_cl50_g2020norm_05d.fffphot',usecols=(0),unpack=True)

        ##### for varyBB
        def varyBB():
            fig,axs = plt.subplots(5,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep nH (TBNEW) fixed at (' + str(par_dict_varyBB['nH'][0][0]) + ' +- ' + str(par_dict_varyBB['nH'][0][1]) + ') E22 atoms cm^-2' + '\n' + 'Total C-stat: ' + str(total_cstat_varyBB[0]) + '/' + str(total_cstat_varyBB[1]))
            Tin_varyBB = np.array([par_dict_varyBB['Tin'][i][0] for i in range(len(par_dict_varyBB['Tin']))])
            Tin_err_varyBB = np.array([par_dict_varyBB['Tin'][i][1] for i in range(len(par_dict_varyBB['Tin']))])
            norm_varyBB = np.array([par_dict_varyBB['norm'][i][0] for i in range(len(par_dict_varyBB['norm']))])
            norm_err_varyBB = np.array([par_dict_varyBB['norm'][i][1] for i in range(len(par_dict_varyBB['norm']))])

            cstatarr_varyBB = np.array([cstat_varyBB[i][0]/cstat_varyBB[i][1] for i in range(len(cstat_varyBB)) ])

            axs[0].errorbar(x=mjds,y=Tin_varyBB,yerr=Tin_err_varyBB,color='r',fmt='x')
            axs[1].errorbar(x=mjds,y=norm_varyBB,yerr=norm_err_varyBB,color='r',fmt='x')
            axs[2].plot(mjds,lumin_varyBB/1e38,'rx')
            axs[3].plot(mjds,flux_varyBB/1e-12,'rx')
            axs[4].plot(mjds,cstatarr_varyBB,'rx')
            axs[4].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('Tin (keV)')
            axs[1].set_ylabel('Norm')
            axs[2].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[3].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[4].set_ylabel('Reduced \n C-stat')

            axs[4].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        ##### for varynH
        def varynH():
            fig,axs = plt.subplots(4,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep Tin fixed at (' + str(par_dict_varynH['Tin'][0][0]) + ' +- ' + str(par_dict_varynH['Tin'][0][1]) + ') keV' + '\n' + 'Keep norm fixed at (' + str(par_dict_varynH['norm'][0][0]) + ' +- ' + str(par_dict_varynH['norm'][0][1]) + '\n' + 'Total C-stat: ' + str(total_cstat_varynH[0]) + '/' + str(total_cstat_varynH[1]))
            nH_varynH = np.array([par_dict_varynH['nH'][i][0] for i in range(len(par_dict_varynH['nH']))])
            nH_err_varynH = np.array([par_dict_varynH['nH'][i][1] for i in range(len(par_dict_varynH['nH']))])

            cstatarr_varynH = np.array([cstat_varynH[i][0]/cstat_varynH[i][1] for i in range(len(cstat_varynH)) ])

            axs[0].errorbar(x=mjds,y=nH_varynH,yerr=nH_err_varynH,color='r',fmt='x')
            axs[0].set_yscale('log')
            axs[1].plot(mjds,lumin_varynH/1e38,'rx')
            axs[2].plot(mjds,flux_varynH/1e-12,'rx')
            axs[3].plot(mjds,cstatarr_varynH,'rx')
            axs[3].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('nH (E22 atoms/cm^2)')
            axs[1].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[2].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[3].set_ylabel('Reduced \n C-stat')

            axs[3].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        ##### for varyall
        def varyall():
            fig,axs = plt.subplots(6,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('All of nH, Tin, and norm are untied' + '\n' + 'Total C-stat: ' + str(total_cstat_varyall[0]) + '/' + str(total_cstat_varyall[1]))
            Tin_varyall = np.array([par_dict_varyall['Tin'][i][0] for i in range(len(par_dict_varyall['Tin']))])
            Tin_err_varyall = np.array([par_dict_varyall['Tin'][i][1] for i in range(len(par_dict_varyall['Tin']))])
            norm_varyall = np.array([par_dict_varyall['norm'][i][0] for i in range(len(par_dict_varyall['norm']))])
            norm_err_varyall = np.array([par_dict_varyall['norm'][i][1] for i in range(len(par_dict_varyall['norm']))])
            nH_varyall = np.array([par_dict_varyall['nH'][i][0] for i in range(len(par_dict_varyall['nH']))])
            nH_err_varyall = np.array([par_dict_varyall['nH'][i][1] for i in range(len(par_dict_varyall['nH']))])

            cstatarr_varyall = np.array([cstat_varyall[i][0]/cstat_varyall[i][1] for i in range(len(cstat_varyall)) ])

            axs[0].errorbar(x=mjds,y=Tin_varyall,yerr=Tin_err_varyall,color='r',fmt='x')
            axs[1].errorbar(x=mjds,y=norm_varyall,yerr=norm_err_varyall,color='r',fmt='x')
            axs[2].errorbar(x=mjds,y=nH_varyall,yerr=nH_err_varyall,color='r',fmt='x')
            axs[3].plot(mjds,lumin_varyall/1e38,'rx')
            axs[4].plot(mjds,flux_varyall/1e-12,'rx')
            axs[5].plot(mjds,cstatarr_varyall,'rx')
            axs[5].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('Tin')
            axs[1].set_ylabel('Norm')
            axs[2].set_ylabel('nH \n 1E22 \n atoms/cm^2')
            axs[3].set_ylabel('Luminosity \n 1e38 \n ergs/s)')
            axs[4].set_ylabel('Flux \n 1e-12 \n ergs/cm^2/s')
            axs[5].set_ylabel('Reduced \n C-stat')

            axs[5].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        #varyall()

    def simpl_diskbb():
        par_dict_varyGamma = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyGamma_par.txt',['nH','Gamma','FracSctr','Tin','norm'])
        cstat_varyGamma,total_cstat_varyGamma = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyGamma_cstat.txt')
        lumin_varyGamma = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyGamma_lumin.txt')
        flux_varyGamma = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyGamma_flux.txt')

        par_dict_varyBB = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB_par.txt',['nH','Gamma','FracSctr','Tin','norm'])
        cstat_varyBB,total_cstat_varyBB = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB_cstat.txt')
        lumin_varyBB = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB_lumin.txt')
        flux_varyBB = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB_flux.txt')

        par_dict_varyBB_Gamma = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB-Gamma_par.txt',['nH','Gamma','FracSctr','Tin','norm'])
        cstat_varyBB_Gamma,total_cstat_varyBB_Gamma = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB-Gamma_cstat.txt')
        lumin_varyBB_Gamma = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB-Gamma_lumin.txt')
        flux_varyBB_Gamma = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB-Gamma_flux.txt')

        par_dict_varyBB_PL = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB-PL_par.txt',['nH','Gamma','FracSctr','Tin','norm'])
        cstat_varyBB_PL,total_cstat_varyBB_PL = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB-PL_cstat.txt')
        lumin_varyBB_PL = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB-PL_lumin.txt')
        flux_varyBB_PL = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/tbnew-simpl-diskbb-varyBB-PL_flux.txt')

        mjds = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.xbgsub_cl50_g2020norm_05d.fffphot',usecols=(0),unpack=True)

        def varyGamma():
            fig,axs = plt.subplots(4,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep nH (TBNEW) fixed at (' + str(par_dict_varyGamma['nH'][0][0]) + ' +- ' + str(par_dict_varyGamma['nH'][0][1]) + ') E22 atoms cm^-2' + '\n' + 'FracSctr fixed at ' + str(par_dict_varyGamma['FracSctr'][0][0]) + ' +- ' + str(par_dict_varyGamma['FracSctr'][0][1])  + '\n' + 'Tin fixed at ' + str(par_dict_varyGamma['Tin'][0][0]) + ' +- ' + str(par_dict_varyGamma['Tin'][0][1]) + '\n' + 'diskbb norm fixed at ' + str(par_dict_varyGamma['norm'][0][0]) + ' +- ' + str(par_dict_varyGamma['norm'][0][1]) + '\n' + 'Total C-stat: ' + str(total_cstat_varyGamma[0]) + '/' + str(total_cstat_varyGamma[1]))
            Gamma_varyGamma = np.array([par_dict_varyGamma['Gamma'][i][0] for i in range(len(par_dict_varyGamma['Gamma']))])
            Gamma_err_varyGamma = np.array([par_dict_varyGamma['Gamma'][i][1] for i in range(len(par_dict_varyGamma['Gamma']))])

            cstatarr_varyGamma = np.array([cstat_varyGamma[i][0]/cstat_varyGamma[i][1] for i in range(len(cstat_varyGamma)) ])

            axs[0].errorbar(x=mjds,y=Gamma_varyGamma,yerr=Gamma_err_varyGamma,color='r',fmt='x')
            axs[0].set_yscale('log')
            axs[1].plot(mjds,lumin_varyGamma,'rx')
            axs[1].set_yscale('log')
            axs[2].plot(mjds,flux_varyGamma,'rx')
            axs[2].set_yscale('log')
            axs[3].plot(mjds,cstatarr_varyGamma,'rx')
            axs[3].set_yscale('log')
            axs[3].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('Gamma')
            axs[1].set_ylabel('Luminosity \n (ergs/s)')
            axs[2].set_ylabel('Flux \n (ergs/cm^2/s)')
            axs[3].set_ylabel('Reduced \n C-stat')

            axs[3].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        def varyBB():
            fig,axs = plt.subplots(5,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep nH (TBNEW) fixed at (' + str(par_dict_varyBB['nH'][0][0]) + ' +- ' + str(par_dict_varyBB['nH'][0][1]) + ') E22 atoms cm^-2' + '\n' + 'FracSctr fixed at ' + str(par_dict_varyBB['FracSctr'][0][0]) + ' +- ' + str(par_dict_varyBB['FracSctr'][0][1])  + '\n' + 'Gamma fixed at ' + str(par_dict_varyBB['Gamma'][0][0]) + ' +- ' + str(par_dict_varyBB['Gamma'][0][1]) + '\n' + 'Total C-stat: ' + str(total_cstat_varyBB[0]) + '/' + str(total_cstat_varyBB[1]))
            Tin_varyBB = np.array([par_dict_varyBB['Tin'][i][0] for i in range(len(par_dict_varyBB['Tin']))])
            Tin_err_varyBB = np.array([par_dict_varyBB['Tin'][i][1] for i in range(len(par_dict_varyBB['Tin']))])
            norm_varyBB = np.array([par_dict_varyBB['norm'][i][0] for i in range(len(par_dict_varyBB['norm']))])
            norm_err_varyBB = np.array([par_dict_varyBB['norm'][i][1] for i in range(len(par_dict_varyBB['norm']))])

            cstatarr_varyBB = np.array([cstat_varyBB[i][0]/cstat_varyBB[i][1] for i in range(len(cstat_varyBB)) ])

            axs[0].errorbar(x=mjds,y=Tin_varyBB,yerr=Tin_err_varyBB,color='r',fmt='x')
            axs[0].set_yscale('log')
            axs[1].errorbar(x=mjds,y=norm_varyBB,yerr=norm_err_varyBB,color='r',fmt='x')
            axs[1].set_yscale('log')
            axs[2].plot(mjds,lumin_varyBB,'rx')
            axs[2].set_yscale('log')
            axs[3].plot(mjds,flux_varyBB,'rx')
            axs[3].set_yscale('log')
            axs[4].plot(mjds,cstatarr_varyBB,'rx')
            axs[4].set_yscale('log')
            axs[4].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('Tin (keV)')
            axs[1].set_ylabel('norm')
            axs[2].set_ylabel('Luminosity \n (ergs/s)')
            axs[3].set_ylabel('Flux \n (ergs/cm^2/s)')
            axs[4].set_ylabel('Reduced \n C-stat')

            axs[4].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        def varyBB_Gamma():
            fig,axs = plt.subplots(6,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep nH (TBNEW) fixed at (' + str(par_dict_varyBB_Gamma['nH'][0][0]) + ' +- ' + str(par_dict_varyBB_Gamma['nH'][0][1]) + ') E22 atoms cm^-2' + '\n' + 'FracSctr fixed at ' + str(par_dict_varyBB_Gamma['FracSctr'][0][0]) + ' +- ' + str(par_dict_varyBB_Gamma['FracSctr'][0][1]) + '\n' + 'Total C-stat: ' + str(total_cstat_varyBB_Gamma[0]) + '/' + str(total_cstat_varyBB_Gamma[1]))
            Gamma_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['Gamma'][i][0] for i in range(len(par_dict_varyBB_Gamma['Gamma']))])
            Gamma_err_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['Gamma'][i][1] for i in range(len(par_dict_varyBB_Gamma['Gamma']))])
            Tin_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['Tin'][i][0] for i in range(len(par_dict_varyBB_Gamma['Tin']))])
            Tin_err_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['Tin'][i][1] for i in range(len(par_dict_varyBB_Gamma['Tin']))])
            norm_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['norm'][i][0] for i in range(len(par_dict_varyBB_Gamma['norm']))])
            norm_err_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['norm'][i][1] for i in range(len(par_dict_varyBB_Gamma['norm']))])

            cstatarr_varyBB_Gamma = np.array([cstat_varyBB_Gamma[i][0]/cstat_varyBB_Gamma[i][1] for i in range(len(cstat_varyBB_Gamma)) ])

            axs[0].errorbar(x=mjds,y=Gamma_varyBB_Gamma,yerr=Gamma_err_varyBB_Gamma,color='r',fmt='x')
            axs[0].axhline(y=4.5,lw=0.5,alpha=0.5)
            axs[0].axhline(y=1,lw=0.5,alpha=0.5)
            axs[1].errorbar(x=mjds,y=Tin_varyBB_Gamma,yerr=Tin_err_varyBB_Gamma,color='r',fmt='x')
            axs[2].errorbar(x=mjds,y=norm_varyBB_Gamma,yerr=norm_err_varyBB_Gamma,color='r',fmt='x')
            axs[2].set_yscale('log')
            axs[3].plot(mjds,lumin_varyBB_Gamma,'rx')
            axs[3].set_yscale('log')
            axs[4].plot(mjds,flux_varyBB_Gamma,'rx')
            axs[4].set_yscale('log')
            axs[5].plot(mjds,cstatarr_varyBB_Gamma,'rx')
            axs[5].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('Gamma')
            axs[1].set_ylabel('Tin (keV)')
            axs[2].set_ylabel('norm')
            axs[3].set_ylabel('Luminosity \n (ergs/s)')
            axs[4].set_ylabel('Flux \n (ergs/cm^2/s)')
            axs[5].set_ylabel('Reduced \n C-stat')

            axs[5].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        def varyBB_PL():
            fig,axs = plt.subplots(7,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep nH (TBNEW) fixed at (' + str(par_dict_varyBB_PL['nH'][0][0]) + ' +- ' + str(par_dict_varyBB_PL['nH'][0][1]) + ') E22 atoms cm^-2' + '\n' + 'Total C-stat: ' + str(total_cstat_varyBB_PL[0]) + '/' + str(total_cstat_varyBB_PL[1]))
            Gamma_varyBB_PL = np.array([par_dict_varyBB_PL['Gamma'][i][0] for i in range(len(par_dict_varyBB_PL['Gamma']))])
            Gamma_err_varyBB_PL = np.array([par_dict_varyBB_PL['Gamma'][i][1] for i in range(len(par_dict_varyBB_PL['Gamma']))])
            FracSctr_varyBB_PL = np.array([par_dict_varyBB_PL['FracSctr'][i][0] for i in range(len(par_dict_varyBB_PL['FracSctr']))])
            FracSctr_err_varyBB_PL = np.array([par_dict_varyBB_PL['FracSctr'][i][1] for i in range(len(par_dict_varyBB_PL['FracSctr']))])
            Tin_varyBB_PL = np.array([par_dict_varyBB_PL['Tin'][i][0] for i in range(len(par_dict_varyBB_PL['Tin']))])
            Tin_err_varyBB_PL = np.array([par_dict_varyBB_PL['Tin'][i][1] for i in range(len(par_dict_varyBB_PL['Tin']))])
            norm_varyBB_PL = np.array([par_dict_varyBB_PL['norm'][i][0] for i in range(len(par_dict_varyBB_PL['norm']))])
            norm_err_varyBB_PL = np.array([par_dict_varyBB_PL['norm'][i][1] for i in range(len(par_dict_varyBB_PL['norm']))])

            cstatarr_varyBB_PL = np.array([cstat_varyBB_PL[i][0]/cstat_varyBB_PL[i][1] for i in range(len(cstat_varyBB_PL)) ])

            axs[0].errorbar(x=mjds,y=Gamma_varyBB_PL,yerr=Gamma_err_varyBB_PL,color='r',fmt='x')
            axs[0].axhline(y=4.5,lw=0.5,alpha=0.5)
            axs[0].set_yscale('log')
            axs[0].set_ylim([0.1,5])
            axs[1].errorbar(x=mjds,y=FracSctr_varyBB_PL,yerr=FracSctr_err_varyBB_PL,color='r',fmt='x')
            axs[1].set_yscale('log')
            axs[1].set_ylim([0.1,1])
            axs[2].errorbar(x=mjds,y=Tin_varyBB_PL,yerr=Tin_err_varyBB_PL,color='r',fmt='x')
            axs[2].set_yscale('log')
            axs[3].errorbar(x=mjds,y=norm_varyBB_PL,yerr=norm_err_varyBB_PL,color='r',fmt='x')
            axs[3].set_yscale('log')
            axs[4].plot(mjds,lumin_varyBB_PL,'rx')
            axs[4].set_yscale('log')
            axs[5].plot(mjds,flux_varyBB_PL,'rx')
            axs[5].set_yscale('log')
            axs[6].plot(mjds,cstatarr_varyBB_PL,'rx')
            axs[6].set_yscale('log')
            axs[6].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('Gamma')
            axs[1].set_ylabel('FracSctr')
            axs[2].set_ylabel('Tin (keV)')
            axs[3].set_ylabel('norm')
            axs[4].set_ylabel('Luminosity \n (ergs/s)')
            axs[5].set_ylabel('Flux \n (ergs/cm^2/s)')
            axs[6].set_ylabel('Reduced \n C-stat')

            axs[6].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        #varyGamma()
        varyBB_Gamma()
        #varyBB_PL()

    def xs_simpl_diskbb():
        par_dict_varyBB_Gamma = get_par(Lv0_dirs.NGC300_2020 + 'spectra_05d/xs_E_tbnew-simpl-diskbb-varyBB-Gamma_par.txt',['nH','Gamma','FracSctr','Tin','norm'])
        cstat_varyBB_Gamma,total_cstat_varyBB_Gamma = get_Cstat(Lv0_dirs.NGC300_2020 + 'spectra_05d/xs_E_tbnew-simpl-diskbb-varyBB-Gamma_cstat.txt')
        lumin_varyBB_Gamma = get_lumin(Lv0_dirs.NGC300_2020 + 'spectra_05d/xs_E_tbnew-simpl-diskbb-varyBB-Gamma_lumin.txt')
        flux_varyBB_Gamma = get_flux(Lv0_dirs.NGC300_2020 + 'spectra_05d/xs_E_tbnew-simpl-diskbb-varyBB-Gamma_flux.txt')

        mjds = np.genfromtxt(Lv0_dirs.NGC300_2020 + 'n300_ulx.xbgsub_cl50_g2020norm_05d.fffphot',usecols=(0),unpack=True)

        def varyBB_Gamma():
            fig,axs = plt.subplots(6,sharex=True,gridspec_kw={'hspace':0})
            fig.suptitle('Keep nH (TBNEW) fixed at (' + str(par_dict_varyBB_Gamma['nH'][0][0]) + ' +- ' + str(par_dict_varyBB_Gamma['nH'][0][1]) + ') E22 atoms cm^-2' + '\n' + 'FracSctr fixed at ' + str(par_dict_varyBB_Gamma['FracSctr'][0][0]) + ' +- ' + str(par_dict_varyBB_Gamma['FracSctr'][0][1]) + '\n' + 'Total C-stat: ' + str(total_cstat_varyBB_Gamma[0]) + '/' + str(total_cstat_varyBB_Gamma[1]))
            Gamma_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['Gamma'][i][0] for i in range(len(par_dict_varyBB_Gamma['Gamma']))])
            Gamma_err_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['Gamma'][i][1] for i in range(len(par_dict_varyBB_Gamma['Gamma']))])
            Tin_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['Tin'][i][0] for i in range(len(par_dict_varyBB_Gamma['Tin']))])
            Tin_err_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['Tin'][i][1] for i in range(len(par_dict_varyBB_Gamma['Tin']))])
            norm_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['norm'][i][0] for i in range(len(par_dict_varyBB_Gamma['norm']))])
            norm_err_varyBB_Gamma = np.array([par_dict_varyBB_Gamma['norm'][i][1] for i in range(len(par_dict_varyBB_Gamma['norm']))])

            cstatarr_varyBB_Gamma = np.array([cstat_varyBB_Gamma[i][0]/cstat_varyBB_Gamma[i][1] for i in range(len(cstat_varyBB_Gamma)) ])

            axs[0].errorbar(x=mjds,y=Gamma_varyBB_Gamma,yerr=Gamma_err_varyBB_Gamma,color='r',fmt='x')
            axs[0].axhline(y=4.5,lw=0.5,alpha=0.5)
            axs[0].axhline(y=1.0,lw=0.5,alpha=0.5)
            axs[1].errorbar(x=mjds,y=Tin_varyBB_Gamma,yerr=Tin_err_varyBB_Gamma,color='r',fmt='x')
            axs[2].errorbar(x=mjds,y=norm_varyBB_Gamma,yerr=norm_err_varyBB_Gamma,color='r',fmt='x')
            axs[2].set_yscale('log')
            axs[3].plot(mjds,lumin_varyBB_Gamma,'rx')
            axs[3].set_yscale('log')
            axs[4].plot(mjds,flux_varyBB_Gamma,'rx')
            axs[4].set_yscale('log')
            axs[5].plot(mjds,cstatarr_varyBB_Gamma,'rx')
            axs[5].axhline(y=1.0,lw=0.5,alpha=0.5)

            axs[0].set_ylabel('Gamma')
            axs[1].set_ylabel('Tin (keV)')
            axs[2].set_ylabel('norm')
            axs[3].set_ylabel('Luminosity \n (ergs/s)')
            axs[4].set_ylabel('Flux \n (ergs/cm^2/s)')
            axs[5].set_ylabel('Reduced \n C-stat')

            axs[5].set_xlabel('Time (MJD)',fontsize=12)

            plt.show()

        varyBB_Gamma()

    #powerlaw()
    #bbodyrad()
    #diskbb()
    #simpl_diskbb()
    #xs_simpl_diskbb()

#data 1:1 bgsub_cl50_100.pha 2:2 bgsub_cl50_101.pha 3:3 bgsub_cl50_102.pha
#response 1 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf 2 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf 3 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf
#arf 1 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf 2 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf 3 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf
#ignore **:0.0-0.285,12.01-**
