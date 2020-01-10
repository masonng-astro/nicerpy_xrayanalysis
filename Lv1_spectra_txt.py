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
binspectra_bgsub = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/jsgrp_*_cl_*.pha'))
spectra_all = '/Volumes/Samsung_T5/NGC300_ULX/jsgrp_cl_spectra.txt'
spectra = np.array(open(spectra_all,'r').read().split('\n'))

separator = list(np.where(spectra=='NO NO NO NO')[0])
separator.insert(0,2)
separator.insert(len(separator),len(spectra)-1)

for i in tqdm(range(len(binspectra_bgsub))):
    spectrum_file = binspectra_bgsub[i][:-3] + 'txt'
    spectrum = open(spectrum_file,'w')

    for j in range(separator[i]+1,separator[i+1]):
        spectrum.write(spectra[j] + '\n')

    spectrum.close()
"""

#not_to_use = ['58279','58334','58399','58464','58514','58599','58604']
not_to_use = []
#total_E = []
##### time to plot all 41 spectra and save them onto PDFs...!
binspectra_cl = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/jsgrp_58*_cl_*.txt'))
binspectra_bg = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/grp_58*_bg_*.txt'))
binspectra_bgsub = sorted(glob.glob('/Volumes/Samsung_T5/NGC300_ULX/grp_58*_bgsub_*.txt'))

"""
with PdfPages('jsgrp_binned_spectra.pdf') as pdf:
    for i in tqdm(range(len(binspectra_cl))):
        if binspectra_cl[i][-21:-16] not in not_to_use:
            E_cl,E_cl_unc,flux_cl,flux_cl_unc = np.genfromtxt(binspectra_cl[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)
            E_bg,E_bg_unc,flux_bg,flux_bg_unc = np.genfromtxt(binspectra_bg[i],usecols=(0,1,2,3),unpack=True,skip_footer=1)
            #total_E += list(E)
            plt.errorbar(E_cl,flux_cl,xerr=E_cl_unc,yerr=flux_cl_unc,fmt='r+')
            plt.errorbar(E_bg,flux_bg,xerr=E_bg_unc,yerr=flux_bg_unc,fmt='b+')
            plt.title('File name: ' + binspectra_cl[i] + '\n' + 'Background: ' + binspectra_bg[i],fontsize=12)
            plt.xlabel('Energy (keV)',fontsize=10)
            plt.ylabel('Normalized counts /s/keV',fontsize=10)
            plt.legend(('cl','bg'),loc='best')
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


################ 9/30 - to implement the 'area under the curve' idea...
#mjds = ['58239','58244','58249','58254','58259','58264','58269','58274','58279',
#        '58284','58289','58294','58309','58314','58324','58329','58334','58339',
#        '58344','58349','58384','58389','58394','58399','58409','58449','58454',
#        '58459','58464','58469','58474','58479','58484','58489','58494','58499',
#        '58504','58509','58514','58599','58604']

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

#data 1:1 bgsub_cl50_100.pha 2:2 bgsub_cl50_101.pha 3:3 bgsub_cl50_102.pha
#response 1 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf 2 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf 3 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.rmf
#arf 1 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf 2 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf 3 /Users/masonng/Documents/MIT/Research/nicer-data/nicer.arf
#ignore **:0.0-0.285,12.01-**
