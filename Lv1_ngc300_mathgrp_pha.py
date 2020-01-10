#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19th 10:58am 2019

Getting exposure time weighted spectra for NGC300 ULX, by using MATHPHA and GRPPHA

Filename format is:

MJD_05d_bgsub_cl50.pha

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import subprocess
import glob
from tqdm import tqdm
import Lv0_dirs

Lv0_dirs.global_par()

bin_size = '05d' #bins of 1 day!

def mathpha(bin_size,filetype):
    """
    Function that takes in a bin size, and does MATHPHA on the set of pha files.
    The file names are already saved in the binned .ffphot files. The function
    will output pha files of the format 'MJD_binsize_' + filetype + '_cl50.pha'!

    bin_size - bin size in days
    filetype - either 'bgsub' or 'bg' or 'cl'!
    """

    normfile = Lv0_dirs.NGC300 + 'n300_ulx.bgsub_cl50_RGnorm_' + bin_size + '.ffphot'
    mjds = np.genfromtxt(normfile,usecols=(0),unpack=True)
    spectra_files = np.genfromtxt(normfile,dtype='str',usecols=(8),unpack=True)

    for i in range(len(spectra_files)): #so for each BINNED MJD/counts
        exposures = []
        spectra_files_list = spectra_files[i].split(',') #get a list where each entry = one pha file

        if filetype == 'bg':
            filetype_spectra_files_list = [spectra_files_list[j][0:2] + spectra_files_list[j][5:] for j in range(len(spectra_files_list))] #takes away the 'sub' part of the name
        elif filetype == 'cl':
            filetype_spectra_files_list = [spectra_files_list[j][6:] for j in range(len(spectra_files_list))] #takes away the bgsub_ prefix
        elif filetype == 'bgsub':
            filetype_spectra_files_list = spectra_files_list
        for j in range(len(filetype_spectra_files_list)): #for each individual pha file
            fits_file = fits.open(filetype_spectra_files_list[j])
            exptime = fits_file[1].header['EXPOSURE'] #get exposure time
            exposures.append(exptime) #now we have a list of exposure times corresponding to each pha file for a given (binned) MJD/count rate

        total_exptime = sum(exposures) #get total exposure time for THAT MJD/count rate bin

        #### now to build the mathpha command
        mathpha = 'mathpha'
        calculation = ''
        for j in range(len(exposures)):
            if filetype == 'bg' or filetype == 'bgsub':
                if j == 0:
                    calculation += '('
                if j != len(exposures)-1:
                    if int(exposures[j]) == exposures[j]:
                        calculation += str(int(exposures[j])) + '*' + filetype_spectra_files_list[j] + '+'
                    else:
                        calculation += str(exposures[j]) + '*' + filetype_spectra_files_list[j] + '+'
                else:
                    calculation += str(exposures[j]) + '*' + filetype_spectra_files_list[j] + ') / ' + str(total_exptime)

            if filetype == 'cl':
                if j == 0:
                    calculation += '('
                if j != len(exposures)-1:
                    calculation += filetype_spectra_files_list[j] + '+'
                else:
                    calculation += filetype_spectra_files_list[j] + ')'

        if filetype == 'cl':
            unit = 'C'
        if filetype == 'bg' or filetype == 'bgsub':
            unit = 'R'
        outfile = 'outfil=' + str(int(mjds[i])) + '_' + bin_size + '_' + filetype + '_cl50.pha'
        exposure = 'exposure=' + str(total_exptime)
        errmeth = 'errmeth=gaussian'
        properr = 'properr=yes'
        ncomments = 'ncomments=0'
        areascal = 'areascal=NULL'
        clobber = 'clobber=YES'

        logfile = str(int(mjds[i])) + '_' + bin_size + '_' + filetype + '_mathpha.log'
        #print(mathpha+' "'+calculation+ '" '+unit+' '+outfile+' '+exposure+' '+errmeth+' '+properr+' '+ncomments+' '+areascal+' '+clobber)
        with open(logfile,'w') as logtextfile:
            logtextfile.write(subprocess.check_output([mathpha,calculation,unit,outfile,exposure,errmeth,properr,ncomments,areascal,clobber]))
    logtextfile.close()

    return

def grppha(bin_size,filetype):
    """
    Function that takes in a bin size, and does GRPPHA on a set of pha files.
    The input file names will be "$MJD_$binsize_bgsub_cl50.pha". The function
    will output pha files of the format 'grp_$MJD_$binsize_$filetype_cl50.pha'!

    bin_size - bin size in days
    filetype - either 'bgsub' or 'bg' or 'cl'!
    """
    binned_phas = sorted(glob.glob('*'+bin_size+'_'+filetype+'_*pha'))

    command_file = Lv0_dirs.NGC300 + 'grppha_' + bin_size + '_' + filetype + '_commands.go'
    writing = open(command_file,'w')
    #### now to build the grppha command
    grppha = 'grppha'
    chatter = 'chatter=0'
    for i in range(len(binned_phas)):
        backfile = binned_phas[i][:-11] + 'bg_cl50.pha'
        #comm = 'comm="group nicer_channels_to_group.dat & systematics 30-1200 0.02 & chkey ANCRFILE /Volumes/Samsung_T5/nixtiaveonaxis20170601v002.arf & chkey RESPFILE /Volumes/Samsung_T5/nicer_upd_d52.rmf & chkey BACKFILE ' + backfile + ' & exit"'
        comm = 'comm="group nicer_channels_to_group.dat & systematics 30-1200 0.02 & chkey ANCRFILE /Volumes/Samsung_T5/nixtiaveonaxis20170601v002.arf & chkey RESPFILE /Volumes/Samsung_T5/nicer_upd_d52.rmf & exit"'
        infile = 'infile="' + binned_phas[i] + '"'
        outfile = 'outfile="grp_' + binned_phas[i] + '"'
        writing.write(grppha+' '+infile+' '+outfile+' '+chatter+' '+comm+'\n')
    writing.close()

    #subprocess.check_output(['source',command_file,'>','grppha_'+bin_size+'.log'],shell=True)

if __name__ == "__main__":
    #mathpha(bin_size,'bgsub')
    grppha(bin_size,'bgsub')
    #bg_mathpha(bin_size)
    #bg_grppha(bin_size)

################################## DEPRECATED ##################################

"""
def js_groupspec(bin_size):

    Function that takes in a bin size, and runs js_groupspec.pro on a set of pha files.
    This is to: 1) get the RIGHT energy resolution, to get appropriately independent bins.
    The current data oversamples!

    binned_phas = sorted(glob.glob('*'+bin_size+'*pha'))
    js_groupspec = 'js_groupspec,'
    for i in range(len(binned_phas)):
        infile = 'infile="' + binned_phas[i] + '",'
        outfile = 'outfile="rb_' + binned_phas[i] + '",'
"""

"""
def bg_mathpha(bin_size):
    Function that takes in a bin size, and does MATHPHA on the set of BACKGROUND pha files.
    The file names are already saved in the binned .ffphot files. The function
    will output pha files of the format 'MJD_binsize_bg_cl50.pha'!


    normfile = Lv0_dirs.NGC300 + 'n300_ulx.bgsub_cl50_RGnorm_' + bin_size + '.ffphot'
    mjds = np.genfromtxt(normfile,usecols=(0),unpack=True)
    spectra_files = np.genfromtxt(normfile,dtype='str',usecols=(8),unpack=True)

    for i in range(len(spectra_files)): #so for each BINNED MJD/counts
        exposures = []
        spectra_files_list = spectra_files[i].split(',') #get a list where each entry = one pha file

        bgspectra_files_list = ['bg' + spectra_files_list[k][5:] for k in range(len(spectra_files_list))]
        for j in range(len(bgspectra_files_list)): #for each individual pha file
            fits_file = fits.open(bgspectra_files_list[j])
            exptime = fits_file[1].header['EXPOSURE'] #get exposure time
            exposures.append(exptime) #now we have a list of exposure times corresponding to each pha file for a given (binned) MJD/count rate

        total_exptime = sum(exposures) #get total exposure time for THAT MJD/count rate bin

        #### now to build the mathpha command
        mathpha = 'mathpha'
        calculation = ''
        for j in range(len(exposures)):
            if j == 0:
                calculation += '('
            if j != len(exposures)-1:
                if int(exposures[j]) == exposures[j]:
                    calculation += str(int(exposures[j])) + '*' + bgspectra_files_list[j] + '+'
                else:
                    calculation += str(exposures[j]) + '*' + bgspectra_files_list[j] + '+'
            else:
                calculation += str(exposures[j]) + '*' + bgspectra_files_list[j] + ')'
        calculation += ' / ' + str(total_exptime)
        unit = 'C' #for rate ; use 'C' for counts
        outfile = 'outfil=' + str(int(mjds[i])) + '_' + bin_size + '_bg_cl50.pha'
        exposure = 'exposure=' + str(total_exptime)
        errmeth = 'errmeth=gaussian'
        properr = 'properr=yes'
        ncomments = 'ncomments=0'
        areascal = 'areascal=NULL'
        clobber = 'clobber=YES'
        #print(mathpha,calculation,unit,outfile,exposure,errmeth,properr,ncomments,areascal,clobber)
        logfile = str(int(mjds[i])) + '_' + bin_size + '_mathpha.log'
        print(mathpha+' "'+calculation+ '" '+unit+' '+outfile+' '+exposure+' '+errmeth+' '+properr+' '+ncomments+' '+areascal+' '+clobber)
        with open(logfile,'w') as logtextfile:
            logtextfile.write(subprocess.check_output([mathpha,calculation,unit,outfile,exposure,errmeth,properr,ncomments,areascal,clobber]))
    logtextfile.close()

    return

def bg_grppha(bin_size):

    Function that takes in a bin size, and does GRPPHA on a set of pha files.
    The input file names will be "MJD_binsize_bgsub_cl50.pha". The function
    will output pha files of the format 'grp_MJD_binsize_bgsub_cl50.pha'!

    binned_phas = sorted(glob.glob('*'+bin_size+'_bg_*pha'))

    command_file = 'grppha_bg_' + bin_size + '_commands.go'
    writing = open(command_file,'w')
    #### now to build the grppha command
    grppha = 'grppha'
    chatter = 'chatter=0'
    comm = 'comm="group nicer_channels_to_group.dat & systematics 30-1200 0.02 & chkey ANCRFILE /Volumes/Samsung_T5/nixtiaveonaxis20170601v002.arf & chkey RESPFILE /Volumes/Samsung_T5/nicer_upd_d52.rmf & exit"'
    for i in range(len(binned_phas)):
        infile = 'infile="' + binned_phas[i] + '"'
        outfile = 'outfile="grp_' + binned_phas[i] + '"'
        writing.write(grppha+' '+infile+' '+outfile+' '+chatter+' '+comm+'\n')
    writing.close()

    return
"""
