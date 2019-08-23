#!/usr/bin/env python2
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

def mathpha(bin_size):
    """
    Function that takes in a bin size, and does MATHPHA on the set of pha files.
    The file names are already saved in the binned .ffphot files. The function
    will output pha files of the format 'MJD_binsize_bgsub_cl50.pha'!
    """

    normfile = Lv0_dirs.NGC300 + 'n300_ulx.bgsub_cl50_RGnorm_' + bin_size + '.ffphot'
    mjds = np.genfromtxt(normfile,usecols=(0),unpack=True)
    spectra_files = np.genfromtxt(normfile,dtype='str',usecols=(8),unpack=True)

    for i in range(len(spectra_files)): #so for each BINNED MJD/counts
        exposures = []
        spectra_files_list = spectra_files[i].split(',') #get a list where each entry = one pha file
        for j in range(len(spectra_files_list)): #for each individual pha file
            fits_file = fits.open(spectra_files_list[j])
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
                    calculation += str(int(exposures[j])) + '*' + spectra_files_list[j] + '+'
                else:
                    calculation += str(exposures[j]) + '*' + spectra_files_list[j] + '+'
            else:
                calculation += str(exposures[j]) + '*' + spectra_files_list[j] + ')'
        calculation += ' / ' + str(total_exptime)
        unit = 'R' #for rate ; use 'C' for counts
        outfile = 'outfil=' + str(int(mjds[i])) + '_' + bin_size + '_bgsub_cl50.pha'
        exposure = 'exposure=' + str(total_exptime)
        errmeth = 'errmeth=gaussian'
        properr = 'properr=yes'
        ncomments = 'ncomments=0'
        areascal = 'areascal=NULL'
        clobber = 'clobber=YES'

        logfile = str(int(mjds[i])) + '_' + bin_size + '_mathpha.log'
        #print(mathpha+' "'+calculation+ '" '+unit+' '+outfile+' '+exposure+' '+errmeth+' '+properr+' '+ncomments+' '+areascal+' '+clobber)
        with open(logfile,'w') as logtextfile:
            logtextfile.write(subprocess.check_output([mathpha,calculation,unit,outfile,exposure,errmeth,properr,ncomments,areascal,clobber]))
    logtextfile.close()

    return

def grppha(bin_size):
    """
    Function that takes in a bin size, and does GRPPHA on a set of pha files.
    The input file names will be "MJD_binsize_bgsub_cl50.pha". The function
    will output pha files of the format 'grp_MJD_binsize_bgsub_cl50.pha'!
    """
    binned_phas = sorted(glob.glob('*'+bin_size+'*pha'))

    command_file = 'grppha_' + bin_size + '_commands.go'
    writing = open(command_file,'w')
    #### now to build the grppha command
    grppha = 'grppha'
    chatter = 'chatter=0'
    comm = 'comm="group 30 1200 3 & systematics 30-1200 0.02 & chkey ANCRFILE /Volumes/Samsung_T5/nixtiaveonaxis20170601v002.arf & chkey RESPFILE /Volumes/Samsung_T5/nicer_upd_d52.rmf & exit"'
    for i in range(len(binned_phas)):
        infile = 'infile="' + binned_phas[i] + '"'
        outfile = 'outfile="grp_' + binned_phas[i] + '"'
        writing.write(grppha+' '+infile+' '+outfile+' '+chatter+' '+comm+'\n')
    writing.close()

    #subprocess.check_output(['source',command_file,'>','grppha_'+bin_size+'.log'],shell=True)

    return

if __name__ == "__main__":
    #mathpha(bin_size)
    grppha(bin_size)
