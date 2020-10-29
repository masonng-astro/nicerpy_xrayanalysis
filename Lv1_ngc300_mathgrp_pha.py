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
import pathlib
import subprocess
import glob
from tqdm import tqdm
import Lv0_dirs

Lv0_dirs.global_par()

bin_size = '03d' #bins of n day!

def orb_phase_rate(orb_phase):
    """
    Returns the corresponding rate from the 20-binned folded profile
    """

    top = 0.022177
    bottom = 0.005351

    if (orb_phase > 0.25 and orb_phase <= 0.80): #off-eclipse
        rate = top

    elif (orb_phase > 0.9 and orb_phase <= 1) or (orb_phase >= 0.0 and orb_phase <= 0.1): #on-eclipse
        rate = bottom

    elif (orb_phase > 0.1 and orb_phase <= 0.15):
        ##### the rate is 0.011361
        rate = 0.011361

    elif (orb_phase > 0.15 and orb_phase <= 0.20):
        ##### the rate is 0.014107
        rate = 0.014107

    elif (orb_phase > 0.20 and orb_phase <= 0.25):
        ##### the rate is 0.019135
        rate = 0.019135

    elif (orb_phase > 0.80 and orb_phase <= 0.85):
        ##### the rate is 0.015688
        rate = 0.015688

    elif (orb_phase > 0.85 and orb_phase <= 0.90):
        ##### the rate is 0.011368
        rate = 0.011368

    return rate

def combine_back_scal(init_back_array,fake_spec,T0,Porb):
    """
    Given a timing model/ephemeris, figure out which orbital phase the events in a
    background spectrum/event file are in (defined by some centroid time), then
    write a mathpha command combining those files appropriately in rate space.
    Unlike combine_back, this SCALES the spectra/GTIs in the ingress/egress sections
    based on the folded PN light curve!

    Remember that:
    025-080 -> off
    080-090 -> ingress
    090-010 -> on
    010-025 -> egress

    init_back_array - array of input background file
    fake_spec - array of files split by orbital phase
    T0 - reference time ; deepest point of eclipse
    Porb - orbital period in days
    """
    command_file = open(Lv0_dirs.NGC300_2020 + '3C50_X1scale_mathpha.go','w')
    on_e = fake_spec[0]
    off_e = fake_spec[1]

    top = 0.022177
    bottom = 0.005351

    counter = 0
    for i in tqdm(range(len(init_back_array))):
        init_back = init_back_array[i]
        gti_no = init_back[-8:-4]
        init_cl = Lv0_dirs.NGC300_2020 + 'cl50/pha/cl50_' + gti_no + '.pha'

        backfile = str(pathlib.Path(init_back).name)
        comb_spec_folder = Lv0_dirs.NGC300_2020 + 'bg_cl50/xspha/' #x for 'extra', s for 'scale'

        tstart = fits.open(init_cl)[1].header['TSTART']
        tstop = fits.open(init_cl)[1].header['TSTOP']
        centroid_t = fits.open(init_cl)[1].header['MJDREFI'] + fits.open(init_cl)[1].header['MJDREFF'] + ((tstart+tstop)/2)/86400

        #first assuming the timing model, determine how much time (in d) has passed
        orb_phase = (centroid_t-T0)%Porb/Porb

        calculation = ''
        rate = orb_phase_rate(orb_phase)

        on_factor = round((top-rate)/(top-bottom),6)
        off_factor = round((rate-bottom)/(top-bottom),6)
        calculation += "('" + init_back + "'+" + str(off_factor) + "*'" + off_e + "'+" + str(on_factor) + "*'" + on_e + "')"

        bg_exp = fits.open(init_back)[1].header['EXPOSURE']
        eclipse_exp = 1e6 #1 Ms was used in the faked spectrum
        #### now to build the mathpha command
        mathpha = 'mathpha'

        unit = 'R' #for rate ; use 'C' for counts
        outfile = "outfil='" + comb_spec_folder + "xs" + backfile + "'"
        exposure = 'exposure=' + str(bg_exp + eclipse_exp)
        errmeth = 'errmeth=gaussian'
        properr = 'properr=yes'
        ncomments = 'ncomments=0'
        areascal = 'areascal=NULL'
        clobber = 'clobber=YES'

        mathpha_line = mathpha+' "'+calculation+ '" '+unit+' '+outfile+' '+exposure+' '+errmeth+' '+properr+' '+ncomments+' '+areascal+' '+clobber
        command_file.write(mathpha_line + '\n')

    command_file.close()

def combine_back(init_back_array,fake_spec,T0,Porb):
    """
    Given a timing model/ephemeris, figure out which orbital phase the events in a
    background spectrum/event file are in (defined by some centroid time), then
    write a mathpha command combining those files appropriately in rate space

    init_back_array - array of input background file
    fake_spec - array of files split by orbital phase
    T0 - reference time ; deepest point of eclipse
    Porb - orbital period in days
    """
    command_file = open(Lv0_dirs.NGC300_2020 + '3C50_X1_mathpha.go','w')

    for i in tqdm(range(len(init_back_array))):
        init_back = init_back_array[i]
        gti_no = init_back[-8:-4]
        init_cl = Lv0_dirs.NGC300_2020 + 'cl50/pha/cl50_' + gti_no + '.pha'

        backfile = str(pathlib.Path(init_back).name)
        comb_spec_folder = Lv0_dirs.NGC300_2020 + 'bg_cl50/xpha/'

        tstart = fits.open(init_cl)[1].header['TSTART']
        tstop = fits.open(init_cl)[1].header['TSTOP']
        centroid_t = fits.open(init_cl)[1].header['MJDREFI'] + fits.open(init_cl)[1].header['MJDREFF'] + ((tstart+tstop)/2)/86400

        #first assuming the timing model, determine how much time (in d) has passed
        orb_phase = (centroid_t-T0)%Porb/Porb

        if (orb_phase >= 0.0 and orb_phase <= 0.20) or (orb_phase >= 0.85 and orb_phase <= 1):
            eclipse = fake_spec[0]
        elif (orb_phase > 0.2 and orb_phase < 0.85):
            eclipse = fake_spec[1]

        bg_exp = fits.open(init_back)[1].header['EXPOSURE']
        eclipse_exp = fits.open(eclipse)[1].header['EXPOSURE']
        #### now to build the mathpha command
        mathpha = 'mathpha'
        calculation = ''
        calculation += "('" + init_back + "'+'" + eclipse + "')"

        unit = 'R' #for rate ; use 'C' for counts
        outfile = "outfil='" + comb_spec_folder + "x" + backfile + "'"
        exposure = 'exposure=' + str(bg_exp + eclipse_exp)
        errmeth = 'errmeth=gaussian'
        properr = 'properr=yes'
        ncomments = 'ncomments=0'
        areascal = 'areascal=NULL'
        clobber = 'clobber=YES'

        mathpha_line = mathpha+' "'+calculation+ '" '+unit+' '+outfile+' '+exposure+' '+errmeth+' '+properr+' '+ncomments+' '+areascal+' '+clobber
        command_file.write(mathpha_line + '\n')

    command_file.close()

    return

def mathpha(bin_size,filetype):
    """
    Function that takes in a bin size, and does MATHPHA on the set of pha files.
    The file names are already saved in the binned .ffphot files. The function
    will output pha files of the format 'MJD_binsize_' + filetype + '_cl50.pha'!

    bin_size - bin size in days
    filetype - either 'bgsub' or 'bg' or 'cl' or 'x' or 'xbgsub'
    x means "extra background"
    """

    normfile = Lv0_dirs.NGC300_2020 + 'n300_ulx.bgsub_cl50_g2020norm_' + bin_size + '.fffphot'
    mjds = np.genfromtxt(normfile,usecols=(0),unpack=True)
    spectra_files = np.genfromtxt(normfile,dtype='str',usecols=(9),unpack=True)

    for i in range(len(spectra_files)): #so for each BINNED MJD/counts
        exposures = []
        truncated_name = []
        spectra_files_list = spectra_files[i].split(',') #get a list where each entry = one pha file
        cl50_files = [spectra_files_list[j][:34] + spectra_files_list[j][40:49] + spectra_files_list[j][55:] for j in range(len(spectra_files_list))] #takes away the bgsub_ prefix

        if filetype == 'bg':
            filetype_spectra_files_list = [spectra_files_list[j][:36] + spectra_files_list[j][39:51] + spectra_files_list[j][54:] for j in range(len(spectra_files_list))] #takes away the 'sub' part of the name
        elif filetype == 'cl':
            filetype_spectra_files_list = [spectra_files_list[j][:34] + spectra_files_list[j][40:49] + spectra_files_list[j][55:] for j in range(len(spectra_files_list))] #takes away the bgsub_ prefix
        elif filetype == 'bgsub':
            filetype_spectra_files_list = spectra_files_list
        elif filetype == 'x':
            filetype_spectra_files_list = [spectra_files_list[j][:36] + spectra_files_list[j][39:45] + 'xpha/xbg' + spectra_files_list[j][54:] for j in range(len(spectra_files_list))]
        elif filetype == 'xbgsub':
            filetype_spectra_files_list = [spectra_files_list[j][:45] + 'xpha/xbgsub_' + spectra_files_list[j][55:] for j in range(len(spectra_files_list))]
        elif filetype == 'xsbg':
            filetype_spectra_files_list = [spectra_files_list[j][:36] + spectra_files_list[j][39:45] + 'xspha/xsbg' + spectra_files_list[j][54:] for j in range(len(spectra_files_list))]
        elif filetype == 'xsbgsub':
            filetype_spectra_files_list = [spectra_files_list[j][:45] + 'xspha/xsbgsub_' + spectra_files_list[j][55:] for j in range(len(spectra_files_list))]

        for j in range(len(filetype_spectra_files_list)): #for each individual pha file
            fits_file = fits.open(filetype_spectra_files_list[j])
            exptime = fits_file[1].header['EXPOSURE'] #get exposure time
            if filetype=='bg' or filetype=='x' or filetype=='xbgsub' or filetype=='xsbg' or filetype=='xsbgsub':
                exptime = fits.open(cl50_files[j])[1].header['EXPOSURE']
            exposures.append(exptime) #now we have a list of exposure times corresponding to each pha file for a given (binned) MJD/count rate

        total_exptime = sum(exposures) #get total exposure time for THAT MJD/count rate bin

        for j in range(len(filetype_spectra_files_list)):
            long_name = filetype_spectra_files_list[j]
            short_name = str(pathlib.Path(long_name).name)
            subprocess.run(['cp',long_name,short_name])
            truncated_name.append(short_name)

        #### now to build the mathpha command
        ################## NOTE, FOR BG/BGSUB FILES, MAKE SURE YOU USE EXPOSURE
        ################## TIMES FROM THE CL50 FILES... SINCE THAT'S HOW MUCH THE
        ################## OBSERVATION WAS EXPOSED FOR...
        mathpha = 'mathpha'
        calculation = ''
        for j in range(len(exposures)):
            if filetype == 'bg' or filetype == 'bgsub' or filetype == 'x' or filetype == 'xbgsub' or filetype == 'xsbg' or filetype == 'xsbgsub':
                if j == 0:
                    calculation += '('
                if j != len(exposures)-1:
                    if int(exposures[j]) == exposures[j]:
                        calculation += str(int(exposures[j])) + '*' + truncated_name[j] + '+'
                    else:
                        calculation += str(exposures[j]) + '*' + truncated_name[j] + '+'
                else:
                    calculation += str(exposures[j]) + '*' + truncated_name[j] + ') / ' + str(total_exptime)

            if filetype == 'cl':
                if j == 0:
                    calculation += '('
                if j != len(exposures)-1:
                    calculation += truncated_name[j] + '+'
                else:
                    calculation += truncated_name[j] + ')'

        if mjds[i] == 58449:
            print(calculation)

        if filetype == 'cl':
            unit = 'C'
        if filetype == 'bg' or filetype == 'bgsub' or filetype == 'x' or filetype == 'xbgsub' or filetype == 'xsbg' or filetype == 'xsbgsub':
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
            output = subprocess.run([mathpha,calculation,unit,outfile,exposure,errmeth,properr,ncomments,areascal,clobber],capture_output=True,text=True)
            logtextfile.write(output.stdout)
            logtextfile.write('*------------------------------* \n')
            logtextfile.write(output.stderr)
            logtextfile.close()

    ### removing the temporary spectra
    for i in range(len(truncated_name)):
        subprocess.run(['rm','-r',truncated_name[i]])

    return

def grppha(bin_size,filetype):
    """
    Function that takes in a bin size, and does GRPPHA on a set of pha files.
    The input file names will be "$MJD_$binsize_bgsub_cl50.pha". The function
    will output pha files of the format 'grp_$MJD_$binsize_$filetype_cl50.pha'!

    bin_size - bin size in days
    filetype - either 'bgsub' or 'bg' or 'cl' or 'x' or 'xbgsub' or 'xsbg' or 'xsbgsub'!
    """
    binned_phas = sorted(glob.glob('*'+bin_size+'_'+filetype+'_*pha'))

    command_file = Lv0_dirs.NGC300_2020 + 'grppha_' + bin_size + '_' + filetype + '_commands.go'
    writing = open(command_file,'w')
    #### now to build the grppha command
    grppha = 'grppha'
    chatter = 'chatter=0'
    for i in range(len(binned_phas)):
        backfile = binned_phas[i][:-11] + 'bg_cl50.pha'
        #comm = 'comm="group nicer_channels_to_group.dat & systematics 30-1200 0.02 & chkey ANCRFILE /Volumes/Samsung_T5/nicer-consim135p-teamonly-array50.arf & chkey RESPFILE /Volumes/Samsung_T5/nicer-rmf6s-teamonly-array50.rmf & chkey BACKFILE ' + backfile + ' & exit"'
        comm = 'comm="group nicer_channels_to_group.dat & systematics 30-1200 0.02 & chkey ANCRFILE /Volumes/Samsung_T5/nicer-consim135p-teamonly-array50.arf & chkey RESPFILE /Volumes/Samsung_T5/nicer-rmf6s-teamonly-array50.rmf & exit"'
        infile = 'infile="' + binned_phas[i] + '"'
        outfile = 'outfile="grp_' + binned_phas[i] + '"'
        writing.write(grppha+' '+infile+' '+outfile+' '+chatter+' '+comm+'\n')
    writing.close()

    #subprocess.run(['source',command_file,'>','grppha_'+bin_size+'.log'],shell=True)

if __name__ == "__main__":

    init_back_array = sorted(glob.glob(Lv0_dirs.NGC300_2020 + 'bg_cl50/pha/*.pha'))
    on_eclipse = Lv0_dirs.NGC300_XMM + 'ngc300x1_oneclipse.fak'
    off_eclipse = Lv0_dirs.NGC300_XMM + 'ngc300x1_offeclipse.fak'
    fake_spec = np.array([on_eclipse,off_eclipse])
    T0 = 58239.3498
    Porb = (1/8.4712e-6)/86400

    #combine_back(init_back_array,fake_spec,T0,Porb)
    #combine_back_scal(init_back_array,fake_spec,T0,Porb)

    #mathpha(bin_size,'bg')
    #mathpha(bin_size,'cl')
    #mathpha(bin_size,'bgsub')

    #mathpha('03d','xbgsub')
    #mathpha('05d','xbgsub')
    #mathpha('10d','xbgsub')

    #mathpha('05d','xsbgsub')

    #mathpha('05d','xsbg')

    #grppha('03d','bg')
    #grppha('05d','bg')
    #grppha('10d','bg')

    #grppha('03d','cl')
    #grppha('05d','cl')
    #grppha('10d','cl')

    #grppha('03d','bgsub')
    #grppha('05d','bgsub')
    #grppha('10d','bgsub')

    #grppha('03d','x')
    #grppha('05d','x')
    #grppha('10d','x')

    #grppha('05d','xsbg')

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
            logtextfile.write(subprocess.run([mathpha,calculation,unit,outfile,exposure,errmeth,properr,ncomments,areascal,clobber]))
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
