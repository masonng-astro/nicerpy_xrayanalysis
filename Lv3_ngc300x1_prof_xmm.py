#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 7 11:14am 2020

Analyzing the folded orbital profiles of NGC 300 X-1 with XMM-Newton data

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import subprocess
import pathlib
import math
#from stingray.pulse.pulsar import pulse_phase

##### The data

eventfile_pn = '/Volumes/Samsung_T5/NGC300_XMMdata/ngc300x1_pn.evt'
offset_pn = 0.09249852833769791
times_pn = fits.open(eventfile_pn)[1].data['TIME']
#phases_pn = pulse_phase(times_pn-times_pn[0],8.4712e-6,ph0=offset_pn)
gtis_pn = fits.open(eventfile_pn)[59].data #59 for pn, 15 for mos1, 19 for mos2
T_pn = sum([ gtis_pn[i]['STOP']-gtis_pn[i]['START'] for i in range(len(gtis_pn)) ]) #exposure time

eventfile_mos1 = '/Volumes/Samsung_T5/NGC300_XMMdata/ngc300x1_mos1.evt'
offset_mos1 = 0.08930142107215873
times_mos1 = fits.open(eventfile_mos1)[1].data['TIME']
#phases_mos1 = pulse_phase(times_mos1-times_mos1[0],8.4712e-6,ph0=offset_mos1)
gtis_mos1 = fits.open(eventfile_mos1)[15].data #59 for pn, 15 for mos1, 19 for mos2
T_mos1 = sum([ gtis_mos1[i]['STOP']-gtis_mos1[i]['START'] for i in range(len(gtis_mos1)) ]) #exposure time

eventfile_mos2 = '/Volumes/Samsung_T5/NGC300_XMMdata/ngc300x1_mos2.evt'
offset_mos2 = 0.06837599396512264
times_mos2 = fits.open(eventfile_mos2)[1].data['TIME']
#phases_mos2 = pulse_phase(times_mos2-times_mos2[0],8.4712e-6,ph0=offset_mos2)
gtis_mos2 = fits.open(eventfile_mos2)[19].data #59 for pn, 15 for mos1, 19 for mos2
T_mos2 = sum([ gtis_mos2[i]['STOP']-gtis_mos2[i]['START'] for i in range(len(gtis_mos2)) ]) #exposure time

MJDREF = fits.open(eventfile_pn)[1].header['MJDREF']

phases = np.arange(0.025,2,0.05)
##### FOR PN:
pn_prof = np.array([1604.0,2164.0,1937.4058099752326,2676.0,3300.0,3839.637730914639,4038.0,4056.0,3540.7859071725998,2140.5, 3804.0,3668.75965331829,3826.5,3006.0,
1034.213403644357,1404.0,483.0,703.8586733495172,474.2092322833178,487.5]*2)/T_pn
pn_inter = interp1d(phases,pn_prof,kind='cubic')
pn_prof_err = np.array([40.049968789001575,46.51881339845203,55.65052783180561,89.5991071384085,99.49874371066204,107.4471256562774,110.06361796706489,110.30865786510145,100.06852620303427,
56.66348030257233,75.53806987208498,74.28709253420669,75.76113779504647,67.14908785679812,39.400392759523484,45.89117562233503,26.916537667389523,32.50117373800949,
26.67635348746897,27.04163456597996]*2)/T_pn

##### FOR MOS1:
mos1_prof = np.array([448.0,572.0,645.7875788096308,864.0,1152.0,1032.3053915630078,1167.0,1236.0,978.8970623790145,676.5,1014.0,1073.182141931018,1099.5,943.5,
444.50831287810894,507.0,318.0,210.0,208.69091755823553,178.58161485859577]*2)/T_mos1
mos1_inter = interp1d(phases,mos1_prof,kind='cubic')
mos1_prof_err = np.array([21.166010488516726,23.916521486202797,29.055356803825024,50.911688245431435,58.78775382679629,55.73925190095661,59.169248769948084,60.89334939055334,
48.883786528584615,31.855140872392965,38.99999999999997,40.190987067477344,40.610959112042714,37.61980861195333,25.8365208870395,27.57716446627533,21.840329667841537,
17.748239349298878, 17.700932599765082, 16.37054979377456]*2)/T_mos1

##### FOR MOS2:
mos2_prof = np.array([336.0,651.0,588.6128641386009,855.3569611778419,1308.0,1011.0,1218.0,1116.0,1323.0,804.5623109294818,1038.0,1129.4531069776747,1100.2580003973176,1054.5,
759.8367288864442,486.4282601354391,439.5,171.0,304.56707598060035,150.0]*2)/T_mos2
mos2_inter = interp1d(phases,mos2_prof,kind='cubic')
mos2_prof_err = np.array([18.33030277982336,25.514701644346147,24.9179694998494,44.77142958345511,62.64183905346333,55.07267925205741,60.44832503882968,57.86190456595775,63.0,36.57101413315826,
39.458839313897684,41.21432922222213,40.63895649552821,39.77122075068852,33.77886277372683,27.023792229746615,25.675864152935514,16.015617378046993,21.376418084077713,15.0]*2)/T_mos2

def step_n_ramp(phase,prof,prof_err,start_phase,end_phase,pguess):
    """
    Fitting 6-model step-and-ramp parameters to the folded profile
    """
    phase_model = np.linspace(start_phase,end_phase,101)
    x = phase[(phase>=start_phase)&(phase<=end_phase)]
    y = prof[(phase>=start_phase)&(phase<=end_phase)]
    y_err = prof_err[(phase>=start_phase)&(phase<=end_phase)]

    def piecewise_linear(x,b1,b2,b3,b4,top,bottom):
        return np.piecewise(x, [(x>=start_phase)&(x<=b1), (x>b1)&(x<=b2), (x>b2)&(x<=b3), (x>b3)&(x<=b4), (x>b4)&(x<=end_phase)], [lambda x:top, lambda x:((bottom-top)/(b2-b1)*x+bottom-(bottom-top)/(b2-b1)*b2), lambda x:bottom, lambda x:((top-bottom)/(b4-b3)*x+top-(top-bottom)/(b4-b3)*b4), lambda x:top])

    plt.figure()
    popt,pcov = curve_fit(piecewise_linear,x,y,p0=pguess,sigma=y_err,absolute_sigma=True)
    pars = popt
    pars_err = np.diag(np.sqrt(pcov))
    print('Top: ' + str(pars[4]) + ' +- ' + str(pars_err[4]))
    print('Bottom: ' + str(pars[5]) + ' +- ' + str(pars_err[5]))
    print('Vertex 1: ' + str(pars[0]) + ' +- ' + str(pars_err[0]))
    print('Vertex 2: ' + str(pars[1]) + ' +- ' + str(pars_err[1]))
    print('Vertex 3: ' + str(pars[2]) + ' +- ' + str(pars_err[2]))
    print('Vertex 4: ' + str(pars[3]) + ' +- ' + str(pars_err[3]))
    plt.plot(phase_model,piecewise_linear(phase_model,*popt),'k-')

    ##### Plotting the folded profiles themselves
    plt.errorbar(x=phase,y=prof,yerr=prof_err,color='r',drawstyle='steps-mid')
    plt.title('Exposure-corrected (profiles from Stingray)',fontsize=12)
    plt.xlabel('Phase',fontsize=12)
    plt.ylabel('Counts/s',fontsize=12)
    plt.legend(('Piecewise fit','Exposure-corrected profile'),loc='best',fontsize=12)

def plot_profs(interpolate):
    spacing = 0.05/20
    phases_long = np.arange(0.025,1.975+spacing,spacing)

    plt.figure()
    if interpolate == False:
        plt.errorbar(x=phases,y=pn_prof,yerr=pn_prof_err,color='r',drawstyle='steps-mid')
        plt.errorbar(x=phases,y=mos1_prof,yerr=mos1_prof_err,color='b',drawstyle='steps-mid')
        plt.errorbar(x=phases,y=mos2_prof,yerr=mos2_prof_err,color='k',drawstyle='steps-mid')
        plt.legend(('PN','MOS1','MOS2'),loc='best',fontsize=12)
    if interpolate == True:
        plt.plot(phases_long,pn_inter(phases_long),'r--')
        plt.plot(phases_long,mos1_inter(phases_long),'b--')
        plt.plot(phases_long,mos2_inter(phases_long),'k--')
        plt.legend(('PN interp1d','MOS1 interp1d','MOS2 interp1d'),loc='best',fontsize=12)
    elif interpolate == 'both':
        plt.errorbar(x=phases,y=pn_prof,yerr=pn_prof_err,color='r',drawstyle='steps-mid')
        plt.plot(phases_long,pn_inter(phases_long),'r--')
        plt.errorbar(x=phases,y=mos1_prof,yerr=mos1_prof_err,color='b',drawstyle='steps-mid')
        plt.plot(phases_long,mos1_inter(phases_long),'b--')
        plt.errorbar(x=phases,y=mos2_prof,yerr=mos2_prof_err,color='k',drawstyle='steps-mid')
        plt.plot(phases_long,mos2_inter(phases_long),'k--')
        plt.legend(('PN','PN interp1d','MOS1','MOS1 interp1d','MOS2','MOS2 interp1d'),loc='best',fontsize=12)

    plt.title('Exposure-corrected (using Stingray fold_events)',fontsize=12)
    plt.xlabel('Phase',fontsize=12)
    plt.ylabel('Counts/s',fontsize=12)
    plt.show()

def shift_profs(nbins,prof,prof_err,offset):
    """
    Shifting the profiles by a constant factor in Fourier space
    """
    ##### Shifting pulse profiles through a shifted FT (see Deepto's 7/20/2020 email)
    if nbins % 2 == 0:
        fft_x = np.array(list(np.arange(int(nbins/2)+1)) + list(np.arange(int(nbins/2)-1) - (int(nbins/2)-1)))
    else:
        fft_x = np.array(list(np.arange(int(nbins/2)+1)) + list(np.arange(int(nbins/2)) - int(nbins/2)))

    shift = np.exp(-2j*np.pi*fft_x*offset/nbins)
    shifted_prof = np.real(np.fft.ifft(np.fft.fft(prof)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile
    shifted_err = np.real(np.fft.ifft(np.fft.fft(prof_err)*shift)) #taking the real component of the inverse transform of the shifted Fourier transform of the original folded profile

    return shift,shifted_prof,shifted_err

def update_fits(eventfile,phases):
    """
    Updating the event file with the orbital phases
    """
    events = fits.open(eventfile,mode='update')
    if 'ORBPHASE' in events[1].columns.names:
        print("ORBPHASE exists already; so overwriting...")
        events[1].data['ORBPHASE'] = phases
    else:
        print("ORBPHASE does not exist, so adding it into the event file")
        datacol = [fits.ColDefs( [fits.Column(name='ORBPHASE',format="D",array=phases) ])]
        cols = events[1].columns
        cols = cols + datacol[0] #adding in the info for ORBPHASE
        bt = fits.BinTableHDU.from_columns(cols,header=events[1].header,name=events[1].name)
        events[1] = bt

    events.flush()

def ftselect_phase(instrument):
    #phase_range = np.arange(0,1.01,0.05)
    phase_range = np.array([0.25,0.8,0.8,0.9,0.9,0.1,0.1,0.25])
    ftselect_txt = open('/Volumes/Samsung_T5/NGC300_XMMdata/ftselect_onoffiegress.txt','w')
    for i in range(0,len(phase_range)-1,2):
        input_events = '/Volumes/Samsung_T5/NGC300_XMMdata/ngc300x1_' + instrument + '.evt'
        output_file = '/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/ngc300x1_' + instrument + '_' + str(int(phase_range[i]*100)).zfill(3) + '-' + str(int(phase_range[i+1]*100)).zfill(3) + '.evt'
        #subprocess.run(['ftselect',input_events+'[events]',output_file,'ORBPHASE >= ' + str(phase_range[i]) + ' && ORBPHASE < ' + str(phase_range[i+1]),'clobber=yes'])
        ftselect_txt.write('ftselect ' + input_events+'[events]' + ' ' + output_file + ' ' + 'ORBPHASE >= ' + str(phase_range[i]) + ' && ORBPHASE < ' + str(phase_range[i+1]) + ' ' + 'clobber=yes' + '\n')
    ftselect_txt.close()

    return

def evselect_instructions(instrument):
    """
    Create a set of instructions to extract spectra from the phase-resolved event files
    """
    if instrument != 'mos1' and instrument != 'mos2' and instrument != 'pn':
        raise ValueError('Make sure that the instrument is either mos1, mos2, or pn!')

    instructions = open('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/evselect_' + instrument + '_instructions.txt','w')
    eventfiles = glob.glob('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/*' + instrument+'*.evt')

    for i in range(len(eventfiles)):
        if instrument == 'mos1' or instrument == 'mos2':
            instructions.write('evselect table='+eventfiles[i]+' withspectrumset=yes spectrumset=' + eventfiles[i][:-3]+'pha energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999')
        elif instrument == 'pn':
            instructions.write('evselect table='+eventfiles[i]+' withspectrumset=yes spectrumset=' + eventfiles[i][:-3]+'pha energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479')

        instructions.write('\n')

def grppha_instructions(instrument):
    """
    Create a set of instructions to introduce a grouping column onto the phase-resolved spectral files
    """
    if instrument != 'mos1' and instrument != 'mos2' and instrument != 'pn':
        raise ValueError('Make sure that the instrument is either mos1, mos2, or pn!')

    instructions = open('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/grppha_' + instrument + '_instructions.txt','w')
    spectra = glob.glob('/Volumes/Samsung_T5/NGC300_XMMdata/phase_res_spec/ngc300x1*' + instrument+'*.pha')

    for i in range(len(spectra)):
        spectrum_file = str(pathlib.Path(spectra[i]).name)
        if instrument == 'pn':
            instructions.write('grppha infile="' + spectrum_file + '" outfile="grp_' + spectrum_file + '" chatter=0 comm="group pn_channels_to_group.dat & systematics 60 2000 0.02 & exit"')
        elif instrument == 'mos1' or instrument == 'mos2':
            instructions.write('grppha infile="' + spectrum_file + '" outfile="grp_' + spectrum_file + '" chatter=0 comm="group mos_channels_to_group.dat & systematics 60 2000 0.02 & exit"')

        instructions.write('\n')

if __name__ == "__main__":
    #plot_profs('both')

    print('Step and Ramp model for PN:')
    step_n_ramp(phases,pn_prof,pn_prof_err,0.225,1.675,np.array([0.65,0.75,1,1.25,0.016,0.0035]))

    #print('Step and Ramp model for MOS1:')
    #step_n_ramp(phases,mos1_prof,mos1_prof_err,0.225,1.675,np.array([0.67,0.8,1,1.25,0.0045,0.001]))

    #print('Step and Ramp model for MOS2:')
    #step_n_ramp(phases,mos2_prof,mos2_prof_err,0.225,1.675,np.array([0.625,0.75,1,1.25,0.0045,0.001]))

    plt.show()

    pn_v2 = 0.735155769390959
    pn_v2_err = 0.0013104756938545274
    pn_v3 = 0.9734922065666258
    pn_v3_err = 0.0025331177762865694

    eclipse_ph = (pn_v2 + pn_v3)/2
    eclipse_ph_err = np.sqrt( (pn_v2_err/2)**2 + (pn_v3_err/2)**2 )
    #print(eclipse_ph,eclipse_ph_err)

    Porb_days = (1/8.4712e-6)/86400
    Porb_days_err = (0.00296e-6/(8.4712e-6)**2)/86400 #make sure to change!! Recall this is from FR/RSS simulations with 40% completeness
    #print(Porb_days)
    #print(Porb_days_err)
    #print(Porb_days,Porb_days_err)
    T0_xmm_MJD = fits.open(eventfile_pn)[1].header['MJDREF'] +(fits.open(eventfile_pn)[1].data['TIME'][0]+eclipse_ph*1/8.4712e-6)/86400
    T0_xmm_err = (eclipse_ph_err*1/8.4712e-6)/86400

    T0_swift_MJD = 58239.3498
    T0_swift_err = 0.0208

    no_cycles = (T0_swift_MJD - T0_xmm_MJD)/Porb_days
    no_cycles_err = np.sqrt( (T0_swift_err/Porb_days)**2 + (T0_xmm_err/Porb_days)**2 + ((T0_swift_err-T0_xmm_err)/Porb_days**2)*Porb_days_err**2 )
    #print(no_cycles,no_cycles_err)

    shift,shift_pn_prof,shift_pn_err = shift_profs(20,pn_prof[:20],pn_prof_err[:20],eclipse_ph)
    #shift,shift_mos1_prof,shift_mos1_err = shift_profs(20,mos1_prof[:20],mos1_prof_err[:20],eclipse_ph)
    #shift,shift_mos2_prof,shift_mos2_err = shift_profs(20,mos2_prof[:20],mos2_prof_err[:20],eclipse_ph)

    """
    plt.errorbar(x=phases,y=np.array(list(shift_pn_prof)*2),yerr=np.array(list(shift_pn_err)*2),drawstyle='steps-mid')
    #plt.errorbar(x=phases,y=np.array(list(shift_mos1_prof)*2),yerr=np.array(list(shift_mos1_err)*2),drawstyle='steps-mid')
    #plt.errorbar(x=phases,y=np.array(list(shift_mos2_prof)*2),yerr=np.array(list(shift_mos2_err)*2),drawstyle='steps-mid')
    plt.title('Exposure-corrected (using Stingray fold_events), shifted by factor of ' + str(eclipse_ph),fontsize=12)
    plt.xlabel('Phase',fontsize=12)
    plt.ylabel('Counts/s',fontsize=12)
    #plt.legend(('PN','MOS1','MOS2'),loc='best',fontsize=12)
    plt.show()
    """

    #### Updating the FITS files with the orbital phase!
    #update_fits(eventfile_pn,phases_pn)
    #update_fits(eventfile_mos1,phases_mos1)
    #update_fits(eventfile_mos2,phases_mos2)

    #ftselect_phase('pn')
    #ftselect_phase('mos1')
    #ftselect_phase('mos2')

    #### CHECKING THE INCLUSION OF ORBPHASE INTO THE EVENT FILES
    #newphases_pn = fits.open(eventfile_pn)[1].data['ORBPHASE']
    #newphases_mos1 = fits.open(eventfile_mos1)[1].data['ORBPHASE']
    #newphases_mos2 = fits.open(eventfile_mos2)[1].data['ORBPHASE']

    #newprof_pn,bins = np.histogram(newphases_pn,bins=np.linspace(0,1,21))
    #newprof_mos1,bins = np.histogram(newphases_mos1,bins=np.linspace(0,1,21))
    #newprof_mos2,bins = np.histogram(newphases_mos2,bins=np.linspace(0,1,21))
    #plt.errorbar(bins[:-1],newprof_pn,drawstyle='steps-mid')
    #plt.errorbar(bins[:-1],newprof_mos1,drawstyle='steps-mid')
    #plt.errorbar(bins[:-1],newprof_mos2,drawstyle='steps-mid')
    #plt.show()

    ##### Get the phase offset between Swift eclipse time and XMM's first event time:
    xmm_first_pn = MJDREF + times_pn[0]/86400
    xmm_ecl = T0_swift_MJD - int(no_cycles)*Porb_days #time of the mid-eclipse BEFORE the first XMM event
    if xmm_ecl > xmm_first_pn:
        xmm_ecl -= Porb_days
    phaseoff_pn = (xmm_first_pn-xmm_ecl)/Porb_days
    #print(phaseoff_pn)

    ##### Appending background file,
    #evselect_instructions('pn')
    #evselect_instructions('mos1')
    #evselect_instructions('mos2')

    ##### Adding the grouping column to the spectral files
    #grppha_instructions('pn')
    #grppha_instructions('mos1')
    #grppha_instructions('mos2')
