#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 8 9:30pm 2019

Analyzing the parameter values that come out of the XSPEC fitting.
This is a re-write of the old script which magically disappeared after
restarting the laptop a few times...

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from mpldatacursor import datacursor
import Lv0_dirs
import subprocess
import glob
from tqdm import tqdm

Lv0_dirs.global_par()

def model_par(model):
    """
    Given the model, return a list of associated Parameters

    model - name of the model used
    """
    if model == 'powerlaw':
        return ['PhoIndex','norm']
    if model == 'bbodyrad':
        return ['kT','norm']
    if model == 'ezdiskbb':
        return ['T_max','norm']
    if model == 'diskbb':
        return ['Tin','norm']
    if model == 'cutoffpl':
        return ['PhoIndex','HighECut','norm']
    if model == 'fdcut':
        return ['cutoffE','foldE']
    if model == 'relxill':
        return ['gamma','logxi','refl_frac','norm']
    if model == 'gauss':
        return ['LineE','Sigma','norm']
    if model == 'laor':
        return ['lineE','norm']
    else:
        raise ValueError("This is not a model pre-entered into model_par!")

def xspec_par(model,E1,E2):
    """
    From the model.txt files, obtain the XSPEC parameter fit values

    model - name of the model used
    E1 - lower bound for energy (4-digit PI string)
    E2 - upper bound for energy (4-digit PI string)
    """
    if 'simpl' in model:
        models = model.split('-')[2:]
    else:
        models = model.split('-')[1:] #sliced to take away the absorption model

    textfile = Lv0_dirs.NGC300 + 'spectral_fit_' + E1+'-'+E2+'/'+model+'.txt'
    contents = open(textfile,'r').read().split('\n')[4:]

    par_values = {}
    for k in range(len(models)):
        model_pars = model_par(models[k])
        for i in range(len(model_pars)): #for each parameter in the model
            par_values[models[k]+'-'+model_pars[i]] = []
            par_values[models[k]+'-'+model_pars[i]+'_unc'] = []
            for j in range(len(contents)):
                if (model_pars[i] in contents[j].split()) and (models[k] in contents[j]):
                    line = contents[j].split()
                    par_values[models[k]+'-'+model_pars[i]].append(float(line[-3])) #parameter value
                    try:
                        par_values[models[k]+'-'+model_pars[i]+'_unc'].append(float(line[-1])) #corresponding uncertainty
                    except ValueError:
                        par_values[models[k]+'-'+model_pars[i]+'_unc'].append(0) #corresponding uncertainty
    return par_values

def plot_ufspectra(model,MJDs,E1,E2):
    """
    Plotting the outputs from ufspectra!

    model - name of the model used
    MJDs - list of MJDs used
    E1 - lower bound for energy (4-digit PI string)
    E2 - upper bound for energy (4-digit PI string)
    """
    no_models = model.count('-') #number of models, e.g., tbnew-bbodyrad-powerlaw has TWO models!
    if 'simpl' in model:
        no_models -= 1 #subtract the convolved models... still need to generalize this somehow
        models = model.split('-')[2:]
    else:
        models = model.split('-')[1:] #sliced to take away the absorption model

    xspec_fits = xspec_par(model,E1,E2)
    no_vars = []
    for i in range(len(models)):
        no_vars.append(len(model_par(models[i])))
    no_cols = np.max(np.array(no_vars)) #maximum number of columns for the plotting

    input_file = Lv0_dirs.NGC300 + 'spectral_fit_'+E1+'-'+E2 + '/'+ model + '_ufspectra.txt'
    data = np.array(open(input_file,'r').read().split('\n'))

    if len(models) > 1:
        multiplier = 4+len(models) + 1
    else:
        multiplier = 4+len(models)

    separator = list(np.where(data==('NO '*multiplier)[:-1])[0])
    separator.insert(0,2)
    separator.insert(len(separator),len(data)-1)

    for i in range(len(MJDs)):
        spectrum_file = Lv0_dirs.NGC300+'spectral_fit_'+E1+'-'+E2+'/indiv_spectra/'+model+'_ufspectra_'+MJDs[i]+'.txt'
        spectrum = open(spectrum_file,'w')
        for j in range(separator[i]+1,separator[i+1]):
            spectrum.write(data[j]+'\n')
        spectrum.close()

    pdf_name = Lv0_dirs.NGC300 + 'spectral_fit_' + E1+'-'+E2+'/' + model + '_ufspectra.pdf'
    with PdfPages(pdf_name) as pdf:
        print('Plotting the ufspectra per MJD')
        for i in tqdm(range(len(MJDs))):
            fig = plt.figure(figsize=(16,9))
            gs = GridSpec(1+no_models,no_cols,width_ratios=[3]*no_cols,height_ratios=no_models*[2]+[8],wspace=0.25) #that '1' is for the overall spectra

            for j in range(len(models)):
                for k in range(len(model_par(models[j]))):
                    sing_model = models[j]
                    var = model_par(sing_model)[k]
                    ax = fig.add_subplot(gs[j,k])
                    ax.errorbar(x=np.array(MJDs,dtype='float'),y=xspec_fits[sing_model+'-'+var],yerr=xspec_fits[sing_model+'-'+var+'_unc'],fmt='x-')
                    ax.set_ylim([0.8*np.min(xspec_fits[sing_model+'-'+var]),1.2*np.max(xspec_fits[sing_model+'-'+var])])
                    ax.set_ylabel(sing_model + '\n' + var,fontsize=12)
                    ax.set_xlabel('Time (MJD)',fontsize=12)

            uftxt = Lv0_dirs.NGC300 + 'spectral_fit_' + E1+'-'+E2+'/indiv_spectra/'+model+'_ufspectra_' + str(int(MJDs[i])) + '.txt'
            if no_models == 1:
                total_no_models = 4 + no_models
            elif no_models > 1:
                total_no_models = 4 + 1 + no_models #the 1 is for the 'total'
            ufdata = np.genfromtxt(uftxt,skip_footer=1,usecols=(tuple(range(total_no_models))),unpack=True) #first 4 are "E, E_err, flux, flux_err"

            ax_spectrum = fig.add_subplot(gs[-1,:])
            ax_spectrum.errorbar(x=ufdata[0],y=ufdata[2],xerr=ufdata[1],yerr=ufdata[3],fmt='+') #plotting the data flux
            ax_spectrum.set_xlabel('Energy, E (keV)',fontsize=12)
            ax_spectrum.set_ylabel('Flux (photons/cm^2/s/keV)',fontsize=12)
            plt.yscale('log')
            plt.xscale('log')
            ax_spectrum.set_xlim([0.4,5])
            ax_spectrum.set_ylim([1e-6,1e-3])
            if no_models == 1:
                ax_spectrum.errorbar(x=ufdata[0],y=ufdata[4]) #for the total
                ax_spectrum.legend(('Data','Total/'+models[0]),loc='best',fontsize=12)
            elif no_models > 1:
                ax_spectrum.errorbar(x=ufdata[0],y=ufdata[4]) #for the total
                for m in range(5,len(ufdata)):
                    ax_spectrum.errorbar(x=ufdata[0],y=ufdata[m])
                ax_spectrum.legend(('Data','Total')+tuple([models[m] for m in range(len(models))]),loc='best',fontsize=12)

            for j in range(len(models)):
                for k in range(len(model_par(models[j]))):
                    sing_model = models[j]
                    var = model_par(sing_model)[k]
                    ax = fig.add_subplot(gs[j,k])
                    ax.errorbar(x=np.array(MJDs,dtype='float')[i],y=xspec_fits[sing_model+'-'+var][i],fmt='+',markersize=15)
                    ax.set_ylim([0.8*np.min(xspec_fits[sing_model+'-'+var]),1.2*np.max(xspec_fits[sing_model+'-'+var])])
                    ax.set_ylabel(sing_model + '\n' + var)

            pdf.savefig()
            plt.close()

    return

def plot_rrcl(model,MJDs,E1,E2,plottype):
    """
    To plot either of ratio, residual, chi, or lumin

    model - name of the model used
    MJDs - list of MJDs used
    E1 - lower bound for energy (4-digit PI string)
    E2 - upper bound for energy (4-digit PI string)
    plottype - either of 'ratio','resid','chi', or 'lumin'
    """
    if plottype != 'ratio' and plottype != 'resid' and plottype != 'chi' and plottype != 'lumin':
        raise ValueError('plottype should be either of ratio, resid, chi, or lumin!')

    input_file = Lv0_dirs.NGC300 + 'spectral_fit_'+E1+'-'+E2+'/'+model+'_'+plottype+'.txt'

    ##### getting the data
    if plottype == 'chi':
        data = np.array(open(input_file,'r').read().split('\n'))
        separator = list(np.where(data=='NO NO NO')[0])
        separator.insert(0,2)
        separator.insert(len(separator),len(data)-1)

        for i in range(len(MJDs)):
            spectrum_file = Lv0_dirs.NGC300+'spectral_fit_'+E1+'-'+E2+'/indiv_'+plottype+'/'+model+'_'+plottype+'_'+MJDs[i]+'.txt'
            spectrum = open(spectrum_file,'w')
            for j in range(separator[i]+1,separator[i+1]):
                spectrum.write(data[j]+'\n')
            spectrum.close()

    if plottype == 'ratio' or plottype == 'resid':
        data = np.array(open(input_file,'r').read().split('\n'))
        separator = list(np.where(data=='NO NO NO NO')[0])
        separator.insert(0,2)
        separator.insert(len(separator),len(data)-1)

        for i in range(len(MJDs)):
            spectrum_file = Lv0_dirs.NGC300+'spectral_fit_'+E1+'-'+E2+'/indiv_'+plottype+'/'+model+'_'+plottype+'_'+MJDs[i]+'.txt'
            spectrum = open(spectrum_file,'w')
            for j in range(separator[i]+1,separator[i+1]):
                spectrum.write(data[j]+'\n')
            spectrum.close()

    if plottype == 'lumin':
        contents = open(input_file).read().split('\n')
        lumin_lines = np.array([float(contents[i].split(' ')[2]) for i in range(2,len(contents),3)])

    pdf_file = Lv0_dirs.NGC300 + 'spectral_fit_'+E1+'-'+E2+'/' + model + '_' + plottype + '.pdf'
    with PdfPages(pdf_file) as pdf:
        if plottype == 'lumin':
            plt.errorbar(MJDs,lumin_lines/1e39,fmt='x',markersize=12)
            plt.xlabel('Time (MJD)',fontsize=12)
            plt.ylabel('Luminosity in 10^(39) erg/s',fontsize=12)
            plt.axhline(y=1,lw=0.5,alpha=0.8,color='r')
            plt.legend(('NGC ULX-1 Model Luminosity','Putative ULX Luminosity Threshold'),loc='best',fontsize=12)
            pdf.savefig()
            plt.close()
        if plottype == 'ratio' or plottype == 'resid':
            textfiles = sorted(glob.glob(Lv0_dirs.NGC300 + 'spectral_fit_'+E1+'-'+E2+'/indiv_' + plottype + '/' + model + '_' + plottype + '*.txt'))
            for i in range(len(textfiles)):
                E,E_err,r,r_err = np.genfromtxt(textfiles[i],usecols=(0,1,2,3),unpack=True) #extract energy, error in energy, ratio/residual and error in ratio/residual
                plt.errorbar(x=E,y=r,xerr=E_err,yerr=r_err,fmt='+')
                plt.title('MJD: ' + str(MJDs[i]),fontsize=12)
                plt.xlabel('Energy (keV)',fontsize=12)
                if plottype == 'ratio':
                    plt.ylabel('ratio (data/model)',fontsize=12)
                    plt.axhline(y=1,lw=0.5,alpha=0.8)
                if plottype == 'resid':
                    plt.ylabel('residuals',fontsize=12)
                    plt.axhline(y=0,lw=0.5,alpha=0.8)
                pdf.savefig()
                plt.close()
        if plottype == 'chi':
            textfiles = sorted(glob.glob(Lv0_dirs.NGC300 + 'spectral_fit_'+E1+'-'+E2+'/indiv_' + plottype + '/' + model + '_' + plottype + '*.txt'))
            for i in range(len(textfiles)):
                E,E_err,chi = np.genfromtxt(textfiles[i],usecols=(0,1,2),unpack=True)
                plt.step(x=E,y=chi)
                plt.axhline(y=0,lw=0.5,alpha=0.8)
                plt.title('MJD: ' + str(MJDs[i]),fontsize=12)
                plt.xlabel('Enegy (keV)',fontsize=12)
                plt.ylabel('sign(data-model) x Delta(chi^2)',fontsize=12)
                pdf.savefig()
                plt.close()

    return

def plot_HID(MJDs):
    """
    Plot the soft color-intensity diagram for a given set of MJDs

    MJDs - list of MJDs used
    """
    mjd_data,soft,soft_err,intensity,intensity_err = np.genfromtxt(Lv0_dirs.NGC300+'soft_color_HID.txt',skip_header=1,usecols=(0,1,2,3,4),unpack=True)

    soft_trunc = [soft[i] for i in range(len(mjd_data)) if str(int(mjd_data[i])) in MJDs]
    soft_err_trunc = [soft_err[i] for i in range(len(mjd_data)) if str(int(mjd_data[i])) in MJDs]
    intensity_trunc = [intensity[i] for i in range(len(mjd_data)) if str(int(mjd_data[i])) in MJDs]
    intensity_err_trunc = [intensity_err[i] for i in range(len(mjd_data)) if str(int(mjd_data[i])) in MJDs]

    pdf_name = Lv0_dirs.NGC300 + 'spectral_fit_' + E1 + '-' + E2 + '/soft_color_HID.pdf'
    with PdfPages(pdf_name) as pdf:
        for i in tqdm(range(len(soft_trunc))):
            plt.figure(figsize=(16,9))
            plt.xlim([-0.1,2.9])
            plt.ylim([-0.9,1.6])
            plt.xlabel('Soft Color: 1.0-2.0 keV/0.4-1.0 keV',fontsize=12)
            plt.ylabel('Intensity (ct/s)',fontsize=12)
            plt.title('MJD: ' + str(MJDs[i]))
            plt.errorbar(x=soft_trunc,y=intensity_trunc,xerr=soft_err_trunc,yerr=intensity_err_trunc,fmt='^')
            plt.errorbar(x=soft_trunc[i],y=intensity_trunc[i],xerr=soft_err_trunc[i],yerr=intensity_err_trunc[i],fmt='+',color='k',markersize=15)
            pdf.savefig()
            plt.close()

    return

def lumin_plus_par(model,MJDs,E1,E2):
    """
    Plot luminosity against a spectral parameter
    11/17: Need to be able to generalize such that I can use this function for >1 model!

    model - name of the model used
    MJDs - list of MJDs used
    E1 - lower bound for energy (4-digit PI string)
    E2 - upper bound for energy (4-digit PI string)
    """
    input_file = Lv0_dirs.NGC300 + 'spectral_fit_'+E1+'-'+E2+'/'+model+'_lumin.txt'
    contents = open(input_file).read().split('\n')
    lumin_lines = [float(contents[i].split(' ')[2]) for i in range(2,len(contents),3)]

    xspec_fits = xspec_par(model,E1,E2)

    fig = plt.figure(figsize=(16,9))
    fig.suptitle('Parameters from ' + model + ' against luminosity')
    gs = GridSpec(2,1,hspace=0) #that '1' is for the overall spectra

    singular_model = model.split('-')[-1] #PLACEHOLDER FOR NOW
    for k in range(len(model_par(singular_model))):
        var = model_par(singular_model)[k]
        ax = fig.add_subplot(gs[k,0])
        ax.errorbar(x=np.array(lumin_lines)/1e39,y=xspec_fits[singular_model+'-'+var],yerr=xspec_fits[singular_model+'-'+var+'_unc'],fmt='x')#,label=MJDs[j])
        #datacursor(formatter='{label}'.format,bbox=None)
        ax.set_ylim([0.8*np.min(xspec_fits[singular_model+'-'+var]),1.2*np.max(xspec_fits[singular_model+'-'+var])])
        ax.set_ylabel(singular_model + '\n' + var)
        ax.set_xlabel('Luminosity in 10^(39) erg/s',fontsize=12)

    plt.show()


if __name__ == "__main__":
    E1 = '0040'
    E2 = '0500'
    absorption = 'tbnew'
    models = ['-powerlaw','-bbodyrad','-ezdiskbb','-diskbb']
    MJDs_25 = ['58239','58244','58249','58254','58259','58264','58269','58274',
                '58289','58309','58314','58324','58329','58334','58339','58389',
                '58399','58449','58454','58459','58464','58484','58489','58504',
                '58509']
    MJDs_23 = ['58239','58244','58249','58254','58259','58264','58269','58274',
                '58289','58309','58314','58324','58329','58334','58339','58389',
                '58449','58454','58464','58484','58489','58504','58509']
    MJDs_21 = ['58239','58244','58249','58254','58259','58264','58269','58274',
                '58309','58314','58324','58329','58334','58339','58389','58449',
                '58454','58484','58489','58504','58509']

    #for i in range(len(models)):
    #    plot_rrcl(absorption+models[i],MJDs_25,E1,E2,'ratio')
    #plot_rrcl(absorption+'-powerlaw-bbodyrad',MJDs_21,E1,E2,'ratio')
    #plot_rrcl(absorption+'-powerlaw-bbodyrad',MJDs_21,E1,E2,'resid')
    #plot_rrcl(absorption+'-powerlaw-bbodyrad',MJDs_21,E1,E2,'chi')
    #plot_rrcl(absorption+'-powerlaw-bbodyrad',MJDs_21,E1,E2,'lumin')
    plot_ufspectra(absorption+'-powerlaw',MJDs_21,E1,E2)
    #xspec_par('tbnew-powerlaw',E1,E2)

    ##### DO NOT DO PLOT_HID FOR MJDS_25!!! YOU WILL GET THE WRONG ANSWER
    #plot_HID(MJDs_21)
    #lumin_plus_par(absorption+'-ezdiskbb',MJDs_21,E1,E2)
