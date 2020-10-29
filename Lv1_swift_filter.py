#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues May 5 2:43pm 2020

Filtering out Swift data.
5/5: Just doing filters based on region files

"""
from __future__ import division, print_function
import numpy as np
import numpy.ma as ma
from astropy.io import fits
from astropy.wcs import WCS,utils
from reproject import reproject_interp
from astropy.coordinates import SkyCoord
from reproject.mosaicking import find_optimal_celestial_wcs
import matplotlib.pyplot as plt
import Lv0_dirs,Lv1_data_gtis,Lv2_presto_subroutines,Lv2_mkdir
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage import gaussian_filter
import os
from tqdm import tqdm
import subprocess
import pathlib
import glob

Lv0_dirs.global_par() #obtaining the global parameters

def get_ra_dec(eventfile):
    """
    Obtain the RA_OBJ and DEC_OBJ corresponding to the observation!

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    """
    event = fits.open(eventfile)
    event_header = event[1].header

    return event_header['RA_OBJ'], event_header['DEC_OBJ']

def barycorr(eventfile,outfile,refframe,orbit_file,output_folder):
    """
    General function to perform the barycenter corrections for a Swift event file

    eventfile - path to the event file. Will extract ObsID from this for the NICER files.
    outfile - path to the output event file with barycenter corrections applied
    refframe - reference frame for barycenter corrections (usually ICRS)
    orbit_file - path to the orbit file of the observation
    output_folder - path to the folder where the outfile will be
    """
    obsid = eventfile[2:13]

    logfile = output_folder + 'barycorr_notes.txt'
    ra,dec = get_ra_dec(eventfile)

    with open(logfile,'w') as logtextfile:
        output = subprocess.run(['barycorr',eventfile,'outfile='+outfile,'orbitfiles='+orbit_file,'ra='+str(ra),'dec='+str(dec),'refframe='+str(refframe),'clobber=YES'],capture_output=True,text=True)
        logtextfile.write(output.stdout)
        logtextfile.write('*------------------------------* \n')
        logtextfile.write(output.stderr)
        logtextfile.close()

def filtering(eventlist,regtype,instructfile):
    """
    Filtering the original event file based on the defined region files from DS9

    eventlist - list of event files
    regtype - type of region (e.g., for NGC 300, there's ngc300ulx1, ngc300x1, ngc300bg)
    instructfile - where to save the instructions file
    """
    if type(eventlist) != np.array and type(eventlist) != list:
        raise TypeError('eventlist should either be an array or a list!')

    instruct = open(instructfile,'w')

    instruct.write('set mission swift' + '\n')
    instruct.write('set inst xrt' + '\n')

    for i in range(len(eventlist)):
        parent_folder = str(pathlib.Path(eventlist[i]).parent)
        filename = str(pathlib.Path(eventlist[i]).name)

        instruct.write('set datadir ' + parent_folder + '/' + '\n')

        for j in range(len(regtype)):
            instruct.write('read event ' + filename + '\n')
            instruct.write('filter region ' + parent_folder + '/' + regtype[j] + '.reg' + '\n') #generalize this later!!! Maybe use a list
            instruct.write('extract spectrum' + '\n')
            instruct.write('save spectrum ' + parent_folder + '/' + regtype[j] + '/' + filename[:-4] + '_' + regtype[j] + '.pha' + '\n')
            instruct.write('extract event' + '\n')
            instruct.write('save event ' + parent_folder + '/' + regtype[j] + '/' + filename[:-4] + '_' + regtype[j] + '.evt' + '\n')
            instruct.write('no' + '\n')
            instruct.write('clear region' + '\n')
            instruct.write('clear data' + '\n')

        instruct.write('\n')

    instruct.close()

def time_order(eventlist,ordered_text):
    """
    Return a list of event files in time order. The order of ObsIDs in the database
    aren't necessarily organized in order of date.

    eventlist - list of event files
    ordered_text - path to the text file that will list the time-ordered event files
    """
    if type(eventlist) != np.array and type(eventlist) != list:
        raise TypeError('eventlist should either be an array or a list!')

    start_times = [fits.open(eventlist[i])[1].header['TSTART'] for i in range(len(eventlist))]

    time_ordered = np.argsort(start_times)
    ordered_eventlist = np.array(eventlist)[time_ordered]

    ordered_file = open(ordered_text,'w')
    for i in range(len(ordered_eventlist)):
        ordered_file.write(ordered_eventlist[i] + '\n')
    ordered_file.close()

    return

def interm_time_order(eventlist,initialfile,ordered_text):
    """
    Return a list of event files in time order. The order of ObsIDs in the database
    aren't necessarily organized in order of date.

    This is for the intermediate event files.

    eventlist - list of event files
    initialfile - initialfile used for the merging; especially if it was generated from merging just the events
    ordered_text - text file listing the paths and extension for the merging
    """
    if type(eventlist) != np.array and type(eventlist) != list:
        raise TypeError('eventlist should either be an array or a list!')

    start_times = [fits.open(eventlist[i])[1].header['TSTART'] for i in range(len(eventlist))]

    time_ordered = np.argsort(start_times)
    ordered_eventlist = np.array(eventlist)[time_ordered]

    ordered_file = open(ordered_text,'w')
    for i in range(len(ordered_eventlist)):
        if i == 0:
            ordered_file.write(initialfile + '[GTI]' + '\n')
        else:
            ordered_file.write(ordered_eventlist[i] + '[GTI]' + '\n')
    ordered_file.close()

    return

def merging(event_filelist,interm_file,interm_filelist,merged_file):
    """
    Facilitating the merging of the total event file, with individual barycenter-corrected
    event files!

    event_filelist - input text file containing paths to the barycenter-corrected event files
    interm_file - intermediate event file
    interm_filelist - input text file containing paths to the INTERMEDIATE event files.
                      Intermediate means events have been merged, but not the GTI rows
    merged_file - path to the merged file!
    """
    subprocess.run(['ftmerge','@'+event_filelist,interm_file])
    subprocess.run(['ftmerge','@'+interm_filelist,merged_file])

    return

def fparkey(fitsfile,keyword,value):
    """
    Running FTOOLS' fparkey on some FITS file

    fitsfile - some FITS file
    keyword - the key, e.g., ANCRFILE
    value - what the input value is for the keyword
    """
    subprocess.run(['fparkey',value,fitsfile,keyword])

    return

def images(imagelist,product,mode,output_folder):
    """
    Combining images from a list of observations to produce either a) integrated image;
    b) do a time-lapse of the observations of a field (controlled by mode)

    imagelist - list of image files (should be FITS files)
    product - whether to create an "integrated" image or time-lapse "movie"
    mode - whether to "show" or "save" the product
    output_folder - output folder for products
    """
    if type(imagelist) != np.array and type(imagelist) != list:
        raise TypeError('imagelist should either be an array or a list!')
    if product != 'integrated' and product != 'movie':
        raise TypeError('product should either be "integrated" or "movie"! If there are new products, it should be updated.')
    if mode != 'save' and mode != 'show':
        raise TypeError('mode should either be "save" or "show"!')

    #ref_coord = SkyCoord('00h55m09.990s','-37d42m12.16s',frame='icrs')
    ref_coord1 = SkyCoord('00h55m30.0s','-37d47m0.0s',frame='icrs')
    ref_coord2 = SkyCoord('00h54m45.0s','-37d37m0.0s',frame='icrs')

    ngc300x1_coord = SkyCoord('00h55m9.6875s','-37d42m17.859s',frame='icrs')
    ngc300ulx1_coord = SkyCoord('00h55m04.8989s','-37d41m45.579',frame='icrs')
    ngc300bg_coord = SkyCoord('00h55m14.6298s','-37d39m57.287s',frame='icrs')
    nicer_fov_coord = SkyCoord('00h55m04.8989s','-37d41m45.579',frame='icrs') #same as ULX-1, because NICER observations were on the ULX

    #ref_coord1 = SkyCoord('00h55m15.0s','-37d43m0.0s',frame='icrs')
    #ref_coord2 = SkyCoord('00h55m00.0s','-37d41m0.0s',frame='icrs')
    obj_name = fits.open(imagelist[0])[0].header['OBJECT']

    if mode == 'save':
        parent_folder = str(pathlib.Path(imagelist[0]).parent)
        if product == 'integrated':
            ref_image = fits.open(imagelist[0])[0]
            base_image = np.zeros(np.shape(ref_image.data))
            base_exp = np.zeros(np.shape(ref_image.data))
            wcs = WCS(ref_image.header)

            ref_pixel1 = utils.skycoord_to_pixel(ref_coord1,wcs)
            ref_pixel2 = utils.skycoord_to_pixel(ref_coord2,wcs)

            ngc300x1_pixel = utils.skycoord_to_pixel(ngc300x1_coord,wcs)
            ngc300ulx1_pixel = utils.skycoord_to_pixel(ngc300ulx1_coord,wcs)
            ngc300bg_pixel = utils.skycoord_to_pixel(ngc300bg_coord,wcs)
            nicer_fov_pixel = utils.skycoord_to_pixel(nicer_fov_coord,wcs)

            for i in tqdm(range(len(imagelist))):
                imagefile = fits.open(imagelist[i])[0]
                expfile = fits.open(imagelist[i][:-6] + 'ex.img')[0]
                array_im,footprint = reproject_interp(imagefile,ref_image.header)
                array_ex,footprint = reproject_interp(expfile,ref_image.header)

                base_image += array_im
                base_exp += array_ex

            plt.figure()
            plt.subplot(projection=wcs)
            ngc300x1_circle = plt.Circle((ngc300x1_pixel[0],ngc300x1_pixel[1]),radius=30/2.36,color='y',lw=0.5,fill=False,label='NGC300 X-1')
            ngc300ulx1_circle = plt.Circle((ngc300ulx1_pixel[0],ngc300ulx1_pixel[1]),radius=35/2.36,color='b',lw=0.5,fill=False,label='NGC300 ULX-1')
            ngc300bg_circle = plt.Circle((ngc300bg_pixel[0],ngc300bg_pixel[1]),radius=120/2.36,color='m',lw=0.5,fill=False,label='NGC300 bg')
            nicer_fov_circle = plt.Circle((nicer_fov_pixel[0],nicer_fov_pixel[1]),radius=186/2.36,color='w',lw=0.5,fill=False,label='NICER FOV')

            plt.gcf().gca().add_artist(ngc300x1_circle)
            plt.gcf().gca().add_artist(ngc300ulx1_circle)
            plt.gcf().gca().add_artist(ngc300bg_circle)
            plt.gcf().gca().add_artist(nicer_fov_circle)

            #log_image = ma.log10(base_image)
            #plt.imshow(log_image.filled(0),vmin=0,vmax=np.log10(np.nanmax(base_image)),cmap='gist_heat')
            #plt.imshow((base_image/base_exp)/np.nanmax(base_image/base_exp),vmin=0,vmax=1,cmap='gist_heat')
            plt.imshow(gaussian_filter(base_image,1),vmin=0,vmax=10,cmap='gist_heat')
            plt.xlabel('Right Ascension (hh:mm:ss)',fontsize=12)
            plt.ylabel('Declination (deg)',fontsize=12)
            plt.xlim([ref_pixel1[0],ref_pixel2[0]])
            plt.ylim([ref_pixel1[1],ref_pixel2[1]])
            plt.legend([ngc300x1_circle,ngc300ulx1_circle,ngc300bg_circle,nicer_fov_circle],["NGC300 X-1 (30 arcsec)",'NGC300 ULX-1 (35 arcsec)','NGC300 bg (120 arcsec)','nicer_fov_pixel (3.1 arcmin)'])
            plt.colorbar().set_label('Counts')

            plt.show()

        if product == 'movie':
            #### saving all the images into one PDF file
            pdf_filename = output_folder + obj_name + '_field_movie.pdf'
            with PdfPages(pdf_filename) as pdf:
                for i in tqdm(range(len(imagelist))):
                    obsid = str(pathlib.Path(imagelist[i]).name)[:13]

                    fitsfile = fits.open(imagelist[i])[0]
                    expfile = fits.open(imagelist[i][:-6] + 'ex.img')[0]
                    date_obs = fitsfile.header['DATE-OBS']

                    wcs = WCS(fitsfile.header)
                    ref_pixel1 = utils.skycoord_to_pixel(ref_coord1,wcs)
                    ref_pixel2 = utils.skycoord_to_pixel(ref_coord2,wcs)

                    ngc300x1_pixel = utils.skycoord_to_pixel(ngc300x1_coord,wcs)
                    ngc300ulx1_pixel = utils.skycoord_to_pixel(ngc300ulx1_coord,wcs)
                    ngc300bg_pixel = utils.skycoord_to_pixel(ngc300bg_coord,wcs)
                    nicer_fov_pixel = utils.skycoord_to_pixel(nicer_fov_coord,wcs)

                    plt.figure(figsize=(16,9))
                    plt.subplot(projection=wcs)
                    norm = fitsfile.data/expfile.data
                    where_are_nans = np.isnan(norm)
                    norm[where_are_nans] = 0
                    plt.imshow(gaussian_filter(norm/np.nanmax(norm),sigma=0.5),vmin=0,vmax=1,cmap='gist_heat')
                    ngc300x1_circle = plt.Circle((ngc300x1_pixel[0],ngc300x1_pixel[1]),radius=30/2.36,color='y',lw=0.5,fill=False,label='NGC300 X-1')
                    ngc300ulx1_circle = plt.Circle((ngc300ulx1_pixel[0],ngc300ulx1_pixel[1]),radius=35/2.36,color='b',lw=0.5,fill=False,label='NGC300 ULX-1')
                    ngc300bg_circle = plt.Circle((ngc300bg_pixel[0],ngc300bg_pixel[1]),radius=120/2.36,color='m',lw=0.5,fill=False,label='NGC300 bg')
                    nicer_fov_circle = plt.Circle((nicer_fov_pixel[0],nicer_fov_pixel[1]),radius=186/2.36,color='w',lw=0.5,fill=False,label='NICER FOV')

                    plt.gcf().gca().add_artist(ngc300x1_circle)
                    plt.gcf().gca().add_artist(ngc300ulx1_circle)
                    plt.gcf().gca().add_artist(ngc300bg_circle)
                    plt.gcf().gca().add_artist(nicer_fov_circle)

                    plt.title('Observation Date: ' + str(date_obs) + ', ObsID: ' + obsid)
                    plt.xlabel('Right Ascension (hh:mm:ss)',fontsize=12)
                    plt.ylabel('Declination (deg)',fontsize=12)
                    plt.xlim([ref_pixel1[0],ref_pixel2[0]])
                    plt.ylim([ref_pixel1[1],ref_pixel2[1]])
                    plt.legend([ngc300x1_circle,ngc300ulx1_circle,ngc300bg_circle,nicer_fov_circle],["NGC300 X-1 (30 arcsec)",'NGC300 ULX-1 (35 arcsec)','NGC300 bg (120 arcsec)','nicer_fov_pixel (3.1 arcmin)'])
                    plt.colorbar().set_label('Counts/s (relative to maximum)')
                    pdf.savefig()
                    plt.close()

            #### saving each image into an individual file

            for i in tqdm(range(len(imagelist))):
                obsid = str(pathlib.Path(imagelist[i]).name)[:13]

                fitsfile = fits.open(imagelist[i])[0]
                expfile = fits.open(imagelist[i][:-6] + 'ex.img')[0]
                date_obs = fitsfile.header['DATE-OBS']

                wcs = WCS(fitsfile.header)
                ref_pixel1 = utils.skycoord_to_pixel(ref_coord1,wcs)
                ref_pixel2 = utils.skycoord_to_pixel(ref_coord2,wcs)

                ngc300x1_pixel = utils.skycoord_to_pixel(ngc300x1_coord,wcs)
                ngc300ulx1_pixel = utils.skycoord_to_pixel(ngc300ulx1_coord,wcs)
                ngc300bg_pixel = utils.skycoord_to_pixel(ngc300bg_coord,wcs)
                nicer_fov_pixel = utils.skycoord_to_pixel(nicer_fov_coord,wcs)

                plt.figure(figsize=(16,9))
                plt.subplot(projection=wcs)
                #plt.imshow(fitsfile.data,vmin=0,vmax=np.max(fitsfile.data),cmap='gist_heat')
                norm = fitsfile.data/expfile.data
                where_are_nans = np.isnan(norm)
                norm[where_are_nans] = 0
                plt.imshow(gaussian_filter(norm/np.nanmax(norm),sigma=0.5),vmin=0,vmax=1,cmap='gist_heat')
                ngc300x1_circle = plt.Circle((ngc300x1_pixel[0],ngc300x1_pixel[1]),radius=30/2.36,color='y',lw=0.5,fill=False,label='NGC300 X-1')
                ngc300ulx1_circle = plt.Circle((ngc300ulx1_pixel[0],ngc300ulx1_pixel[1]),radius=35/2.36,color='b',lw=0.5,fill=False,label='NGC300 ULX-1')
                ngc300bg_circle = plt.Circle((ngc300bg_pixel[0],ngc300bg_pixel[1]),radius=120/2.36,color='m',lw=0.5,fill=False,label='NGC300 bg')
                nicer_fov_circle = plt.Circle((nicer_fov_pixel[0],nicer_fov_pixel[1]),radius=186/2.36,color='w',lw=0.5,fill=False,label='NICER FOV')

                plt.gcf().gca().add_artist(ngc300x1_circle)
                plt.gcf().gca().add_artist(ngc300ulx1_circle)
                plt.gcf().gca().add_artist(ngc300bg_circle)
                plt.gcf().gca().add_artist(nicer_fov_circle)

                plt.title('Observation Date: ' + str(date_obs) + ', ObsID: ' + obsid)
                plt.xlabel('Right Ascension (hh:mm:ss)',fontsize=12)
                plt.ylabel('Declination (deg)',fontsize=12)
                plt.xlim([ref_pixel1[0],ref_pixel2[0]])
                plt.ylim([ref_pixel1[1],ref_pixel2[1]])
                plt.legend([ngc300x1_circle,ngc300ulx1_circle,ngc300bg_circle,nicer_fov_circle],["NGC300 X-1 (30 arcsec)",'NGC300 ULX-1 (35 arcsec)','NGC300 bg (120 arcsec)','nicer_fov_pixel (3.1 arcmin)'])
                plt.colorbar().set_label('Counts')
                plt.savefig(output_folder + obj_name + '_meanimage' + str(i).zfill(4) + '.png',format='png')
                plt.close()


    if mode == 'show':
        for i in range(len(imagelist)):
            fitsfile = fits.open(imagelist[i])[0]
            wcs = WCS(fitsfile.header)
            ref_pixel1 = utils.skycoord_to_pixel(ref_coord1,wcs)
            ref_pixel2 = utils.skycoord_to_pixel(ref_coord2,wcs)

            plt.subplot(projection=wcs)
            plt.imshow(fitsfile.data,vmin=0,vmax=np.max(fitsfile.data),cmap='gist_heat')
            plt.xlabel('Right Ascension (hh:mm:ss)',fontsize=12)
            plt.ylabel('Declination (deg)',fontsize=12)
            plt.xlim([ref_pixel1[0],ref_pixel2[0]])
            plt.ylim([ref_pixel1[1],ref_pixel2[1]])
            plt.colorbar()
            plt.show()


if __name__ == "__main__":
    swift_xrt_event = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/'
    filetype = 'ngc300_nicerfov'

    ancrfile = '/Volumes/Samsung_T5/swxpc0to12s6_20010101v013.arf'
    respfile = '/Volumes/Samsung_T5/swxpc0to12s6_20130101v014.rmf'

    nicer_ancrfile = '/Volumes/Samsung_T5/nicer-consim135p-teamonly-array50.arf'
    nicer_respfile = '/Volumes/Samsung_T5/nicer-rmf6s-teamonly-array50.rmf'

    ngc300_events = sorted(glob.glob(swift_xrt_event + '/sw*pc*po*cl*evt'))
    instructfile = swift_xrt_event + 'instructions.txt'

    #for i in tqdm(range(len(ngc300_events))):
    #    inputfile = ngc300_events[i]
    #    outputfile = str(pathlib.Path(ngc300_events[i]).parent) + '/' + str(pathlib.Path(ngc300_events[i]).name)[:20] + '_bary_cl.evt'
    #    orbitfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/auxil/sw' + str(pathlib.Path(inputfile).name)[2:13] + 'sao.fits'

    #    barycorr(inputfile,outputfile,'ICRS',orbitfile,swift_xrt_event)

    ngc300_events_bary = sorted(glob.glob(swift_xrt_event + 'sw*pc*po*bary*cl*evt'))
    #filtering(ngc300_events_bary,['ngc300ulx1','ngc300x1','ngc300bg'],instructfile)
    #filtering(ngc300_events_bary,['ngc300_nicerfov'],swift_xrt_event + 'nicerfov_instructions.txt')

    ngc300_type_events = sorted(glob.glob(swift_xrt_event + filetype + '/sw*pc*po*bary*cl*evt'))
    #time_order(ngc300_type_events,swift_xrt_event + filetype + '/eventfiles.list')
    #interm_time_order(ngc300_type_events,swift_xrt_event + filetype + '/' + filetype + '_intermediate.evt',swift_xrt_event + filetype + '/eventfiles_intermediate.list')

    #merging(swift_xrt_event + filetype + '/eventfiles.list',swift_xrt_event + filetype + '/' + filetype + '_intermediate.evt',swift_xrt_event + filetype + '/eventfiles_intermediate.list',swift_xrt_event + filetype + '/' + filetype + '_merge.evt')

    #merging('/Volumes/Samsung_T5/n300_ulx_2020/nicerdata_spectra/eventfiles.list','/Volumes/Samsung_T5/n300_ulx_2020/nicerdata_spectra/nicerdata_spectra_intermediate.evt','/Volumes/Samsung_T5/n300_ulx_2020/nicerdata_spectra/eventfiles_intermediate.list','/Volumes/Samsung_T5/n300_ulx_2020/nicerdata_spectra/nicerdata_spectra_merge.evt')

    #fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge.pha','BACKFILE','ngc300bg_merge.pha')
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge.pha','BACKFILE','ngc300bg_merge.pha')
    #fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge.pha','RESPFILE',respfile)

    #fparkey('/Volumes/Samsung_T5/n300_ulx_2020/nicerdata_spectra/nicerdata_spectra_merge.pha','ANCRFILE','/Volumes/Samsung_T5/nicer-consim135p-teamonly-array50.arf')
    #fparkey('/Volumes/Samsung_T5/n300_ulx_2020/nicerdata_spectra/nicerdata_spectra_merge.pha','RESPFILE','/Volumes/Samsung_T5//Volumes/Samsung_T5/nicer-rmf6s-teamonly-array50.rmf')

    #### First overlap

    ###merging(swift_xrt_event + 'ngc300x1/' + 'niceroverlap_all.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_all_int.evt', swift_xrt_event + 'ngc300x1/' + 'niceroverlap_all_int.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_all.evt')
    ###merging(swift_xrt_event + 'ngc300x1/' + 'niceroverlap_spec1.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_spec1_int.evt', swift_xrt_event + 'ngc300x1/' + 'niceroverlap_spec1_int.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_spec1.evt' )
    ###merging(swift_xrt_event + 'ngc300x1/' + 'niceroverlap_spec2.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_spec2_int.evt', swift_xrt_event + 'ngc300x1/' + 'niceroverlap_spec2_int.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_spec2.evt' )

    #merging(swift_xrt_event + 'ngc300bg/' + 'niceroverlap_all.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_all_int.evt', swift_xrt_event + 'ngc300bg/' + 'niceroverlap_all_int.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_all.evt')
    #merging(swift_xrt_event + 'ngc300bg/' + 'niceroverlap_spec1.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_spec1_int.evt', swift_xrt_event + 'ngc300bg/' + 'niceroverlap_spec1_int.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_spec1.evt' )
    #merging(swift_xrt_event + 'ngc300bg/' + 'niceroverlap_spec2.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_spec2_int.evt', swift_xrt_event + 'ngc300bg/' + 'niceroverlap_spec2_int.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_spec2.evt' )

    #merging(swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap_all.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_all_int.evt', swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap_all_int.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_all.evt')
    #merging(swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap_spec1.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec1_int.evt', swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap_spec1_int.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec1.evt' )
    #merging(swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap_spec2.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec2_int.evt', swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap_spec2_int.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec2.evt' )

    #merging(swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap_all.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_all_int.evt', swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap_all_int.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_all.evt')
    #merging(swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap_spec1.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_spec1_int.evt', swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap_spec1_int.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_spec1.evt' )
    #merging(swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap_spec2.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_spec2_int.evt', swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap_spec2_int.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_spec2.evt' )

    #merging(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap_all.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_all_intermediate.evt',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap_all_int.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_all.evt')
    #merging(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap_spec1.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_spec1_intermediate.evt',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap_spec1_int.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_spec1.evt')
    #merging(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap_spec2.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_spec2_intermediate.evt',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap_spec2_int.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_spec2.evt')

    #fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_all.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_spec1.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_spec2.pha','ANCRFILE',ancrfile)

    #fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_all.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_spec1.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap_spec2.pha','RESPFILE',respfile)

    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_all.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec1.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec2.pha','ANCRFILE',ancrfile)

    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_all.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec1.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec2.pha','RESPFILE',respfile)

    #fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_all.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_spec1.pha','RESPFILE',respfile)
    #fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_spec2.pha','RESPFILE',respfile)

    #fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_all.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_spec1.pha','ANCRFILE',ancrfile)
    #fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap_spec2.pha','ANCRFILE',ancrfile)

    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_all.pha','ANCRFILE',nicer_ancrfile)
    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_all.pha','RESPFILE',nicer_respfile)
    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_spec1.pha','ANCRFILE',nicer_ancrfile)
    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_spec1.pha','RESPFILE',nicer_respfile)
    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_spec2.pha','ANCRFILE',nicer_ancrfile)
    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap_spec2.pha','RESPFILE',nicer_respfile)

    #fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_all.pha','BACKFILE','ngc300bg_merge_niceroverlap_all.pha')
    #fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_spec1.pha','BACKFILE','ngc300bg_merge_niceroverlap_spec1.pha')
    #fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap_spec2.pha','BACKFILE','ngc300bg_merge_niceroverlap_spec2.pha')

    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_all.pha','BACKFILE','ngc300bg_merge_niceroverlap_all.pha')
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec1.pha','BACKFILE','ngc300bg_merge_niceroverlap_spec1.pha')
    #fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap_spec2.pha','BACKFILE','ngc300bg_merge_niceroverlap_spec2.pha')

    #### Second overlap
    #merging(swift_xrt_event + 'ngc300x1/' + 'niceroverlap2_spec1.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec1_int.evt', swift_xrt_event + 'ngc300x1/' + 'niceroverlap2_spec1_int.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec1.evt' )
    #merging(swift_xrt_event + 'ngc300x1/' + 'niceroverlap2_spec2.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec2_int.evt', swift_xrt_event + 'ngc300x1/' + 'niceroverlap2_spec2_int.list', swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec2.evt' )

    #merging(swift_xrt_event + 'ngc300bg/' + 'niceroverlap2_spec1.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap2_spec1_int.evt', swift_xrt_event + 'ngc300bg/' + 'niceroverlap2_spec1_int.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap2_spec1.evt' )
    #merging(swift_xrt_event + 'ngc300bg/' + 'niceroverlap2_spec2.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap2_spec2_int.evt', swift_xrt_event + 'ngc300bg/' + 'niceroverlap2_spec2_int.list', swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap2_spec2.evt' )

    #merging(swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap2_spec1.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec1_int.evt', swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap2_spec1_int.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec1.evt' )
    #merging(swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap2_spec2.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec2_int.evt', swift_xrt_event + 'ngc300ulx1/' + 'niceroverlap2_spec2_int.list', swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec2.evt' )

    #merging(swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap2_spec1.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap2_spec1_int.evt', swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap2_spec1_int.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap2_spec1.evt' )
    #merging(swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap2_spec2.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap2_spec2_int.evt', swift_xrt_event + 'ngc300_nicerfov/' + 'niceroverlap2_spec2_int.list', swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap2_spec2.evt' )

    #merging(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap2_spec1.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap2_spec1_intermediate.evt',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap2_spec1_int.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap2_spec1.evt')
    #merging(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap2_spec2.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap2_spec2_intermediate.evt',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/niceroverlap2_spec2_int.list',Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap2_spec2.evt')

    """
    fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec1.pha','ANCRFILE',ancrfile)
    fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec2.pha','ANCRFILE',ancrfile)

    fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec1.pha','RESPFILE',respfile)
    fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec2.pha','RESPFILE',respfile)

    fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap2_spec1.pha','ANCRFILE',ancrfile)
    fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap2_spec2.pha','ANCRFILE',ancrfile)

    fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap2_spec1.pha','RESPFILE',respfile)
    fparkey(swift_xrt_event + 'ngc300bg/ngc300bg_merge_niceroverlap2_spec2.pha','RESPFILE',respfile)

    fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec1.pha','ANCRFILE',ancrfile)
    fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec2.pha','ANCRFILE',ancrfile)

    fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec1.pha','RESPFILE',respfile)
    fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec2.pha','RESPFILE',respfile)

    fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap2_spec1.pha','RESPFILE',respfile)
    fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap2_spec2.pha','RESPFILE',respfile)

    fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap2_spec1.pha','ANCRFILE',ancrfile)
    fparkey(swift_xrt_event + 'ngc300_nicerfov/ngc300_nicerfov_merge_niceroverlap2_spec2.pha','ANCRFILE',ancrfile)
    """

    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap2_spec1.pha','ANCRFILE',nicer_ancrfile)
    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap2_spec1.pha','RESPFILE',nicer_respfile)
    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap2_spec2.pha','ANCRFILE',nicer_ancrfile)
    #fparkey(Lv0_dirs.NGC300_2020 + 'nicerdata_spectra/nicerdata_spectra_overlap2_spec2.pha','RESPFILE',nicer_respfile)

    """
    fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec1.pha','BACKFILE','ngc300bg_merge_niceroverlap2_spec1.pha')
    fparkey(swift_xrt_event + 'ngc300x1/ngc300x1_merge_niceroverlap2_spec2.pha','BACKFILE','ngc300bg_merge_niceroverlap2_spec2.pha')

    fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec1.pha','BACKFILE','ngc300bg_merge_niceroverlap2_spec1.pha')
    fparkey(swift_xrt_event + 'ngc300ulx1/ngc300ulx1_merge_niceroverlap2_spec2.pha','BACKFILE','ngc300bg_merge_niceroverlap2_spec2.pha')
    """

    """
    image_files = glob.glob('/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/products/sw*xpc*_sk.img')
    image_ordered_text = '/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/products/image_ordered.list'
    #time_order(image_files,image_ordered_text)

    time_ordered_images = open(image_ordered_text).read().split('\n')[:-1]
    images(time_ordered_images,'integrated','save','/Volumes/Samsung_T5/NGC300_ULX_Swift/xrt/event/images/')
    """
