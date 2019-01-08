#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 19:41:48 2018

@author: masonng

To make a movie (using FFMpeg) by stitching together images of:
    
1) Count rate from each of the FPMs - need to use a colorbar to show counts! 
See if I can add in the image from the NICER Mission Guide

https://heasarc.gsfc.nasa.gov/docs/nicer/mission_guide/

https://pillow.readthedocs.io/en/5.3.x/reference/Image.html 
    
https://stackoverflow.com/questions/35692507/plot-several-image-files-in-matplotlib-subplots
    
"""
from __future__ import division
import datetime
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('Agg')
from tqdm import tqdm
import os
import matplotlib.image as mpimg
import pyqtgraph as pg
import time
from scipy import stats
from scipy import signal

timestart = time.time()

print(datetime.datetime.now())

work_dir = '/Users/masonng/Documents/MIT/Research/Nicer-Logistics/'
obsid = '1034070104'  #Cen X-3
obsname = 'Cen_X-3'

def open_fits(work_dir, obsid):
    ## use barycentered data
    
    event = work_dir + obsid + '/xti/event_cl/ni' + obsid + '_0mpu7_cl_bary.evt'
    event = fits.open(event) #see event.info() for each card
    #event[1].header for events; event[2].header for GTIs

    return event

def get_data(work_dir,obsid):
    
    event = open_fits(work_dir,obsid)
    
    #recall: PI/PI_FAST are energy discriminants; PI = slow chain, able to register
    #the data from when the photons are incident on the anodes to create an electron
    #cloud; PI = fast chain, 'takes in' the data over a much shorter timescale
    #PI_FAST is really secondary data - only useful for background discriminant
    pi = event[1].data['PI']
#    print 'Done with PI_FAST'
    times = event[1].data['TIME']
#    print 'Done with TIME'
    detid_data = event[1].data['DET_ID']
    counts = np.ones(len(pi))
    
    return times, pi, counts, detid_data

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

#def get_subregion_indices(work_dir,obsid):
#    """
#    Get the starting/ending indices (or for 1s bins, the time) for the subregions
#    in the light curves. Subregions = where you get >5 counts.
#    """
#    
#    times, pi, counts, detid_data = get_data(work_dir,obsid)
#    shifted_t = times-times[0] 
#    
#    ## rebin the counts into 1-sec intervals
#    t_bins = np.linspace(0,int(shifted_t[-1])+1,int(shifted_t[-1]+2)) #get 1-second bins
#    intensity_1s, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins)
#
#    ### get indices where the counts are > 5:
#    index_subreg = [0]
#    threshold = np.where(intensity_1s>=5)
#    for i in range(len(threshold[0])-1):
#        if threshold[0][i+1]-threshold[0][i] > 50:
#            index_subreg.append(threshold[0][i])
#            index_subreg.append(threshold[0][i+1])
#    index_subreg.append(threshold[0][-1])
#    
#    return index_subreg
        
def get_subregion_indices(work_dir,obsid):
    ### FOR NOW, CONSIDER 1s BINS! 
    gtis = open_fits(work_dir,obsid)[2].data
    data_starts = np.array([gtis[i][0] for i in range(len(gtis))])
    data_stops = np.array([gtis[i][1] for i in range(len(gtis))])
    
    starts = data_starts-data_starts[0]
    stops = data_stops-data_stops[0]
    
    ### Use the list of GTIs to get subregion indices!
    index_subreg = [0]
    for i in range(len(starts)-1):
        if starts[i+1]-stops[i] >= 50:
            #if the time between start of next subregion and the end of the
            #previous subregion is > 50s, consider that a new subregion!
            index_subreg.append(int(round(stops[i])))
            index_subreg.append(int(round(starts[i+1])))
    index_subreg.append(int(round(stops[-1])))
    
    return index_subreg

def detid():
    ## get the DET_IDs from 0-7, 10-17, 20-27, 30-37, 40-47, 50-57, 60-67
    detids = []
    for i in range(7):
        dets = np.linspace(10*i,10*i+7,8)
        for j in range(len(dets)):
            detids.append(int(dets[j]))
    
    return detids

def detid_str():
    detids = detid()
    detids_str = []
    for i in range(len(detids)):
        if detids[i]<10:
            detids_str.append('0'+str(detids[i]))
        else:
            detids_str.append(str(detids[i]))
            
    return detids_str

def counts_per_s(work_dir,obsid,doplot):
    """
    Output is a dictionary, which has as keys, the DET_ID, and correspondingly
    an array where each entry = counts per second at any given second.
    
    counts_dict = {'detector':[t=1,t=2,t=3,...], 'detector':[t=1,t=2,...]}
    """ 
    
    counts_dict = {}
    
    times,pi,counts,detid_data = get_data(work_dir,obsid)
    
    shifted_t = times-times[0]
    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])+1)
    detids = detid()
    
    for i in range(len(detids)): #for each FPM
        times,pi,counts,detid_data = get_data(work_dir,obsid)
        np.place(counts,detid_data!=np.int8(detids[i]),0) #replace not-detector by 0
        summed_data, bin_edges, binnumber = stats.binned_statistic(shifted_t,counts,statistic='sum',bins=t_bins)
        if detids[i] < 10:
            index = '0' + str(detids[i])
            counts_dict[index] = summed_data
        else:
            counts_dict[str(detids[i])] = summed_data

    if doplot==True:
        boundaries = get_subregion_indices(work_dir,obsid)
        
        from matplotlib.backends.backend_pdf import PdfPages
        count = 0
        for i in tqdm(range(0,len(boundaries),2)):
            xlim1 = boundaries[i]
            xlim2 = boundaries[i+1]
            count+=1
            subr_dir = work_dir+obsid+'/detectors/'
            ensure_dir(subr_dir)
            filename = subr_dir + 'SubR'+str(count)+'.pdf'
            with PdfPages(filename) as pdf:
                max_y = 0
                for i in range(len(detids)):
                    if detids[i] < 10:
                        index = '0' + str(detids[i])
                    else:
                        index = str(detids[i])
                    if max(counts_dict[index][xlim1:xlim2]) > max_y:
                        max_y = max(counts_dict[index][xlim1:xlim2])
               
                for i in range(len(detids)):    
                    if detids[i] < 10:
                        index = '0' + str(detids[i])
                    else:
                        index = str(detids[i])
                    plt.figure(figsize=(10,8))
                    plt.plot(t_bins[:-1][xlim1:xlim2],counts_dict[index][xlim1:xlim2],'-')
                    plt.title(obsname + ' ' + obsid + ': Counts/s vs time (s), DETID ' + index,fontsize=15)
                    plt.xlabel('Time (s)',fontsize=15)
                    plt.xlim([xlim1,xlim2])
                    plt.ylim([0,max_y])
                    plt.ylabel('Count/s',fontsize=15)
                    pdf.savefig()
                    plt.close()
            
    return counts_dict

#cenx3_obsid_list = ['0034070101','0034070102','0034070103','0034070104','1034070101','1034070102','1034070103','1034070104','1034070105','1034070106']
#for i in tqdm(range(len(cenx3_obsid_list))):
#    counts_per_s(work_dir,cenx3_obsid_list[i],True)

def det_coords(detector_id):
    """
    In direct, 1-1 correspondence as to how the NASA NICER team defined its detector array!
    """
    coord_dict = {'06':(0,0),'07':(0,1),'16':(0,2),'17':(0,3),'27':(0,4),'37':(0,5),'47':(0,6),'57':(0,7),
                  '05':(1,0),'15':(1,1),'25':(1,2),'26':(1,3),'35':(1,4),'36':(1,5),'46':(1,6),'56':(1,7),
                  '04':(2,0),'14':(2,1),'24':(2,2),'34':(2,3),'44':(2,4),'45':(2,5),'54':(2,6),'55':(2,7),
                  '03':(3,0),'13':(3,1),'23':(3,2),'33':(3,3),'43':(3,4),'53':(3,5),'66':(3,6),'67':(3,7),
                  '02':(4,0),'12':(4,1),'22':(4,2),'32':(4,3),'42':(4,4),'52':(4,5),'64':(4,6),'65':(4,7),
                  '01':(5,0),'11':(5,1),'21':(5,2),'31':(5,3),'41':(5,4),'51':(5,5),'62':(5,6),'63':(5,7),
                  '00':(6,0),'10':(6,1),'20':(6,2),'30':(6,3),'40':(6,4),'50':(6,5),'60':(6,6),'61':(6,7)}

    return coord_dict[detector_id]
    
def det_array(work_dir,obsid):
    """
    Time = length of observation in s
    Detector_x = number of FPMs in the 'x-direction'
    Detector_y = number of FPMs in the 'y-direction'
    
    Generates a time x detector_x x detector_y array! 
    """
    
    counts_dict = counts_per_s(work_dir,obsid,False)
    counts_all_time = np.zeros((len(counts_dict['00']),8,7))
    detids_str = detid_str()
    
    for i in tqdm(range(len(detids_str))):
        x = det_coords(detids_str[i])[1]
        y = det_coords(detids_str[i])[0]
        counts_all_time[:,x,y] = counts_dict[detids_str[i]] 
        
    return counts_all_time

#### MIGHT BE TOO COMPUTATIONALLY INTENSIVE? AND NOT AS VISUALLY USEFUL?
def plot_array(work_dir,obsid):
    counts_array = det_array(work_dir,obsid)
    boundaries = get_subregion_indices(work_dir,obsid)
    im = mpimg.imread(work_dir+"Nicer_Detector.png")

    count=0
    for i in range(0,len(boundaries)-1,2):
        count+=1
        subr_dir = work_dir+obsid+'/SubR'+str(count) +'/'
        ensure_dir(subr_dir)
        xlim1 = boundaries[i]
        xlim2 = boundaries[i+1]
        
        fig = plt.figure(figsize=(10,8))
        ax1 = fig.add_subplot(121)
        ax1.set_xlabel('x', fontsize=15)
        ax1.set_ylabel('y', fontsize=15)
        ax2 = fig.add_subplot(122)
        ax2.imshow(im)
        ax2.axis('off')
        img = ax1.imshow(counts_array[xlim1],cmap='gist_heat',vmin=0,vmax=np.max(counts_array[xlim1:xlim2]))
        plt.colorbar(img,ax=ax1,fraction=0.046,pad=0.04)

        for j in tqdm(range(1,len(counts_array[xlim1:xlim2]))):
            ax1.set_title('NICER Detector Array Heatmap showing \n count rate per second \n for t ='+str(xlim1)+'-'+str(xlim2)+'s. Now t='+str(xlim1+j), fontsize=12)
            img.set_data(counts_array[xlim1+j])#,cmap='gist_heat',vmin=0,vmax=np.max(counts_array[xlim1:xlim2]))
            fig.canvas.draw()
            
            plt.savefig(subr_dir+'/'+obsid+'_'+str(j),dpi=600,fmt='png')
            plt.clf()
            
            
#count = counts_per_s(work_dir,obsid,False)
#boundaries = get_subregion_indices(work_dir,obsid)
cube_data = det_array(work_dir,obsid)
###imv = pg.ImageView()
###imv.show()
###imv.setImage(cube_data)
##            
from PyQt4.QtGui import QApplication
from PyQt4.QtCore import QTimer
def startApp():
    import m1
    import m2
    wnd = createWindow()
    wnd.show()
import sys
app = QApplication(sys.argv)
#splash = createSplashScreen()
#splash.show()
#            
import numpy as np
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
#data = np.ones((230,10,10))
imv = pg.ImageView()
imv.setImage(cube_data)
imv.show()
##
##
roi = pg.ROI([0,0],[1,1],pen=pg.mkPen('r',width=2))
imv.addItem(roi)
def getcoordinates(roi):
    data2,xdata = roi.getArrayRegion(data,imv.imageItem,returnMappedCoords=True)
    print(np.round(xdata,1))
roi.sigRegionChangeFinished.connect(getcoordinates)
#
pg.QtGui.QApplication.exec_()
QTimer.singleShot(1, startApp) # call startApp only after the GUI is ready
#sys.exit(app.exec_())

#if __name__ == '__main__':
#    import sys
#    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
#        QtGui.QApplication.instance().exec_()


#import pyqtgraph as pg
#import numpy as np
#x = np.random.normal(size=1000)
#y = np.random.normal(size=1000)
#pg.plot(x, y, pen=None, symbol='o')  ## setting pen=None disables line drawing
#pg.QtGui.QApplication.exec_()
#a = counts_per_s(work_dir,obsid,False)


#event = open_fits(work_dir,obsid)
#min(event[1].data['DET_ID']),max(event[1].data['DET_ID']

timeend = time.time()

print (timeend-timestart)

### 
### OLD AND WAY TOO SLOW SORTING FOR DET_ARRAY

        #ask this on stackoverflow tomorrow?
        #say I have data of counts for each detector, at every second
        #how to reshape such that i get a cube where i have something like
        #indices telling me [time, x position, y position] = counts?

#    for i in range(len(detids_str)):
        
    
    
#    times,pi_fast,detid_data = get_data(work_dir,obsid)    
#    shifted_t = times-times[0]
#    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])+1
#    
#    binned_counts = counts_per_s(work_dir,obsid,False)
#    counts_all_time = np.zeros((len(t_bins),7,8))
#    for i in tqdm(range(len(counts_all_time))): #for every second in the subregion
#        for j in tqdm(range(len(detids_str))):
#            counts_all_time[i][det_coords(detids_str[j])] = counts_per_s(work_dir,obsid,False)[detids_str[j]][i]

#def det_array(work_dir,obsid):
#    """
#    Generate a 7x8 (R x C) empty array first, then replace each entry by the
#    counts later. MUST comply with how NICER defined its arrays!
#    """
#    
#    ## maybe instead of this slow AF for loop, try using binned_statistic??
#    detids_str = detid_str()
#    times,pi_fast,detid_data = get_data(work_dir,obsid)    
#    shifted_t = times-times[0]
#    t_bins = np.linspace(0,int(shifted_t[-1]),int(shifted_t[-1])+1
#    
#    binned_counts = counts_per_s(work_dir,obsid,False)
#    counts_all_time = np.zeros((len(t_bins),7,8))
#    for i in tqdm(range(len(counts_all_time))): #for every second in the subregion
#        for j in tqdm(range(len(detids_str))):
#            counts_all_time[i][det_coords(detids_str[j])] = counts_per_s(work_dir,obsid,False)[detids_str[j]][i]
#
#    return counts_all_time