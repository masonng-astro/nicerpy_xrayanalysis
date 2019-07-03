#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 1 2.06pm 2019

Obtaining the GTIs that are > 1000s.

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs,Lv0_call_eventcl,Lv0_call_ufa,Lv1_data_gtis
import subprocess

Lv0_dirs.global_par() #obtaining the global parameters

def desired_gtis(obsid,bary,gap,desired_length):
    """
    Obtain a list of tuples for GTIs, both shifted (i.e., starts at 0) and the actual
    GTIs (MJD), that have intervals longer than a desired length (e.g. 1000s)

    obsid - Observation ID of the object of interest (10-digit str)
    bary - Whether the data is barycentered. True/False
    gap - gap between data that would still be considered data
    -- for example, if you have data for t=1000-1040 and t=1045-1100, then if
    gap < 5, you get TWO GTIs. If gap >= 5, then you have 1 GTI from 1000-1100.
    This is different to desired_length - this is just meant to 'join up data'
    which are separated by mere seconds!
    desired_length - if GTIs > some desired length, use that data. This is used to
    weed out short GTIs.
    """
    if type(obsid) != str:
        raise TypeError("ObsID should be a string!")
    if bary != True and bary != False:
        raise ValueError("bary should either be True or False!")

    shifted_gtis,actual_gtis = Lv1_data_gtis.get_gtis(obsid,bary,gap)

    shifted_gtis_desired = [] #list of desired (shifted) GTIs that are > desired length
    actual_gtis_desired = [] #list of actual GTIs (MJD) that are > desired length
    for i in range(0,len(shifted_gtis),2):
        if (shifted_gtis[i+1]-shifted_gtis[i])>desired_length:
            shifted_gtis_desired.append((shifted_gtis[i],shifted_gtis[i+1]))
            actual_gtis_desired.append((actual_gtis[i],actual_gtis[i+1]))

    outputfile = open(Lv0_dirs.NICERSOFT_DATADIR+obsid+"_pipe/interval_times.txt","w")
    for i in range(len(actual_gtis_desired)):
        outputfile.write(str(actual_gtis_desired[i][0]) + ' ' + str(actual_gtis_desired[i][1]) + '\n')
    ### Using subprocess to break up the big .evt file!
#    output_file = open(Lv0_dirs.NICERSOFT_DATADIR+obsid+'_pipe/' + "ni"+str(obsid)+"_GTIs_"+str(desired_length)+"s.txt","w")
#    for i in range(len(actual_gtis_desired)):
#        output_file.write(actual_gtis_desired[i][0] + ' ' + actual_gtis_desired[i][1]+'\n')
#        infile = '"ni'+str(obsid)+'_nicersoft_bary.evt[TIME=' + str(actual_gtis_desired[i][0])+':'+str(actual_gtis_desired[i][1])+']"'
#        outfile = 'ni'+str(obsid)+'_nicersoft_bary_GTI'+str(i+1)+'.evt'

#        subprocess.call('heainit',shell=True,executable='/usr/local/bin/interactive_bash')
#        subprocess.call(['niextract-events',infile,outfile],shell=True)

    return shifted_gtis_desired,actual_gtis_desired

if __name__ == "__main__":
    shifted,actual = desired_gtis('1034090111',True,50,1000)
    for i in range(len(shifted)):
        print(shifted[i],actual[i])
