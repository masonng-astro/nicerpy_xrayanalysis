#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 8 2:11pm 2019

Extracting the GTIs from the FITS files. Use the event_cl files.

"""
from __future__ import division, print_function

def obsid_to_obj(obsid):
    """
    'Convert' the ObsID into the corresponding object name.

    obsid - Observation ID of the object of interest (10-digit str)
    """
    if type(obsid) != str:
        raise ValueError("Make sure the ObsID is given as a string!")
    cen_x3 = ['0034070101','0034070102','0034070103','0034070104',
              '1034070101','1034070102','1034070103','1034070104',
              '1034070105','1034070106'] #for Centaurus X-3!
    crab = ['1013010125'] #for Crab Pulsar!
    grs1915p105 = ['1103010117'] #for GRS1915+105
    gx349p2 = ['1034090101','1034090102','1034090103','1034090104','1034090105',
               '1034090106','1034090107','1034090108','1034090109','1034090110',
               '1034090111','1034090112']
    j0243_6124 = ['1050390'+str(i) for i in range(101,161)]
    at2018cow = ['12002501' + str(i+1).zfill(2) for i in range(26)]
    j0030_0451 = ['1060020'+str(i) for i in range(101,438)]
    j1231_1411 = ['00600601' + str(i+1).zfill(2) for i in range(13)] + ['1060060'+str(i) for i in range(101,374)] + ['2060060'+str(i) for i in range(301,374)]

    if obsid in cen_x3:
        return 'Cen X-3'
    elif obsid in crab:
        return 'Crab'
    elif obsid in grs1915p105:
        return 'GRS 1915+105'
    elif obsid in gx349p2:
        return 'GX 349+2'
    elif obsid in j0243_6124:
        return 'J0243.6+6124'
    elif obsid in at2018cow:
        return 'AT2018cow'
    elif obsid in j0030_0451:
        return 'J0030+0451'
    elif obsid in j1231_1411:
        return 'J1231-1411'
    else:
        raise ValueError("Object name not found. Make sure you've uploaded the ObsID + corresponding object name to the module Lv2_sources!")

if __name__ == "__main__":
    obsid_to_obj['1034070104']
