#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 2 10:45am 2019

Generating averaged power spectra for W = 1 to W = 300 for not demodulated
J1231-1411 data!

"""
from __future__ import division, print_function
import numpy as np
import time
from tqdm import tqdm
import glob
import subprocess

import Lv0_dirs
import Lv2_average_ps_methods,Lv2_average_merge_ps_methods
import Lv3_detection_level

import matplotlib.pyplot as plt

"""
for W in range(1,301):

    print('W = ' + str(W))

    pyfile = 'Lv3_average_ps_main.py'
    pyfile_contents = open(pyfile,'r')
    contents = pyfile_contents.read().split('\n')
    pyfile_contents.close()

    newstring_W = "    W = " + str(W) + " #number of consecutive Fourier bins to average over"
    newstring_filename = "    pngname = '" + Lv0_dirs.NICERSOFT_DATADIR + 'merged_events/merged000005/W_dir/W_' + str(W).zfill(3) + ".png'"
    pyfile_contents = open('Lv3_average_ps_main.py','w')
    for j in range(len(contents)):
        if j != 51 and j != 123:
            pyfile_contents.write(contents[j] + '\n')
        if j == 51:
            pyfile_contents.write(newstring_W + '\n')
        if j == 123:
            pyfile_contents.write(newstring_filename + '\n')

    pyfile_contents.close()

    execfile("Lv3_average_ps_main.py")
"""

from scipy.optimize import curve_fit
def neg_sqrt(x,a,n):
    return a*x**(n)

stds = np.array([0.17151563,0.13250076,0.11201556,0.09895773,0.08971697,0.08257810,
                0.07702619,0.07231257,0.06853172,0.06528318,0.06225732,0.05970226,
                0.05745892,0.05547180,0.05373177,0.05203332,0.05061853,0.04919654,
                0.04793032,0.04674205,0.04573258,0.04471104,0.04380930,0.04283883,
                0.04198711,0.04135160,0.04051614,0.03974083,0.03923577,0.03858371,
                0.03797808,0.03738553,0.03672342,0.03609305,0.03566748,0.03519816,
                0.03472995,0.03436365,0.03391587,0.03346646,0.03316299,0.03278463,
                0.03240944,0.03212927,0.03173375,0.03134329,0.03099390,0.03064317,
                0.03035390,0.02997892,0.02994112,0.02945002,0.02934270,0.02907256,
                0.02881691,0.02857167,0.02828932,0.02804486,0.02786123,0.02764023,
                0.02733023,0.02713993,0.02689832,0.02676551,0.02647976,0.02632849,
                0.02609446,0.02599551,0.02584631,0.02556697,0.02529887,0.02517167,
                0.02501944,0.02495649,0.02467688,0.02468062,0.02440038,0.02431525,
                0.02403045,0.02387878,0.02370795,0.02368280,0.02347985,0.02331529,
                0.02332113,0.02317480,0.02305960,0.02309585,0.02284412,0.02265442,
                0.02261264,0.02252684,0.02238650,0.02238772,0.02206224,0.02200986,
                0.02188044,0.02183891,0.02172427,0.02149658])

N = np.arange(1,101)

relation = stds[0]/np.sqrt(N)

p,cov = curve_fit(neg_sqrt,N,stds,p0=[0.2,-0.5])
print(p[0],p[1],np.sqrt(cov))


plt.plot(N,stds,'rx-')
plt.plot(N,relation,'bx-')
plt.plot(N,neg_sqrt(N,p[0],p[1]),'kx-')
plt.xlabel('N',fontsize=12)
plt.ylabel('Standard Dev.',fontsize=12)
plt.legend(('Actual','Expected','curvefit'),loc='best')
plt.show()
