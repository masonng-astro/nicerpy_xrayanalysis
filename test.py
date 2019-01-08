#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Jan 5th

@author: masonng


"""
from __future__ import division, print_function
import datetime
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy import signal

timestart = time.time()

print(datetime.datetime.now())

class Complex:
     def __init__(self, realpart, imagpart):
         self.r = realpart
         self.i = imagpart

x = Complex(4.0,3.0)
print(x.r,x.i)

timeend = time.time()

print(timeend-timestart)
