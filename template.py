#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 7 12:02pm 2019


"""
from __future__ import division
import datetime
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from scipy import stats
from scipy import signal

timestart = time.time()

print(datetime.datetime.now())

timeend = time.time()

print (timeend-timestart)
