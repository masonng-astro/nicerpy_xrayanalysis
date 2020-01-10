#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:30am 2019

Program for gunzip - to unzip the files that I got from ciri!

"""
from __future__ import division, print_function
import numpy as np
import Lv0_dirs
import subprocess

def unzip_all(obsdir):
    subprocess.run(['gunzip','-r',obsdir])

if __name__ == "__main__":
    unzip_all(Lv0_dirs.NICER_DATADIR+'1034070101')
