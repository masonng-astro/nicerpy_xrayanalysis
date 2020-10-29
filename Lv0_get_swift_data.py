#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues April 21 7:10pm 2020

Program that takes in a FITS table, from the HEASARC archive, of observation IDs
and pointing modes (for Swift), then uses the Perl script download_wget.pl to
obtain the data files!

"""
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
import Lv0_dirs,Lv1_data_gtis,Lv2_presto_subroutines,Lv2_mkdir
import os
from tqdm import tqdm
import subprocess
import pathlib

Lv0_dirs.global_par()

def download_txt(txtfile):
    """
    Given the text file of download instructions, take the URLs of where the data
    are stored within the HEASARC archive, and use download_wget.pl to retrieve them
    """
    contents = open(txtfile,'r').read().split('\n')
    urls = [contents[i].split(' ')[-1] for i in range(len(contents)-1)]

    for i in tqdm(range(len(urls))):
        subprocess.run(['perl','/Volumes/Samsung_T5/download_wget.pl',urls[i]])

if __name__ == "__main__":
    textfile = '/Volumes/Samsung_T5/NGC300_ULX_Swift/ngc300x1-ulx1_photon_pointing.txt'
    download_txt(textfile)
