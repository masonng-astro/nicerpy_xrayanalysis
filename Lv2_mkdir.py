#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thurs Jan 10 3:50pm 2019

Creating a folder if it does not exist in the directory.

https://gist.github.com/keithweaver/562d3caa8650eefe7f84fa074e9ca949
"""
from __future__ import division, print_function
import os

def makedir(dir):
    """ 
    Creating a folder if it does not exist in the directory.
 
    dir - desired directory (provide FULL path!)
    """

    try:
        if not os.path.exists(dir):
            os.makedirs(dir)
    except OSError:
        print('This directory did not exist - creating ' + dir + ' now!')
