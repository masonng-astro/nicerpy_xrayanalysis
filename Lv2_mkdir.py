#!/usr/bin/env python
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

    if os.path.exists(dir):
        print('The path already exists!')
        return
    else:
        print('This directory did not exist - creating ' + dir + ' now!')
        os.makedirs(dir)

if __name__ == "__main__":
    print('hi')
    makedir('/Volumes/Samsung_T5/hahaha')
