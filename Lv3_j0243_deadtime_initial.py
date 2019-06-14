#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 8 10:32pm 2019

Looking for the ObsID with minimal deadtime corrections.

Note: Did this script in a rush. Make sure to rename!

"""
from __future__ import division, print_function
import numpy as np

deadtime = np.array([2643.4, 2782.6, 16645.4, 19320.5, 2198.3, 0, 420.9,
                    615.7, 0, 735.4, 0, 495.8, 3126.0, 257.4,
                    1276.5, 1915.8, 0, 374.9, 2171.9, 2240.9,
                    1517.3, 3676.9, 1016.0, 1590.0, 1854.5,
                    109.0, 640.9, 263.5, 594.0, 279.7, 618.5,
                    208.1, 645.3, 3193.7, 2029.3, 2686.1, 3368.9,
                    3207.7, 2127.0, 2104.8, 3492.7, 2755.4, 1564.1,
                    1805.5, 6964.3, 1442.7, 1077.7, 1649.0, 3966.1,
                    1917.4, 1699.0, 1111.4, 370.4, 16463.2, 12004.6,
                    6281.5, 5513.6, 2621.7, 459.3, 149.1])

on_source_t = np.array([3099, 4171, 21590, 24029, 2497, 41, 1495,
                        1394, 2334, 900, 689, 1003, 4566, 311, 1495,
                        3816, 18, 405, 3917, 3916, 1610, 3756, 1291,
                        1704, 1887, 210, 733, 482, 1034, 528, 718,
                        210, 990, 5319, 3133, 3288, 4541, 3230,
                        2266, 2114, 3636, 2768, 2326, 2175, 7058,
                        1535, 1251, 1705, 5003, 2025, 2278, 1606,
                        500, 19208, 14521, 7735, 5955, 2785, 492,
                        643])

percent = deadtime/on_source_t*100

for i in range(len(percent)):
    print(str(i+1),on_source_t[i],percent[i])