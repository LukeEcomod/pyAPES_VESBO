# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 09:57:38 2020

@author: 03081268
"""


import xarray
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

fpath = r'results/Scenarios/S1/*.nc'
files = glob.glob(fpath)

data = []
for k in range(0,5):
    data.append(xarray.open_dataset(files[k]))
print('done')