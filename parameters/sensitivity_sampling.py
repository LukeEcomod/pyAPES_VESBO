# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 18:05:03 2018

@author: L1656
"""
import numpy as np
import itertools

# sensitivity sampling for LAI (pine, spruce, decid)
frac = np.linspace(0,1,5)
frac = frac.tolist()
combinations = list(itertools.product(*[frac,frac,frac]))
combinations = filter(lambda x: sum(x) == 1, combinations)

LAI = np.linspace(0,6,7)
LAI = [0, 1, 2, 3, 5, 6]
LAIcombinations = []
for lai in LAI:
    LAIcombinations += (lai * np.array(combinations)).tolist()