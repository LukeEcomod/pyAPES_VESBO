# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:59:47 2019

@author: L1656
"""
from tools.plotting import plot_lad_profiles
from matplotlib import pyplot as plt

plot_lad_profiles("Krycklan_C2.txt",normed=True)
plt.xlim([0.0,0.15])
plt.ylim([0.0,25])

plt.savefig('figures/lad_profiles.png',dpi=500, transparent=False)
