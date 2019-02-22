# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 15:26:24 2019

@author: L1656
"""
# %% Plot vegetation inventory



# %% Plot lad profile
from pyAPES_utilities.plotting import plot_lad_profiles
import matplotlib.pyplot as plt

def plot_lad():
    pos = (-0.1,1.03)
    plt.figure(figsize=(8,4))
    plot_lad_profiles("letto2009.txt",subplot=1,subplots=3)
    plt.annotate('a', (-0.25,1.03), xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.xlim([0.0,0.75])
    plt.ylim([0.0,25])
    plot_lad_profiles("letto2014.txt",subplot=2,subplots=3)
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.xlim([0.0,0.75])
    plt.ylim([0.0,25])
    plt.ylabel('')
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plot_lad_profiles("letto2016_partial.txt",subplot=3,subplots=3)
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.xlim([0.0,0.75])
    plt.ylim([0.0,25])
    plt.ylabel('')
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plt.tight_layout()
#    plt.savefig('figures/case_Lettosuo/lad_profiles.png',dpi=300, transparent=False)
