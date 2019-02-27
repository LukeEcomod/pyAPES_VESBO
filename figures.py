# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 15:26:24 2019

@author: L1656
"""

# %% Plot vegetation inventory
import numpy as np
import pandas as pd

def plot_vegetation():

    def plot_regression(x, y, color='k', title='', axislabels={'x':'', 'y':''},alpha=0.5):
        """
        Plot x,y scatter with linear regression line.
        Args:
            x,y (array): arrays for x and y data
            color (str or tuple): color code
            title (str): title of plot
            axislabels (dict)
                'x' (str): x-axis label
                'y' (str): y-axis label
        """
        plt.scatter(x, y, marker='o', color=color, alpha=alpha)
        idx = np.isfinite(x) & np.isfinite(y)
        model, resid, _, _ = np.linalg.lstsq(x[idx,np.newaxis], y[idx,np.newaxis])
        R2 = 1 -resid / (y[idx,np.newaxis].size*y[idx,np.newaxis].var())
        plt.annotate("y = %.2fx \nR$^2$ = %.2f" % (model[0], R2), (0.45, 0.85), xycoords='axes fraction', ha='center', va='center', fontsize=9)
        limx = [0.0, 1.2*max(x[np.isfinite(x)])]
        limy = [0.0, 1.2*max(y[np.isfinite(y)])]
        plt.plot(limx, [model[0]*limx[0], model[0]*limx[1]], 'r', linewidth=1)
        plt.ylim(limy)
        plt.xlim(limx)
        plt.title(title)
        plt.xlabel(axislabels['x'])
        plt.ylabel(axislabels['y'])

    fp = r'H:\Lettosuo\aluskasvillisuus\regression_data.txt'
    dat = pd.read_csv(fp)
    pos = (-0.2,1.05)
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,4,1)
    plt.annotate('a', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_regression(dat['gra_cov'], dat['gra_bm'],
                    axislabels={'x':'coverage [%]', 'y':'biomass [g m$^{-2}$]'}, title='Graminoid')
    plt.subplot(1,4,2)
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_regression(dat['for_cov'], dat['for_bm'],
                    axislabels={'x':'coverage [%]', 'y':''}, title='Forbs')
    plt.subplot(1,4,3)
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_regression(dat['shr_cov'], dat['shr_bm'],
                    axislabels={'x':'coverage [%]', 'y':''}, title='Shrubs')
    plt.subplot(1,4,4)
    plt.annotate('d', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_regression(dat['bm'], dat['LAI'], alpha=0.5, title='All',
                    axislabels={'x':'biomass [g m$^{-2}$]', 'y':'LAI [m$^2$ m$^{-2}$]'})
    plt.tight_layout()
#    plt.savefig('figures/case_Lettosuo/vege_regressions.png',dpi=300, transparent=False)

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
