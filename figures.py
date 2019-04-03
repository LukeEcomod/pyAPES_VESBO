# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 15:26:24 2019

@author: L1656
"""

# %% Plot vegetation inventory
import numpy as np
import pandas as pd

vege=['khaki','lightgreen', 'limegreen', 'forestgreen']
ff=['seagreen','peru','saddlebrown','khaki','lightgreen', 'limegreen', 'forestgreen']

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
        plt.annotate("y = %.2fx\nR$^2$ = %.2f" % (model[0], R2), (0.35, 0.8), xycoords='axes fraction', ha='center', va='center', fontsize=9)
        limx = [0.0, 1.2*max(x[np.isfinite(x)])]
        limy = [0.0, 1.2*max(y[np.isfinite(y)])]
        plt.plot(limx, [model[0]*limx[0], model[0]*limx[1]], 'r', linewidth=1)
        plt.ylim(limy)
        plt.xlim(limx)
        plt.title(title, fontsize=10)
        plt.xlabel(axislabels['x'])
        plt.ylabel(axislabels['y'])

    fp = r'H:\Lettosuo\aluskasvillisuus\regression_data.txt'
    dat = pd.read_csv(fp)
    pos = (-0.2,1.05)
    plt.figure(figsize=(8,2.3))
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

    fp = r'H:\Lettosuo\aluskasvillisuus\inventory_data.txt'
    dat = pd.read_csv(fp,index_col=0)
    labels=['SP$_{\mathrm{all},2009}$', 'VP$_{\mathrm{clc},2015}$', 'VP$_{\mathrm{ref},2015}$', 'SP$_{\mathrm{ref},2017}$', 'VP$_{\mathrm{ref},2017}$', 'SP$_{\mathrm{par},2017}$', 'SP$_{\mathrm{par},2018}$', 'VP$_{\mathrm{clc},2017}$', 'VP$_{\mathrm{clc},2018}$']
    pos = (-0.16,1.02)
    width = 0.75
    plt.figure(figsize=(6.5,4))
    ax1=plt.subplot(2,1,1)
    plt.annotate('a', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.annotate('reference', (1.2,70))
    plt.annotate(' partial\nharvest', (4.9,65))
    plt.annotate('clear-cut', (6.9,70))
    dat[['seedlings','shrubs','graminoid','forbs']].plot(kind='bar', stacked=True, ax=ax1, colors=vege, width=width)
    plt.plot([4.5, 4.5],[0,1110],'--k')
    plt.plot([6.5, 6.5],[0,1110],'--k')
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.ylim([0,75])
    ax1.set_yticks([0, 25, 50, 75])
    plt.ylabel('coverage [%]',labelpad=10.5)
    ax1.legend().set_visible(False)
    ax1_1 = ax1.twinx()
    ax1_1.plot(range(9),dat['LAI'],'ok', markersize=4)
    ax1_1.plot(range(9),dat['LAI_meas'],'xk', markersize=4)
    plt.ylabel('LAI [m$^2$ m$^{-2}$]',labelpad=10)
    plt.ylim([0,2])
    ax1.legend().set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1_1.spines['top'].set_visible(False)

    ax2=plt.subplot(2,1,2)
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    dat[['moss','litter','baresoil','seedlings','shrubs','graminoid','forbs']].plot(kind='bar', stacked=True, ax=ax2, colors=ff, width=width)
    plt.plot([1, 2],[110,110],'ok', label='LAI estimated', markersize=4)
    plt.plot([1, 2],[110,110],'xk', label='LAI measured', markersize=4)
    handles, labels1 = ax2.get_legend_handles_labels()
    labels1[2:]=labels1[:1:-1]
    handles[2:]=handles[:1:-1]
    ax2.legend(handles, labels1,bbox_to_anchor=(1.02,0.35), loc="center left", frameon=False, borderpad=0.0)
    plt.plot([4.5, 4.5],[0,1110],'--k')
    plt.plot([6.5, 6.5],[0,1110],'--k')
    plt.ylabel('coverage [%]')
    plt.ylim([0,100])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_xticklabels(labels)
    plt.xticks(rotation=65)
    plt.tight_layout()
#    plt.savefig('figures/case_Lettosuo/vege_inventories.png',dpi=300, transparent=False)

# %% Plot lad profile
from pyAPES_utilities.plotting import plot_lad_profiles
import matplotlib.pyplot as plt

def plot_lad(biomass_function='marklund_mod'):
    pos = (-0.1,1.03)
    plt.figure(figsize=(8,4))
    plot_lad_profiles("letto2009.txt",subplot=1,subplots=3, biomass_function=biomass_function)
    plt.annotate('a', (-0.25,1.03), xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.xlim([0.0,0.6])
    plt.ylim([0.0,27])
    plot_lad_profiles("letto2014.txt",subplot=2,subplots=3, biomass_function=biomass_function)
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.xlim([0.0,0.6])
    plt.ylim([0.0,27])
    plt.ylabel('')
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plot_lad_profiles("letto2016_partial.txt",subplot=3,subplots=3, biomass_function=biomass_function)
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.xlim([0.0,0.6])
    plt.ylim([0.0,27])
    plt.ylabel('')
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plt.tight_layout()
#    plt.savefig('figures/case_Lettosuo/lad_profiles.png',dpi=300, transparent=False)

def plot_wtd(results):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_timeseries_xr

    # Read observed WTD
    WTD = read_forcing("Lettosuo_WTD_pred.csv", cols='all')

    plt.figure()
    plt.fill_between(WTD.index, WTD['control_max'].values, WTD['control_min'].values,
                     facecolor='k', alpha=0.3)
    plt.plot(WTD.index, WTD['control'].values,':k', linewidth=1.0)

    plt.fill_between(WTD.index, WTD['partial_max'].values, WTD['partial_min'].values,
                     facecolor='b', alpha=0.3)
    plt.plot(WTD.index, WTD['partial'].values,':b', linewidth=1.0)

    plt.fill_between(WTD.index, WTD['clearcut_max'].values, WTD['clearcut_min'].values,
                     facecolor='r', alpha=0.3)
    plt.plot(WTD.index, WTD['clearcut'].values,':r', linewidth=1.0)

    plot_timeseries_xr(results, 'soil_ground_water_level', colors=['k','b','r'], xticks=True)
