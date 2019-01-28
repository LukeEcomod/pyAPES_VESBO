# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:59:47 2019

@author: L1656
"""
from tools.plotting import plot_lad_profiles
from matplotlib import pyplot as plt


# %% Plot lad profile
plot_lad_profiles("Krycklan_C2.txt",normed=True)
plt.xlim([0.0,0.15])
plt.ylim([0.0,25])

plt.savefig('figures/lad_profiles.png',dpi=500, transparent=False)

# %% Plot ET component comparison
from tools.iotools import read_forcing, xarray_to_df
from canopy.constants import MOLAR_MASS_H2O
import numpy as np
import pandas as pd
from tools.plotting import plot_xy, plot_timeseries_df
import seaborn as sns
pal = sns.color_palette("hls", 6)

def plot_ETcomponents(results):

    Data = read_forcing("Svarberget_EC_2014_2016_chi.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    results['overstory_transpiration'] = results['canopy_pt_transpiration'][:,0,:3].sum(dim='planttype')
    results['overstory_evaporation'] = results['canopy_evaporation_ml'][:,0,2:].sum(dim='canopy')
    results['understory_transpiration'] = results['canopy_pt_transpiration'][:,0,3:].sum(dim='planttype')
    results['understory_evaporation'] = results['canopy_evaporation_ml'][:,0,:2].sum(dim='canopy')
    results['overstory_throughfall'] = results['canopy_throughfall_ml'][:,0,2]

    variables=['forcing_precipitation',
               'overstory_transpiration',
               'understory_transpiration',
               'overstory_evaporation',
               'understory_evaporation',
               'ffloor_evaporation']

    series = []
    for var in variables:
        series.append(results[var].to_pandas())
    df = pd.concat(series, axis=1)
    df.columns = variables

    Data = Data.merge(df, how='outer', left_index=True, right_index=True)

    for var in variables:
        Data[var] *= 1800*1000  # m/s --> mm/30min

    Data = Data.resample('D').sum()

    Data['ET_mod'] = (Data.overstory_transpiration +
                      Data.understory_transpiration +
                      Data.overstory_evaporation +
                      Data.understory_evaporation +
                      Data.ffloor_evaporation)

    Data['E_mod'] = (Data.overstory_evaporation +
                     Data.understory_transpiration +
                     Data.understory_evaporation +
                     Data.ffloor_evaporation)

    Data['I_mod'] = (Data.overstory_evaporation)

    labels=['Modelled', 'Measured']

    plt.figure(figsize=(8,6))

    # Overstory transpiration
    ay = plt.subplot(345)
    plot_xy(Data.overstory_transpiration, Data.overstory_transpiration, color=pal[1], axislabels={'x': '', 'y': 'Modelled'})

    ax = plt.subplot(3,4,(6,8), sharey=ay)
    plot_timeseries_df(Data, ['overstory_transpiration'], colors=[pal[1],'k'], xticks=False,
                       labels=labels)
    plt.title('Overstory transpiration, T (mm/d)', fontsize=10)
    plt.legend(loc="upper rigth", frameon=False, borderpad=0.0)

    # Residual ET - Tr
    plt.subplot(349, sharex=ay, sharey=ay)
    plot_xy(Data.E_mod, Data.E_mod, color=pal[2], axislabels={'x': 'Measured', 'y': 'Modelled'})

    plt.subplot(3,4,(10,12), sharex=ax, sharey=ay)
    plot_timeseries_df(Data, ['E_mod','I_mod'], colors=[pal[2],'k'], xticks=True,
                       labels=labels)
    plt.title('Residual, ET - T (mm/d)', fontsize=10)
    plt.legend(loc="upper rigth", frameon=False, borderpad=0.0)

    # Evapotranspiration
    ay = plt.subplot(341, sharex=ay, sharey=ay)
    plot_xy(Data.ET_mod, Data.ET_mod, color=pal[0], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(3,4,(2,4), sharex=ax, sharey=ay)
    plot_timeseries_df(Data, ['ET_mod'], colors=[pal[0],'k'], xticks=False,
                       labels=labels)
    plt.title('Evapotranspiration, ET (mm/d)', fontsize=10)
    plt.legend(loc="upper rigth", frameon=False, borderpad=0.0)

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0)
