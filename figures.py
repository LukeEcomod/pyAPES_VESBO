# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:59:47 2019

@author: L1656
"""

from matplotlib import pyplot as plt
import matplotlib.dates

# %% Plot lad profile
from tools.plotting import plot_lad_profiles

def plot_lad():
    plot_lad_profiles("Krycklan_C2.txt",normed=True)
    plt.xlim([0.0,0.15])
    plt.ylim([0.0,25])
    
    plt.savefig('figures/lad_profiles.png',dpi=500, transparent=False)

# %% Plot ET component comparison
from tools.iotools import read_forcing
from canopy.constants import MOLAR_MASS_H2O
import numpy as np
import pandas as pd
from tools.plotting import plot_xy, plot_timeseries_df, plot_diurnal
import seaborn as sns
pal = sns.color_palette("hls", 6)

def plot_ETcomponents(results):

    Data = read_forcing("Partitioning_C2_2016-2017.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    results['overstory_transpiration'] = results['canopy_pt_transpiration'][:,0,:3].sum(dim='planttype')
    results['overstory_evaporation'] = results['canopy_evaporation_ml'][:,0,2:].sum(dim='canopy')
    results['understory_transpiration'] = results['canopy_pt_transpiration'][:,0,3:].sum(dim='planttype')
    results['understory_evaporation'] = results['canopy_evaporation_ml'][:,0,:2].sum(dim='canopy')
    results['overstory_throughfall'] = results['canopy_throughfall_ml'][:,0,2]

    variables=['forcing_precipitation',
               'overstory_throughfall',
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

    Data['I_mod'] = (Data.forcing_precipitation -
                     Data.overstory_throughfall)
    Data['TF_mod'] = (Data.overstory_throughfall)

    plt.figure(figsize=(9,8))

    # Overstory transpiration
    ay = plt.subplot(445)
    plot_xy(Data.Tr, Data.overstory_transpiration, color=pal[1], axislabels={'x': '', 'y': 'Modelled'})

    ax = plt.subplot(4,4,(6,8), sharey=ay)
    plot_timeseries_df(Data, ['Tr', 'overstory_transpiration'], colors=['k', pal[1]], xticks=False,
                       labels=['Measured (sapflow)', 'Modelled'])
    plt.title('Overstory transpiration, T (mm/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Residual ET - Tr
    plt.subplot(449, sharex=ay, sharey=ay)
    plot_xy(Data.E, Data.E_mod, color=pal[3], axislabels={'x': 'Measured', 'y': 'Modelled'})

    plt.subplot(4,4,(10,12), sharex=ax, sharey=ay)
    plot_timeseries_df(Data, ['E','E_mod'], colors=['k',pal[3]], xticks=False,
                       labels=['Residual from measured ET and T', 'Modelled E'])
    plt.title('Evaporation, E (mm/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Interception evaporation
    plt.subplot(4,4,(14,16), sharex=ax, sharey=ay)
    plot_timeseries_df(Data, ['overstory_evaporation'], colors=[pal[4]], xticks=True,
                       labels=['Modelled'])
    plt.title('Interception evaporation, I (mm/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Evapotranspiration
    ay = plt.subplot(441, sharex=ay, sharey=ay)
    plot_xy(Data.ET, Data.ET_mod, color=pal[2], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(4,4,(2,4), sharex=ax, sharey=ay)
    plot_timeseries_df(Data, ['ET','ET_mod'], colors=['k', pal[2]], xticks=False,
                       labels=['Measured (Eddy-covariance)', 'Modelled'])
    plt.title('Evapotranspiration, ET (mm/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0,w_pad=0.1, h_pad=0.1)
#
#    # Prec
#
#    plt.figure(figsize=(9,4))
#    plt.subplot(241)
#    plot_xy(Data.I, Data.I_mod, color=pal[3], axislabels={'x': 'Measured', 'y': 'Modelled'})
#
#    ax=plt.subplot(2,4,(2,4))
#    plot_timeseries_df(Data, ['I_mod','I'], colors=[pal[3],'k'], xticks=False,
#                       labels=labels)
#    plt.title('Interception (mm/d)', fontsize=10)
#    plt.legend(loc="upper right", frameon=False, borderpad=0.0)
#
#    plt.subplot(245)
#    plot_xy(Data.TF, Data.TF_mod, color=pal[3], axislabels={'x': 'Measured', 'y': 'Modelled'})
#
#    plt.subplot(2,4,(6,8), sharex=ax)
#    plot_timeseries_df(Data, ['TF_mod','TF','Prec'], colors=[pal[3],'k',pal[0]], xticks=True,
#                       labels=['Modelled', 'Measured', 'Prec'])
#    plt.title('Throughfall (mm/d)', fontsize=10)
#    plt.legend(loc="upper right", frameon=False, borderpad=0.0)
#
#    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
#    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))
#    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0)

# %% Plot ET component comparison
from tools.iotools import read_forcing
from canopy.constants import MOLAR_MASS_H2O
import numpy as np
import pandas as pd
from tools.plotting import plot_xy, plot_timeseries_df, plot_diurnal
import seaborn as sns
pal = sns.color_palette("hls", 6)

def plot_interception(results):

    grey=(0.4, 0.4, 0.4)

    Data = read_forcing("Interception_rain_events_2016.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    results['overstory_throughfall'] = results['canopy_throughfall_ml'][:,0,2]

    variables=['forcing_precipitation',
               'overstory_throughfall']

    series = []
    for var in variables:
        series.append(results[var].to_pandas())
    df = pd.concat(series, axis=1)
    df.columns = variables

    Data = Data.merge(df, how='outer', left_index=True, right_index=True)

    for var in variables:
        Data[var] *= 1800*1000  # m/s --> mm/30min

    Data['I_mod'] = (Data.forcing_precipitation -
                     Data.overstory_throughfall)

    Data_events = Data.groupby(['start_date']).sum()
    Data_events = Data_events.drop(['0.1.1900 0:00'], axis=0)

    Data_events.index = pd.to_datetime(Data_events.index, dayfirst=True)
    Data_events.sort_index(inplace=True)

    Data2 = Data[['start_date']].copy()
    Data2 = Data2.merge(Data_events, how='outer', left_index=True, right_index=True)
    Data2 = Data2.fillna(0.0)

    data_err = [(Data_events.I_SE_max_mm - Data_events.I_av_mm).values,
                (Data_events.I_av_mm - Data_events.I_SE_min_mm).values]
    mod_err = [(Data_events.I_SE_max_mm - Data_events.I_av_mm).values,
                (Data_events.I_av_mm - Data_events.I_SE_min_mm).values]

    plt.figure(figsize=(8.5,5))
    plt.subplot(231)
    _, _, bars = plt.errorbar(Data_events.I_av_mm, Data_events.I_mod, xerr=data_err, yerr=mod_err, fmt='o',ecolor=grey, color=grey, alpha=0.5)
    [bar.set_alpha(0.2) for bar in bars]
    plot_xy(Data_events.I_av_mm, Data_events.I_mod, color=(1,1,1), axislabels={'x': 'Data derived IL (mm)', 'y': 'Modelled IL (mm)'})
    plt.xlim([0,7.5])
    plt.ylim([0,7.5])

    axx=plt.subplot(232)
    _, _, bars = plt.errorbar(Data_events.ICOS_precip, Data_events.I_av_mm, yerr=data_err, fmt='o',ecolor=pal[0], color=pal[0], alpha=.5)
    [bar.set_alpha(0.2) for bar in bars]
    plt.xlim([0,30])

    plt.ylabel('Data derived IL (mm)')
    plt.xlabel('Precipitation (mm)')

    plt.subplot(233, sharey=axx)
    _, _, bars = plt.errorbar(Data_events.forcing_precipitation, Data_events.I_mod, yerr=mod_err, fmt='o',ecolor=pal[4], color=pal[4], alpha=.5)
    [bar.set_alpha(0.2) for bar in bars]
    plt.xlim([0,30])
    plt.ylabel('Modelled IL (mm)')
    plt.xlabel('Precipitation (mm)')
    plt.ylim([0,7.5])

    ax = plt.subplot(2,3,(4,5))
    plt.fill_between(Data2.index, np.cumsum(Data2['I_SE_max_mm'].values), 
                np.cumsum(Data2['I_SE_min_mm'].values), facecolor=pal[0], alpha=0.2)
    plt.fill_between(Data2.index, 1.2*np.cumsum(Data2['I_mod'].values), 
                0.8*np.cumsum(Data2['I_mod'].values), facecolor=pal[4], alpha=0.2)
    plot_timeseries_df(Data2, ['ICOS_precip','I_av_mm',
                              'I_mod'], colors=['k',pal[0],pal[4]], xticks=True,
                       labels=['Precipitation', 'Data derived IL', 'Modelled IL'], cum=True,
                       unit_conversion={'unit':'mm', 'conversion':1/1800}, legend=False)
    plt.ylabel('(mm)')
    plt.ylim([0,230])

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0)
    plt.legend(bbox_to_anchor=(1.05,0.5), loc="center left", frameon=False, borderpad=0.0)

# %% Plot energy 
from tools.iotools import read_forcing, xarray_to_df
import numpy as np
import pandas as pd
from tools.plotting import plot_xy, plot_timeseries_df, plot_diurnal
import seaborn as sns
pal = sns.color_palette("hls", 6)
import matplotlib.dates

def plot_energy(results):

    Data = read_forcing("Svarberget_EC_2014_2016_chi.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    variables=['canopy_SH','canopy_LE','canopy_SWnet','canopy_LWnet']
    df = xarray_to_df(results, variables, sim_idx=0)
    Data = Data.merge(df, how='outer', left_index=True, right_index=True)

    ixLE = np.where(np.isfinite(Data.LE))[0]
    ixSH = np.where(np.isfinite(Data.SH))[0]
    ixLWnet = np.where(np.isfinite(Data.LWnet))[0]
    ixSWnet = np.where(np.isfinite(Data.SWnet))[0]
    labels=['Modelled', 'Measured']

    # Energy
    plt.figure(figsize=(10,8))
    plt.subplot(441)
    plot_xy(Data.SWnet[ixSWnet], Data.canopy_SWnet[ixSWnet], color=pal[0], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(445)
    plot_xy(Data.LWnet[ixLWnet], Data.canopy_LWnet[ixLWnet], color=pal[1], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(449)
    plot_xy(Data.SH[ixSH], Data.canopy_SH[ixSH], color=pal[3], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(4,4,13)
    plot_xy(Data.LE[ixLE], Data.canopy_LE[ixLE], color=pal[2], axislabels={'x': 'Measured', 'y': 'Modelled'})

    ax = plt.subplot(4,4,(2,3))
    plot_timeseries_df(Data, ['canopy_SWnet', 'SWnet'], colors=[pal[0],'k'], xticks=False,
                       labels=['Modelled', 'Measured'], marker=[None, '.'])
    plt.title('Net shortwave radiation [W m-2]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(4,4,(6,7), sharex=ax)
    plot_timeseries_df(Data, ['canopy_LWnet', 'LWnet'], colors=[pal[1],'k'], xticks=False,
                       labels=['Modelled', 'Measured'], marker=[None, '.'])
    plt.title('Net longwave radiation [W m-2]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(4,4,(10,11), sharex=ax)
    plot_timeseries_df(Data, ['canopy_SH','SH'], colors=[pal[3],'k'], xticks=False,
                       labels=labels, marker=[None, '.'])
    plt.title('Sensible heat flux [W m-2]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(4,4,(14,15), sharex=ax)
    plot_timeseries_df(Data, ['canopy_LE','LE'], colors=[pal[2],'k'], xticks=True,
                       labels=labels, marker=[None, '.'])
    plt.title('Latent heat flux [W m-2]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    ax =plt.subplot(444)
    plot_diurnal(Data.SWnet[ixSWnet], color='k', legend=False)
    plot_diurnal(Data.canopy_SWnet[ixSWnet], color=pal[0], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')

    ax =plt.subplot(448, sharex=ax)
    plot_diurnal(Data.LWnet[ixLWnet], color='k', legend=False)
    plot_diurnal(Data.canopy_LWnet[ixLWnet], color=pal[1], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')
    
    plt.subplot(4,4,12, sharex=ax)
    plot_diurnal(Data.SH[ixSH], color='k', legend=False)
    plot_diurnal(Data.canopy_SH[ixSH], color=pal[3], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')

    plt.subplot(4,4,16, sharex=ax)
    plot_diurnal(Data.LE[ixLE], color='k', legend=False)
    plot_diurnal(Data.canopy_LE[ixLE], color=pal[2], legend=False)

    plt.tight_layout(rect=(0, 0, 0.90, 1),w_pad=0.1)

