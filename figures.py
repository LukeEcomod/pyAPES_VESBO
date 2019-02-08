# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:59:47 2019

@author: L1656
"""

from matplotlib import pyplot as plt
import matplotlib.dates

# %% Plot lad profile
from tools.plotting import plot_lad_profiles

def plot_lad(file="Krycklan_C2.txt"):
    plot_lad_profiles(file,normed=True)
    plt.xlim([0.0,0.15])
    plt.ylim([0.0,25])
    
    plt.savefig('figures/lad_profiles.png',dpi=500, transparent=False)

# %% Plot ET component comparison
from tools.iotools import read_forcing, xarray_to_df
from canopy.constants import MOLAR_MASS_H2O
import numpy as np
import pandas as pd
from tools.plotting import plot_xy, plot_timeseries_df, plot_diurnal
import seaborn as sns
pal = sns.color_palette("hls", 6)

def plot_ETcomponents(results):

    Data = read_forcing("Sapflow_2016-2017.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    Data2 = read_forcing("Svarberget_EC_2014_2016_chi.csv", cols=['ET_gapfilled'],
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data2.ET_gapfilled = Data2.ET_gapfilled * 1e-3 * MOLAR_MASS_H2O * 3600 * 24  # mmol m-2 s-1 - > mm/d

    results['overstory_transpiration'] = results['canopy_pt_transpiration'][:,:,:3].sum(dim='planttype')
    results['overstory_evaporation'] = results['canopy_evaporation_ml'][:,:,2:].sum(dim='canopy')
    results['understory_transpiration'] = results['canopy_pt_transpiration'][:,:,3:].sum(dim='planttype')
    results['understory_evaporation'] = results['canopy_evaporation_ml'][:,:,:2].sum(dim='canopy')

    results['ET_mod'] = (results['overstory_transpiration'] +
                         results['understory_transpiration'] +
                         results['overstory_evaporation'] +
                         results['understory_evaporation'] +
                         results['ffloor_evaporation'])

    results['T_mod'] = results['overstory_transpiration']

    results['E_mod'] = (results['understory_transpiration'] +
                        results['overstory_evaporation'] +
                        results['understory_evaporation'] +
                        results['ffloor_evaporation'])

    results['I_mod'] = results['overstory_evaporation']

    results['Precip'] = results['forcing_precipitation']

    variables = ['ET_mod', 'T_mod', 'E_mod','I_mod','Precip']

    series = []
    varnames = []
    for var in variables:
        series.append(results.sel(simulation=0)[var].to_pandas())
        varnames.append(var)
        series.append(results[var].min(dim='simulation').to_pandas())
        varnames.append(var + '_min')
        series.append(results[var].max(dim='simulation').to_pandas())
        varnames.append(var + '_max')
    df = pd.concat(series, axis=1)
    df.columns = varnames

    Data = Data.merge(df, how='outer', left_index=True, right_index=True)
    Data = Data.merge(Data2, how='outer', left_index=True, right_index=True)
    Data = Data.fillna(method='ffill')

    for var in varnames:
        Data[var] *= 1000*3600*24  # m/s --> mm/d

    Data = Data.resample('D').mean()

    Data['ET'] = Data.ET_gapfilled
    Data['E'] = Data.ET - Data.Transp_av
    Data['E_min'] = Data.ET - Data.Transp_max
    Data['E_max'] = Data.ET - Data.Transp_min

    plt.figure(figsize=(9,8))

    # Overstory transpiration
    ay = plt.subplot(445)
    plot_xy(Data.Transp_av, Data.T_mod, color=pal[0], axislabels={'x': '', 'y': 'Modelled'})

    ax = plt.subplot(4,4,(6,8), sharey=ay)
    plt.fill_between(Data.index, Data['T_mod_max'].values, Data['T_mod_min'].values,
                     facecolor=pal[0], alpha=0.3)
    plt.fill_between(Data.index, Data['Transp_max'].values, Data['Transp_min'].values,
                     facecolor='k', alpha=0.2)
    plot_timeseries_df(Data, ['T_mod','Transp_av'], colors=[pal[0],'k'], xticks=False,
                       labels=['Modelled','Measured (sapflow)'], linestyles=['-',':'])
    plt.title('Overstory transpiration, T (mm/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Residual ET - Tr
    plt.subplot(449, sharex=ay, sharey=ay)
    plot_xy(Data.E, Data.E_mod, color=pal[3], axislabels={'x': 'Measured', 'y': 'Modelled'})

    plt.subplot(4,4,(10,12), sharex=ax, sharey=ay)
    plt.fill_between(Data.index, Data['E_mod_max'].values, Data['E_mod_min'].values,
                     facecolor=pal[3], alpha=0.3)
    plt.fill_between(Data.index, Data['E_max'].values, Data['E_min'].values,
                     facecolor='k', alpha=0.2)
    plot_timeseries_df(Data, ['E_mod','E'], colors=[pal[3],'k'], xticks=False,
                       labels=['Modelled E', 'Residual from measured ET and T'], linestyles=['-',':'])
    plt.title('Evaporation, E (mm/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Interception evaporation
    plt.subplot(4,4,(14,16), sharex=ax, sharey=ay)
    plt.fill_between(Data.index, Data['I_mod_max'].values, Data['I_mod_min'].values,
                     facecolor=pal[4], alpha=0.3)
    plot_timeseries_df(Data, ['I_mod'], colors=[pal[4]], xticks=True,
                       labels=['Modelled'])
    plt.title('Interception evaporation, I (mm/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Evapotranspiration
    plt.subplot(441, sharex=ay, sharey=ay)
    plot_xy(Data.ET, Data.ET_mod, color=pal[2], axislabels={'x': '', 'y': 'Modelled'})
    plt.xlim([0, 5])

    plt.subplot(4,4,(2,4), sharex=ax, sharey=ay)
    plt.fill_between(Data.index, Data['ET_mod_max'].values, Data['ET_mod_min'].values,
                     facecolor=pal[3], alpha=0.3)
    plot_timeseries_df(Data, ['ET_mod','ET'], colors=[pal[2],'k'], xticks=False,
                       labels=['Modelled','Measured (Eddy-covariance)'], linestyles=['-',':'])
    plt.ylim([0, 5])
    plt.title('Evapotranspiration, ET (mm/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0,w_pad=0.1, h_pad=0.1)

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

    index=pd.date_range('07-01-2016','11-01-2016',freq='0.5H')
    data=pd.DataFrame(index=index, columns=[])

    fp = r'C:\Users\L1656\Documents\Git_repos\Modeling_cases\pyAPES_Krycklan_C2\forcing\Interception_events_2016.csv'
    dat = pd.read_csv(fp, sep=';', header='infer')
    dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
    data=data.merge(dat, how='outer', left_index=True, right_index=True)
    data.ix[:,1:] = data.ix[:,1:].fillna(0.0)
    data.ix[:,0] = data.ix[:,0].fillna(method='bfill')
    data.ix[:,0] = data.ix[:,0].fillna('01-01-70 00:00')

    results['overstory_throughfall'] = results['canopy_throughfall_ml'][:,:,2]

    results['I_mod'] = results['forcing_precipitation'] - results['overstory_throughfall']

    results['Precip'] = results['forcing_precipitation']

    variables = ['I_mod','Precip']

    series = []
    varnames = []
    for var in variables:
        series.append(results.sel(simulation=0)[var].to_pandas())
        varnames.append(var)
        series.append(results[var].min(dim='simulation').to_pandas())
        varnames.append(var + '_min')
        series.append(results[var].max(dim='simulation').to_pandas())
        varnames.append(var + '_max')
    df = pd.concat(series, axis=1)
    df.columns = varnames

    data = data.merge(df, how='outer', left_index=True, right_index=True)

    for var in varnames:
       data[var] *= 1800*1000  # m/s --> mm/30min

    Data_events = data.groupby(['Sampling']).sum()

    Data_events.index = pd.to_datetime(Data_events.index, dayfirst=True)
    Data_events.sort_index(inplace=True)
    Data_events = Data_events.drop(Data_events.index[0], axis=0)

    Data_events['I_av'] = Data_events.Interc_wav.values
    Data_events['I_max'] = Data_events.Interc_wav.values + 0.5 * Data_events.Interc_stdev.values
    Data_events['I_min'] = Data_events.Interc_wav.values - 0.5 * Data_events.Interc_stdev.values

    Data2 = data[['Sampling']].copy()
    Data2 = Data2.merge(Data_events, how='outer', left_index=True, right_index=True)
    Data2 = Data2.fillna(0.0)

    data_err = [0.5 * Data_events.Interc_stdev.values,
                0.5 * Data_events.Interc_stdev.values]
    mod_err = [(Data_events.I_mod - Data_events.I_mod_min).values,
               (Data_events.I_mod_max - Data_events.I_mod).values]

    plt.figure(figsize=(8.5,5))
    plt.subplot(231)
    _, _, bars = plt.errorbar(Data_events.I_av, Data_events.I_mod, xerr=data_err, yerr=mod_err, fmt='o',ecolor=grey, color=grey, alpha=0.5)
    [bar.set_alpha(0.2) for bar in bars]
    plot_xy(Data_events.I_av, Data_events.I_mod, color=(1,1,1), axislabels={'x': 'Data derived IL (mm)', 'y': 'Modelled IL (mm)'})
    plt.xlim([0,7.5])
    plt.ylim([0,7.5])

    axx=plt.subplot(232)
    _, _, bars = plt.errorbar(Data_events.precip_open_area_mm, Data_events.I_av, yerr=data_err, fmt='o',ecolor=pal[0], color=pal[0], alpha=.5)
    [bar.set_alpha(0.2) for bar in bars]
    plt.xlim([0,30])

    plt.ylabel('Data derived IL (mm)')
    plt.xlabel('Precipitation* (mm)')

    plt.subplot(233, sharey=axx)
    _, _, bars = plt.errorbar(Data_events.Precip, Data_events.I_mod, yerr=mod_err, fmt='o',ecolor=pal[4], color=pal[4], alpha=.5)
    [bar.set_alpha(0.2) for bar in bars]
    plt.xlim([0,30])
    plt.ylabel('Modelled IL (mm)')
    plt.xlabel('Precipitation (mm)')
    plt.ylim([0,7.5])

    ax = plt.subplot(2,3,(4,5))
    plt.fill_between(Data2.index, np.cumsum(Data2['I_max'].values), 
                np.cumsum(Data2['I_min'].values), facecolor=pal[0], alpha=0.2)
    plt.fill_between(Data2.index, np.cumsum(Data2['I_mod_max'].values), 
                np.cumsum(Data2['I_mod_min'].values), facecolor=pal[4], alpha=0.2)
    plot_timeseries_df(Data2, ['precip_open_area_mm','Precip','I_av','I_mod'],
                       colors=['k',grey,pal[0],pal[4]], xticks=True,
                       labels=['Precipitation*', 'Precipitation', 'Data derived IL', 'Modelled IL'], cum=True,
                       unit_conversion={'unit':'mm', 'conversion':1/1800}, legend=False)
    plt.ylabel('(mm)')
    plt.ylim([0,250])
    plt.xlim(['7.1.2016', '11.1.2016'])

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
    plot_xy(Data.SWnet[ixSWnet], Data.canopy_SWnet[ixSWnet], color=pal[0], alpha=0.2, axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(445)
    plot_xy(Data.LWnet[ixLWnet], Data.canopy_LWnet[ixLWnet], color=pal[1], alpha=0.2, axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(449)
    plot_xy(Data.SH[ixSH], Data.canopy_SH[ixSH], color=pal[3], alpha=0.2, axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(4,4,13)
    plot_xy(Data.LE[ixLE], Data.canopy_LE[ixLE], color=pal[2], alpha=0.2, axislabels={'x': 'Measured', 'y': 'Modelled'})

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
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m'))

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

# %% Plot forcing

def plot_forcing(results):

    Dat = read_forcing("Svarberget_forcing_2014_2016.csv",cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    # VPD
    from canopy.micromet import e_sat
    # vapor pressure
    esat, s = e_sat(Dat['Tair'].values)
    Dat['VPD'] = 1e-3 * (esat - Dat['H2O'].values * Dat['P'].values)
    Dat['Prec'] = Dat['Prec'].values * 1000 * 3600 * 24
    Dat['Wliq'] = Dat['Wliq'].values * 100
    Dat['Par'] = Dat['diffPar'].values + Dat['dirPar'].values
    Data = Dat.resample('D').mean()
    Data['Tmax'] = Dat['Tair'].resample('D').max()
    Data['Tmin'] = Dat['Tair'].resample('D').min()
    Data['VPDmax'] = Dat['VPD'].resample('D').max()
    Data['VPDmin'] = Dat['VPD'].resample('D').min()
    plt.figure(figsize=(7,8))

    # Precipitation
    ax = plt.subplot(411)
    plot_timeseries_df(Data, ['Prec','Wliq'], colors=[pal[4],pal[3]], xticks=False,legend=False,
                       labels=['Preciptation', 'Top soil moisture'])
    plt.title('Precipitation (mm/d) and soil moisture (%)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)
    plt.ylim([0, 35])

    # Temperature
    plt.subplot(412, sharex=ax)
    plt.fill_between(Data.index, Data['Tmax'].values, Data['Tmin'].values, facecolor=pal[0], alpha=0.2)
    plot_timeseries_df(Data, ['Tair','Tsoil'], colors=pal[0:], xticks=False,legend=False,
                       labels=['Air', 'Top soil'])
    plt.title('Temperature (degC)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)
    plt.ylim([-5, 25])

    # VPD
    plt.subplot(413, sharex=ax)
    plt.fill_between(Data.index, Data['VPDmax'].values, Data['VPDmin'].values, facecolor=pal[2], alpha=0.3)
    plot_timeseries_df(Data, ['VPD'], colors=[pal[2]], xticks=False,legend=False)
    plt.title('Vapor pressure deficit (kPa)', fontsize=10)
    plt.ylim([0, 2.0])

    # Par
    plt.subplot(414, sharex=ax)
    plot_timeseries_df(Data, ['Par', 'diffPar'], colors=pal[0:], xticks=True,legend=False,
                       labels=['Total', 'Diffuse'])
    plt.title('Photosynthetically active radiation (W/m2)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)
    plt.ylim([0, 170])

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0)

def plot_interception2(results):

    grey=(0.4, 0.4, 0.4)

    Data = read_forcing("Interception_rain_events_2016.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    results['overstory_throughfall'] = results['canopy_throughfall_ml'][:,:,2]

    results['I_mod'] = results['forcing_precipitation'] - results['overstory_throughfall']

    results['Precip'] = results['forcing_precipitation']

    variables = ['I_mod','Precip']

    series = []
    varnames = []
    for var in variables:
        series.append(results.sel(simulation=0)[var].to_pandas())
        varnames.append(var)
        series.append(results[var].min(dim='simulation').to_pandas())
        varnames.append(var + '_min')
        series.append(results[var].max(dim='simulation').to_pandas())
        varnames.append(var + '_max')
    df = pd.concat(series, axis=1)
    df.columns = varnames

    Data = Data.merge(df, how='outer', left_index=True, right_index=True)

    for var in varnames:
        Data[var] *= 1800*1000  # m/s --> mm/30min

    Data_events = Data.groupby(['start_date']).sum()
    Data_events = Data_events.drop(['0.1.1900 0:00'], axis=0)

    Data_events.index = pd.to_datetime(Data_events.index, dayfirst=True)
    Data_events.sort_index(inplace=True)

    Data2 = Data[['start_date']].copy()
    Data2 = Data2.merge(Data_events, how='outer', left_index=True, right_index=True)
    Data2 = Data2.fillna(0.0)

    data_err = [(Data_events.I_av_mm - Data_events.I_SE_min_mm).values,
                (Data_events.I_SE_max_mm - Data_events.I_av_mm).values]
    mod_err = [(Data_events.I_mod - Data_events.I_mod_min).values,
               (Data_events.I_mod_max - Data_events.I_mod).values]

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
    _, _, bars = plt.errorbar(Data_events.Precip, Data_events.I_mod, yerr=mod_err, fmt='o',ecolor=pal[4], color=pal[4], alpha=.5)
    [bar.set_alpha(0.2) for bar in bars]
    plt.xlim([0,30])
    plt.ylabel('Modelled IL (mm)')
    plt.xlabel('Precipitation (mm)')
    plt.ylim([0,7.5])

    ax = plt.subplot(2,3,(4,5))
    plt.fill_between(Data2.index, np.cumsum(Data2['I_SE_max_mm'].values), 
                np.cumsum(Data2['I_SE_min_mm'].values), facecolor=pal[0], alpha=0.2)
    plt.fill_between(Data2.index, np.cumsum(Data2['I_mod_max'].values), 
                np.cumsum(Data2['I_mod_min'].values), facecolor=pal[4], alpha=0.2)
    plot_timeseries_df(Data2, ['ICOS_precip','I_av_mm','I_mod'],
                       colors=['k',pal[0],pal[4]], xticks=True,
                       labels=['Precipitation', 'Data derived IL', 'Modelled IL'], cum=True,
                       unit_conversion={'unit':'mm', 'conversion':1/1800}, legend=False)
    plt.ylabel('(mm)')
    plt.ylim([0,230])

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0)
    plt.legend(bbox_to_anchor=(1.05,0.5), loc="center left", frameon=False, borderpad=0.0)

def plot_GPP(results):

    Data = read_forcing("Svarberget_EC_2014_2016_chi.csv", cols=['GPP_gapfilled'],
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data.GPP_gapfilled = Data.GPP_gapfilled * 1e-3 * 3600 * 24 * 12.01e-3 # umol m-2 s-1 - > gC m-2 d-1

    variables = ['canopy_GPP']

    series = []
    varnames = []
    for var in variables:
        series.append(results.sel(simulation=0)[var].to_pandas())
        varnames.append(var)
        series.append(results[var].min(dim='simulation').to_pandas())
        varnames.append(var + '_min')
        series.append(results[var].max(dim='simulation').to_pandas())
        varnames.append(var + '_max')
    df = pd.concat(series, axis=1)
    df.columns = varnames

    Data = Data.merge(df, how='outer', left_index=True, right_index=True)
    Data = Data.fillna(method='ffill')

    for var in varnames:
        Data[var] *= 1e-3 * 3600 * 24 * 12.01e-3 # umol m-2 s-1 - > gC m-2 d-1

    Data = Data.resample('D').mean()

    plt.figure(figsize=(9,2.5))

    ax=plt.subplot(141)
    plot_xy(Data.GPP_gapfilled, Data.canopy_GPP, color=pal[5], axislabels={'x': 'Measured', 'y': 'Modelled'})
    plt.xlim([0,15])

    ay=plt.subplot(1,4,(2,4), sharey=ax)
    plt.fill_between(Data.index, Data['canopy_GPP_max'].values, Data['canopy_GPP_min'].values,
                     facecolor=pal[5], alpha=0.3)
    plot_timeseries_df(Data, ['canopy_GPP','GPP_gapfilled'], colors=[pal[5],'k'], xticks=True,
                       labels=['Modelled','Measured (Eddy-covariance)'], linestyles=['-',':'])
    plt.ylim([0,15])
    plt.title('Gross primary production (g C/m2/d)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    ay.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ay.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0)

def plot_Tcomponents(results):

    Data = read_forcing("Sapflow_2016-2017.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    results['pine'] = results['canopy_pt_transpiration'][:,:,0]
    results['spruce'] = results['canopy_pt_transpiration'][:,:,1]
    results['birch'] = results['canopy_pt_transpiration'][:,:,2]
    results['shrubs'] = results['canopy_pt_transpiration'][:,:,3]

    variables = ['pine', 'spruce', 'birch','shrubs']

    series = []
    varnames = []
    for var in variables:
        series.append(results.sel(simulation=0)[var].to_pandas())
        varnames.append(var)
        series.append(results[var].min(dim='simulation').to_pandas())
        varnames.append(var + '_min')
        series.append(results[var].max(dim='simulation').to_pandas())
        varnames.append(var + '_max')
    df = pd.concat(series, axis=1)
    df.columns = varnames

    Data = Data.merge(df, how='outer', left_index=True, right_index=True)
    Data = Data.fillna(method='ffill')

    for var in varnames:
        Data[var] *= 1000*3600*24  # m/s --> mm/d

    Data = Data.resample('D').mean()

    plt.figure(figsize=(9,8))

    # Overstory transpiration
    ay = plt.subplot(445)
    plot_xy(Data.Transp_av, Data.pine, color=pal[0], axislabels={'x': '', 'y': 'Modelled'})

    ax = plt.subplot(4,4,(6,8))
    plot_timeseries_df(Data, ['pine','Transp_av'], colors=[pal[0],'k'], xticks=False,
                       labels=['Modelled','Measured (sapflow)'], linestyles=['-',':'])
    plt.title('Pine', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    plt.subplot(449)
    plot_xy(Data.Transp_av, Data.spruce, color=pal[0], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(4,4,(10,12), sharex=ax)
    plot_timeseries_df(Data, ['spruce','Transp_av'], colors=[pal[0],'k'], xticks=False,
                       labels=['Modelled','Measured (sapflow)'], linestyles=['-',':'])
    plt.title('Spruce', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    plt.subplot(4,4,13)
    plot_xy(Data.Transp_av, Data.birch, color=pal[0], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(4,4,(14,16), sharex=ax)
    plot_timeseries_df(Data, ['birch','Transp_av'], colors=[pal[0],'k'], xticks=True,
                       labels=['Modelled','Measured (sapflow)'], linestyles=['-',':'])
    plt.title('birch', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    plt.subplot(441)
    plot_xy(Data.Transp_av, Data.shrubs, color=pal[0], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(4,4,(2,4), sharex=ax)
    plot_timeseries_df(Data, ['shrubs','Transp_av'], colors=[pal[0],'k'], xticks=False,
                       labels=['Modelled','Measured (sapflow)'], linestyles=['-',':'])
    plt.title('shrubs', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0,w_pad=0.1, h_pad=0.1)