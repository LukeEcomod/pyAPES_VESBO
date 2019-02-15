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
    pos = (-0.3,1.06)
    pos2 = (-0.09,1.02)
    # Overstory transpiration
    ay = plt.subplot(445)
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_xy(Data.Transp_av, Data.T_mod, color=pal[0], axislabels={'x': '', 'y': 'Modelled'})
    plt.ylabel('Modelled',labelpad=7)
    ax = plt.subplot(4,4,(6,8), sharey=ay)
    plt.fill_between(Data.index, Data['T_mod_max'].values, Data['T_mod_min'].values,
                     facecolor=pal[0], alpha=0.3)
    plt.fill_between(Data.index, Data['Transp_max'].values, Data['Transp_min'].values,
                     facecolor='k', alpha=0.2)
    plot_timeseries_df(Data, ['T_mod','Transp_av'], colors=[pal[0],'k'], xticks=False,
                       labels=['Modelled','Measured'], linestyles=['-',':'])
    plt.title('Canopy transpiration [mm d$^{-1}$]', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Residual ET - Tr
    plt.subplot(449, sharex=ay, sharey=ay)
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_xy(Data.E, Data.E_mod, color=pal[3], axislabels={'x': 'Measured', 'y': 'Modelled'})
    plt.ylabel('Modelled',labelpad=7)
    plt.subplot(4,4,(10,12), sharex=ax, sharey=ay)
    plt.fill_between(Data.index, Data['E_mod_max'].values, Data['E_mod_min'].values,
                     facecolor=pal[3], alpha=0.3)
    plt.fill_between(Data.index, Data['E_max'].values, Data['E_min'].values,
                     facecolor='k', alpha=0.2)
    plot_timeseries_df(Data, ['E_mod','E'], colors=[pal[3],'k'], xticks=False,
                       labels=['Modelled', 'Measured'], linestyles=['-',':'])
    plt.title('Evaporation (with understory transpiration) [mm d$^{-1}$]', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Interception evaporation
    plt.subplot(4,4,(14,16), sharex=ax, sharey=ay)
    plt.annotate('d', pos2, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.fill_between(Data.index, Data['I_mod_max'].values, Data['I_mod_min'].values,
                     facecolor=pal[4], alpha=0.3)
    plot_timeseries_df(Data, ['I_mod'], colors=[pal[4]], xticks=True,
                       labels=['Modelled'])
    plt.title('Interception evaporation [mm d$^{-1}$]', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    # Evapotranspiration
    plt.subplot(441, sharex=ay, sharey=ay)
    plt.annotate('a', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_xy(Data.ET, Data.ET_mod, color=pal[2], axislabels={'x': '', 'y': 'Modelled'})
    plt.xlim([0, 5])
    plt.ylabel('Modelled',labelpad=7)
    plt.subplot(4,4,(2,4), sharex=ax, sharey=ay)
    plt.fill_between(Data.index, Data['ET_mod_max'].values, Data['ET_mod_min'].values,
                     facecolor=pal[3], alpha=0.3)
    plot_timeseries_df(Data, ['ET_mod','ET'], colors=[pal[2],'k'], xticks=False,
                       labels=['Modelled','Measured'], linestyles=['-',':'])
    plt.ylim([0, 5])
    plt.title('Evapotranspiration [mm d$^{-1}$]', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    ay.set_xticks([0, 1, 2, 3, 4, 5])
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

    plt.figure(figsize=(9,5.5))
    pos = (-0.25, 1.01)
    pos2 = (-0.11, 1.01)
    plt.subplot(231)
    plt.annotate('a', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    _, _, bars = plt.errorbar(Data_events.I_av, Data_events.I_mod, xerr=data_err, yerr=mod_err, fmt='o',ecolor=grey, color=grey, alpha=0.5)
    [bar.set_alpha(0.2) for bar in bars]
    plot_xy(Data_events.I_av, Data_events.I_mod, color=(1,1,1))
    plt.ylabel('Modelled interception [mm]', labelpad=7)
    plt.xlabel('Measured interception [mm]')
    plt.xlim([0,7.5])
    plt.ylim([0,7.5])

    axx=plt.subplot(232)
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    _, _, bars = plt.errorbar(Data_events.precip_open_area_mm, Data_events.I_av, yerr=data_err, fmt='o',ecolor=pal[0], color=pal[0], alpha=.5)
    [bar.set_alpha(0.2) for bar in bars]
    plt.xlim([0,30])

    plt.ylabel('Measured interception [mm]', labelpad=7)
    plt.xlabel('Precipitation* [mm]')

    plt.subplot(233, sharey=axx)
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    _, _, bars = plt.errorbar(Data_events.Precip, Data_events.I_mod, yerr=mod_err, fmt='o',ecolor=pal[4], color=pal[4], alpha=.5)
    [bar.set_alpha(0.2) for bar in bars]
    plt.xlim([0,30])
    plt.ylabel('Modelled interception [mm]', labelpad=7)
    plt.xlabel('Precipitation [mm]')
    plt.ylim([0,7.5])

    ax = plt.subplot(2,3,(4,5))
    plt.annotate('d', pos2, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.fill_between(Data2.index, np.cumsum(Data2['I_max'].values), 
                np.cumsum(Data2['I_min'].values), facecolor=pal[0], alpha=0.2)
    plt.fill_between(Data2.index, np.cumsum(Data2['I_mod_max'].values), 
                np.cumsum(Data2['I_mod_min'].values), facecolor=pal[4], alpha=0.2)
    plot_timeseries_df(Data2, ['precip_open_area_mm','Precip','I_av','I_mod'],
                       colors=['k',grey,pal[0],pal[4]], xticks=True,
                       labels=['Precipitation*', 'Precipitation', 'Measured interception', 'Modelled interception'], cum=True,
                       unit_conversion={'unit':'mm', 'conversion':1/1800}, legend=False)
    plt.ylabel('[mm]')
    plt.ylim([0,250])
    plt.xlim(['7.1.2016', '11.1.2016'])

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0, h_pad=2, w_pad=2)
    plt.legend(bbox_to_anchor=(1.05,0.5), loc="center left", frameon=False, borderpad=0.0)

# %% Plot energy 
from tools.iotools import read_forcing, xarray_to_df
import numpy as np
import pandas as pd
from tools.plotting import plot_xy, plot_timeseries_df, plot_diurnal
import seaborn as sns
pal = sns.color_palette("hls", 6)
import matplotlib.dates
from tools.timeseries_tools import diurnal_cycle

def plot_energy(results):

    Data = read_forcing("Svarberget_EC_2014_2016_chi.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    variables=['canopy_SH','canopy_LE','canopy_SWnet','canopy_LWnet','canopy_GPP']
    df = xarray_to_df(results, variables, sim_idx=0)
    Data = Data.merge(df, how='outer', left_index=True, right_index=True)
    df_sim=[]
    for sim in range(27):
        df_sim.append(xarray_to_df(results, variables, sim_idx=sim))

    ixLE = np.where(np.isfinite(Data.LE))[0]
    ixSH = np.where(np.isfinite(Data.SH))[0]
    ixLWnet = np.where(np.isfinite(Data.LWnet))[0]
    ixSWnet = np.where(np.isfinite(Data.SWnet))[0]
    ixGPP = np.where(np.isfinite(Data.GPP))[0]
    labels=['Modelled', 'Measured']

    pos= (-0.38, 1.05)

    # Energy
    plt.figure(figsize=(9,10))
    ay1=plt.subplot(541)
    plt.annotate('a', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_xy(Data.SWnet[ixSWnet], Data.canopy_SWnet[ixSWnet], color=pal[0], alpha=0.2, axislabels={'x': '', 'y': 'Modelled'})

    ay2=plt.subplot(545)
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_xy(Data.LWnet[ixLWnet], Data.canopy_LWnet[ixLWnet], color=pal[1], alpha=0.2, axislabels={'x': '', 'y': 'Modelled'})
    plt.ylabel('Modelled',labelpad=-4)

    ay3=plt.subplot(549)
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_xy(Data.SH[ixSH], Data.canopy_SH[ixSH], color=pal[3], alpha=0.2, axislabels={'x': '', 'y': 'Modelled'})
    plt.ylabel('Modelled',labelpad=-4)

    ay4=plt.subplot(5,4,13)
    plt.annotate('d', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_xy(Data.LE[ixLE], Data.canopy_LE[ixLE], color=pal[2], alpha=0.2, axislabels={'x': '', 'y': 'Modelled'})

    ay5=plt.subplot(5,4,17)
    plt.annotate('e', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_xy(Data.GPP[ixGPP], Data.canopy_GPP[ixGPP], color=pal[5], alpha=0.2, axislabels={'x': 'Measured', 'y': 'Modelled'})
    plt.ylabel('Modelled',labelpad=2)

    ax = plt.subplot(5,4,(2,3), sharey=ay1)
    plot_timeseries_df(Data, ['canopy_SWnet', 'SWnet'], colors=[pal[0],'k'], xticks=False,
                       labels=['Modelled', 'Measured'], marker=[None, '.'], legend=False)
    plt.title('Net shortwave radiation [W m$^{-2}$]', fontsize=10)
#    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(5,4,(6,7), sharex=ax, sharey=ay2)
    plot_timeseries_df(Data, ['canopy_LWnet', 'LWnet'], colors=[pal[1],'k'], xticks=False,
                       labels=['Modelled', 'Measured'], marker=[None, '.'], legend=False)
    plt.title('Net longwave radiation [W m$^{-2}$]', fontsize=10)
#    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(5,4,(10,11), sharex=ax, sharey=ay3)
    plot_timeseries_df(Data, ['canopy_SH','SH'], colors=[pal[3],'k'], xticks=False,
                       labels=labels, marker=[None, '.'], legend=False)
    plt.title('Sensible heat flux [W m$^{-2}$]', fontsize=10)
#    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(5,4,(14,15), sharex=ax, sharey=ay4)
    plot_timeseries_df(Data, ['canopy_LE','LE'], colors=[pal[2],'k'], xticks=False,
                       labels=labels, marker=[None, '.'], legend=False)
    plt.title('Latent heat flux [W m$^{-2}$]', fontsize=10)
#    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(5,4,(18,19), sharex=ax, sharey=ay5)
    plot_timeseries_df(Data, ['canopy_GPP','GPP'], colors=[pal[5],'k'], xticks=True,
                       labels=labels, marker=[None, '.'], legend=False)
    plt.title('Gross primary production [$\mu$mol m$^{-2}$ s$^{-1}$]', fontsize=10)
#    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m'))

    ax =plt.subplot(544)
    diurnal=[]
    for sim in range(27):
        diurnal.append(diurnal_cycle(df_sim[sim].canopy_SWnet[ixSWnet])['canopy_SWnet']['median'])
    diurnal=np.array(diurnal)
    plt.fill_between(range(24), np.amin(diurnal,axis=0), np.amax(diurnal,axis=0),
                         color=pal[0], alpha=0.3)
    plot_diurnal(Data.SWnet[ixSWnet], color='k', legend=False)
    plot_diurnal(Data.canopy_SWnet[ixSWnet], color=pal[0], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')

    ax =plt.subplot(548, sharex=ax)
    diurnal=[]
    for sim in range(27):
        diurnal.append(diurnal_cycle(df_sim[sim].canopy_LWnet[ixLWnet])['canopy_LWnet']['median'])
    diurnal=np.array(diurnal)
    plt.fill_between(range(24), np.amin(diurnal,axis=0), np.amax(diurnal,axis=0),
                         color=pal[1], alpha=0.3)
    plot_diurnal(Data.LWnet[ixLWnet], color='k', legend=False)
    plot_diurnal(Data.canopy_LWnet[ixLWnet], color=pal[1], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')
    
    plt.subplot(5,4,12, sharex=ax)
    diurnal=[]
    for sim in range(27):
        diurnal.append(diurnal_cycle(df_sim[sim].canopy_SH[ixSH])['canopy_SH']['median'])
    diurnal=np.array(diurnal)
    plt.fill_between(range(24), np.amin(diurnal,axis=0), np.amax(diurnal,axis=0),
                         color=pal[3], alpha=0.3)
    plot_diurnal(Data.SH[ixSH], color='k', legend=False)
    plot_diurnal(Data.canopy_SH[ixSH], color=pal[3], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')

    plt.subplot(5,4,16, sharex=ax)
    diurnal=[]
    for sim in range(27):
        diurnal.append(diurnal_cycle(df_sim[sim].canopy_LE[ixLE])['canopy_LE']['median'])
    diurnal=np.array(diurnal)
    plt.fill_between(range(24), np.amin(diurnal,axis=0), np.amax(diurnal,axis=0),
                         color=pal[2], alpha=0.3)
    plot_diurnal(Data.LE[ixLE], color='k', legend=False)
    plot_diurnal(Data.canopy_LE[ixLE], color=pal[2], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')

    plt.subplot(5,4,20, sharex=ax)
    diurnal=[]
    for sim in range(27):
        diurnal.append(diurnal_cycle(df_sim[sim].canopy_GPP[ixGPP])['canopy_GPP']['median'])
    diurnal=np.array(diurnal)
    plt.fill_between(range(24), np.amin(diurnal,axis=0), np.amax(diurnal,axis=0),
                         color=pal[5], alpha=0.3)
    plot_diurnal(Data.GPP[ixGPP], color='k', legend=False)
    plot_diurnal(Data.canopy_GPP[ixGPP], color=pal[5], legend=False)

    plt.tight_layout(w_pad=0.1)

    ay1.set_yticks([0, 400, 800])
    ay1.set_xticks([0, 400, 800])
    ay2.set_yticks([-100, -50, 0])
    ay3.set_yticks([-500, 0, 500])
    ay4.set_yticks([0, 300, 600])
    ay4.set_xticks([0, 300, 600])
    ay5.set_yticks([-25, 0, 25])
    ay5.set_xticks([-25, 0, 25])
# %% Plot forcing

def plot_forcing(results):

    Dat = read_forcing("Svarberget_forcing_2014_2016.csv",cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    pos=( -0.08,1.05)

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
    plt.annotate('a', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')

    # Temperature
    plt.subplot(412, sharex=ax)
    plt.fill_between(Data.index, Data['Tmax'].values, Data['Tmin'].values, facecolor=pal[0], alpha=0.2)
    plot_timeseries_df(Data, ['Tair','Tsoil'], colors=pal[0:], xticks=False,legend=False,
                       labels=['Air', 'Top soil'])
    plt.title('Temperature (degC)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)
    plt.ylim([-5, 25])
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')

    # VPD
    plt.subplot(413, sharex=ax)
    plt.fill_between(Data.index, Data['VPDmax'].values, Data['VPDmin'].values, facecolor=pal[2], alpha=0.3)
    plot_timeseries_df(Data, ['VPD'], colors=[pal[2]], xticks=False,legend=False)
    plt.title('Vapor pressure deficit (kPa)', fontsize=10)
    plt.ylim([0, 2.0])
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')

    # Par
    plt.subplot(414, sharex=ax)
    plot_timeseries_df(Data, ['Par', 'diffPar'], colors=pal[0:], xticks=True,legend=False,
                       labels=['Total', 'Diffuse'])
    plt.title('Photosynthetically active radiation (W/m2)', fontsize=10)
    plt.legend(loc="upper right", frameon=False, borderpad=0.0)
    plt.ylim([0, 170])
    plt.annotate('d', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0)

def plot_forcing_manu(results):

    Dat = read_forcing("Svarberget_forcing_2014_2016.csv",cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    pos=( -0.1,1.04)
    fp = r'C:\Users\L1656\Documents\Git_repos\Modeling_cases\pyAPES_Krycklan_C2\forcing\discharge_C2.csv'
    dat = pd.read_csv(fp, sep=';', header='infer')
    dat.index = pd.to_datetime(dat.ix[:,0], yearfirst=True)
    dat = dat[(dat.index >= results.date[0].values) & (dat.index <= results.date[-1].values)]

    # VPD
    from canopy.micromet import e_sat
    # vapor pressure
    esat, s = e_sat(Dat['Tair'].values)
    Dat['VPD'] = 1e-3 * (esat - Dat['H2O'].values * Dat['P'].values)
    Dat['Prec'] = Dat['Prec'].values * 1000 * 3600 * 24
    Dat['Wliq'] = Dat['Wliq'].values * 100
    Dat['Par'] = (Dat['diffPar'].values + Dat['dirPar'].values)
    Data = Dat.resample('D').mean()
    Data=Data.merge(dat, how='outer', left_index=True, right_index=True)
    Data['Tmax'] = Dat['Tair'].resample('D').max()
    Data['Tmin'] = Dat['Tair'].resample('D').min()
    Data['VPDmax'] = Dat['VPD'].resample('D').max()
    Data['VPDmin'] = Dat['VPD'].resample('D').min()

    y=[-100,200]
    x=pd.to_datetime(['7.1.2016','7.1.2016'])
    plt.figure(figsize=(7,6))

    ax = plt.subplot(311)
    plt.fill_between(Data.index, Data['Tmax'].values, Data['Tmin'].values, facecolor=pal[0], alpha=0.2)
    plot_timeseries_df(Data, ['Tair'], colors=[pal[0]], xticks=False,legend=False,
                       labels=['Air temperature'])
    plt.ylabel('Temperature [$^\mathrm{o}$C]')
    plt.ylim([0, 25])
    plt.annotate('a', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    ax2 = ax.twinx()
    ax2.plot(Data.index, Data['Prec'].values * 0.0 - 1.0, color=pal[0], linewidth=1.5, label='Air temperature')
    ax2.plot(Data.index, Data['U'].values, color=pal[4], linewidth=1.5, linestyle='--', label='Wind speed')
    ax2.legend(loc="upper right", frameon=False, borderpad=0.0)
    ax2.plot(x,y, linewidth=1.5, linestyle=':', color='k')
    plt.ylim([0, 10])
    plt.ylabel('Wind speed [m s$^{-1}$]',labelpad=8)

    ax1 = plt.subplot(312, sharex=ax)
    plt.fill_between(Data.index, Data['VPDmax'].values, Data['VPDmin'].values, facecolor=pal[2], alpha=0.3)
    plot_timeseries_df(Data, ['VPD'], colors=[pal[2], pal[0]], xticks=False,legend=False,
                       labels=['VPD'])
    plt.ylabel('VPD [kPa]',labelpad=13)
    plt.ylim([0, 3])
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    ax2 = ax1.twinx()
    ax2.plot(Data.index, Data['Prec'].values * 0.0 - 1.0, color=pal[2], linewidth=1.5, label='VPD')
    ax2.plot(Data.index, Data['Par'].values, color=pal[1], linewidth=1.5, linestyle='--',label='PAR')
    ax2.legend(loc="upper right", frameon=False, borderpad=0.0)
    ax2.plot(x,y, linewidth=1.5, linestyle=':', color='k')
    plt.ylim([0, 160])
    plt.ylabel('PAR [W m$^{-2}$]',labelpad=2)

    Data['Prec'] = -Data['Prec']
    ax3=plt.subplot(313, sharex=ax)
    plot_timeseries_df(Data, ['Prec'], colors=[pal[3]], xticks=True,legend=False)
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.ylim([-40,0])
    plt.yticks([-40,-30,-20,-10,0])
    plt.gca().axes.set_yticklabels([40,30,20,10,0])
    plt.ylabel('Precipitation [mm d$^{-1}$]')
    ax2 = ax3.twinx()
    ax2.plot(Data.index, Data['Prec'].values * 0.0 - 1.0, color=pal[3], linewidth=1.5, label='Precipitation')
    ax2.plot(Data.index, Data['Q_mm'].values, color=pal[5], linewidth=1.5, linestyle='--', label='Runoff')
    ax2.legend(loc="center right", frameon=False, borderpad=0.0)
    ax2.plot(x,y, linewidth=1.5, linestyle=':', color='k')
    plt.ylim([0, 7])
    plt.ylabel('Runoff [mm d$^{-1}$]',labelpad=14)

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

def scatter_plots(results, data_var='Transp_av', mod_var='T_mod', label='Transpiration [mm d$^{-1}$]',color=pal[0]):

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
    Datt = Dat.resample('D').mean()
    Data['Tmax'] = Datt['Tair'].resample('D').max()
    Data['Tmin'] = Datt['Tair'].resample('D').min()
    Data['VPDmax'] = Datt['VPD'].resample('D').max()
    Data['VPDmin'] = Datt['VPD'].resample('D').min()

    Data = Data.merge(Datt, how='outer', left_index=True, right_index=True)
    Data = Data[Data.index >= "6.17.2016"]

    EPS = np.finfo(float).eps
    ix = np.where(Data['Prec'] > 1)[0]
    ixx = np.where(Data['Prec'] < 1)[0]

    plt.figure(figsize=(7,8))
    ax1 = plt.subplot(4,3,1)
    plt.scatter(Data[data_var][ixx], Data.VPD[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[data_var][ix], Data.VPD[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.ylabel('VPD [kPa]')
    plt.title('Measured')
    plt.subplot(4,3,2, sharey=ax1, sharex=ax1)
    plt.scatter(Data[mod_var][ixx], Data.VPD[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[mod_var][ix], Data.VPD[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plt.title('Modelled')
    ax = plt.subplot(4,3,3, sharey=ax1)
    plt.scatter(Data[data_var][ixx] - Data[mod_var][ixx], Data.VPD[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[data_var][ix] - Data[mod_var][ix], Data.VPD[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plt.title('Measured - Modelled')

    ax2 = plt.subplot(4,3,4, sharex=ax1)
    plt.scatter(Data[data_var][ixx], Data.Tair[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[data_var][ix], Data.Tair[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.ylabel('Tair [degC]')
    plt.subplot(4,3,5, sharey=ax2, sharex=ax1)
    plt.scatter(Data[mod_var][ixx], Data.Tair[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[mod_var][ix], Data.Tair[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plt.subplot(4,3,6, sharey=ax2, sharex=ax)
    plt.scatter(Data[data_var][ixx] - Data[mod_var][ixx], Data.Tair[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[data_var][ix] - Data[mod_var][ix], Data.Tair[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)

    ax3 = plt.subplot(4,3,7, sharex=ax1)
    plt.scatter(Data[data_var][ixx], Data.Par[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[data_var][ix], Data.Par[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.ylabel('Par [W m$^{-2}$]')
    plt.subplot(4,3,8, sharey=ax3, sharex=ax1)
    plt.scatter(Data[mod_var][ixx], Data.Par[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[mod_var][ix], Data.Par[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plt.subplot(4,3,9, sharey=ax3, sharex=ax)
    plt.scatter(Data[data_var][ixx] - Data[mod_var][ixx], Data.Par[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[data_var][ix] - Data[mod_var][ix], Data.Par[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)

    ax4 = plt.subplot(4,3,10, sharex=ax1)
    plt.scatter(Data[data_var][ixx], Data.U[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[data_var][ix], Data.U[ix], marker='x', color=color, alpha=0.4)
    plt.ylabel('wind speed [m s$^{-1}$]')
    plt.subplot(4,3,11, sharey=ax4, sharex=ax1)
    plt.scatter(Data[mod_var][ixx], Data.U[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[mod_var][ix], Data.U[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
    plt.xlabel(label)
    plt.subplot(4,3,12, sharey=ax4, sharex=ax)
    plt.scatter(Data[data_var][ixx] - Data[mod_var][ixx], Data.U[ixx], marker='o', color=color, alpha=0.2)
    plt.scatter(Data[data_var][ix] - Data[mod_var][ix], Data.U[ix], marker='x', color=color, alpha=0.4)
    plt.setp(plt.gca().axes.get_yticklabels(), visible=False)

def plot_ETcomponents_stacked(results):

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

    results['ETu_mod'] = (results['understory_transpiration'] +
                        results['understory_evaporation'] +
                        results['ffloor_evaporation'])

    results['I_mod'] = results['overstory_evaporation']

    results['Precip'] = results['forcing_precipitation']

    variables = ['ET_mod','T_mod', 'ETu_mod','I_mod','Precip']

    series = []
    varnames = []
    for var in variables:
        series.append(results.sel(simulation=0)[var].to_pandas())
        varnames.append(var)
    df = pd.concat(series, axis=1)
    df.columns = varnames
    Data = Data.merge(df, how='outer', left_index=True, right_index=True)
    Data = Data.merge(Data2, how='outer', left_index=True, right_index=True)
    Data = Data.fillna(method='ffill')

    for var in variables:
        Data[var] *= 1000*3600*24  # m/s --> mm/d

    Data = Data.resample('D').mean()

    for var in variables[1:]:
        Data[var] = Data[var] / Data['ET_mod']

    plt.figure(figsize=(9,3))

    ax=plt.subplot(111)
    plot_timeseries_df(Data, ['ETu_mod','T_mod','I_mod'], colors=[pal[3],pal[0],pal[1]], xticks=True,
                       labels=['ETu','T','I'], linestyles=['-',':'], stack=True)
    plt.ylim([0,1])
    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y-%m-%d'))

    plt.tight_layout(rect=(0, 0, 1, 1), pad=1.0,w_pad=0.1, h_pad=0.1)