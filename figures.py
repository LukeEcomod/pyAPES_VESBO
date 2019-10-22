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

import seaborn as sns
pal = sns.color_palette("hls", 6)
EPS = np.finfo(float).eps

def plot_vegetation(save=False):

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

    fp = r'O:\H-levy\Lettosuo\aluskasvillisuus\regression_data.txt'
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
                    axislabels={'x':'coverage [%]', 'y':''}, title='Herbs')
    plt.subplot(1,4,3)
    plt.annotate('c', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_regression(dat['shr_cov'], dat['shr_bm'],
                    axislabels={'x':'coverage [%]', 'y':''}, title='Shrubs')
    plt.subplot(1,4,4)
    plt.annotate('d', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plot_regression(dat['bm'], dat['LAI'], alpha=0.5, title='All',
                    axislabels={'x':'biomass [g m$^{-2}$]', 'y':'LAI [m$^2$ m$^{-2}$]'})
    plt.tight_layout()
    if save:
        plt.savefig('figures/case_Lettosuo/vege_regressions.png',dpi=300, transparent=False)

    fp = r'O:\H-levy\Lettosuo\aluskasvillisuus\inventory_data_clc_mounding.txt'
    dat = pd.read_csv(fp,index_col=0)
    labels=['SP$_{\mathrm{all},2009}$', 'VP$_{\mathrm{ref},2015}$', 'VP$_{\mathrm{par},2015}$', 'VP$_{\mathrm{clc},2015}$', 'SP$_{\mathrm{ref},2017}$', 'VP$_{\mathrm{ref},2017}$', 'SP$_{\mathrm{par},2017}$', 'VP$_{\mathrm{par},2017}$', 'SP$_{\mathrm{par},2018}$', 'VP$_{\mathrm{par},2018}$', 'VP$_{\mathrm{clc},2017}$', 'VP$_{\mathrm{clc},2018}$']
    pos = (-0.16,1.02)
    width = 0.75
    plt.figure(figsize=(7.5,4))
    ax1=plt.subplot(2,1,1)
    plt.annotate('a', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    plt.annotate('reference', (1.8,70))
    plt.annotate(' partial harvest', (6.2,70))
    plt.annotate('clear-cut', (9.8,70))
    dat[['seedlings','shrubs','graminoid','herbs']].plot(kind='bar', stacked=True, ax=ax1, colors=vege, width=width)
    plt.plot([5.5, 5.5],[0,1110],'--k')
    plt.plot([9.5, 9.5],[0,1110],'--k')
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.ylim([0,75])
    ax1.set_yticks([0, 25, 50, 75])
    plt.ylabel('coverage [%]',labelpad=10.5)
    ax1.legend().set_visible(False)
    ax1_1 = ax1.twinx()
    ax1_1.plot(range(12),dat['LAI'],'ok', markersize=4)
#    ax1_1.plot([8, 11],[0.88, 1.52],'xk', markersize=4)
    plt.ylabel('LAI [m$^2$ m$^{-2}$]',labelpad=10)
    plt.ylim([0,2])
    ax1.legend().set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1_1.spines['top'].set_visible(False)

    ax2=plt.subplot(2,1,2)
    plt.annotate('b', pos, xycoords='axes fraction', fontsize=12, fontweight='bold')
    dat[['moss','litter','baresoil','seedlings','shrubs','graminoid','herbs']].plot(kind='bar', stacked=True, ax=ax2, colors=ff, width=width)
    plt.plot([1, 2],[110,110],'ok', label='LAI estimated', markersize=4)
#    plt.plot([1, 2],[110,110],'xk', label='LAI measured', markersize=4)
    handles, labels1 = ax2.get_legend_handles_labels()
    labels1[1:]=labels1[:0:-1]
    handles[1:]=handles[:0:-1]
    ax2.legend(handles, labels1,bbox_to_anchor=(1.02,0.35), loc="center left", frameon=False, borderpad=0.0)
    plt.plot([5.5, 5.5],[0,1110],'--k')
    plt.plot([9.5, 9.5],[0,1110],'--k')
    plt.ylabel('coverage [%]')
    plt.ylim([0,100])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_xticklabels(labels)
    plt.xticks(rotation=65)
    plt.tight_layout()
    if save:
        plt.savefig('figures/case_Lettosuo/vege_inventories.png',dpi=300, transparent=False)

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

def plot_wtd(results, fig=True):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_timeseries_xr

    # Read observed WTD
    WTD = read_forcing("Lettosuo_WTD_pred.csv", cols='all')

    if fig:
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

    plot_timeseries_xr(results, 'soil_ground_water_level', colors=['k','b','r','g','m'], xticks=True)

def plot_Tsoil(results, site='Letto1', sim_idx=0,fmonth=5, lmonth=9,l1=True):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_fluxes, xarray_to_df

    # Read observed
    Data = read_forcing("Lettosuo_Tsoil_2010_2018.csv", cols='all')
    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[site]

    df = xarray_to_df(results.isel(soil=4), ['soil_temperature'], sim_idx=sim_idx)
    df['soil_temperature_5cm'] = df['soil_temperature'].copy()
    df2 = xarray_to_df(results.isel(soil=19), ['soil_temperature'], sim_idx=sim_idx)
    df['soil_temperature_30cm'] = df2['soil_temperature'].copy()

    plot_fluxes(df, Data, norain=False,
                res_var=['soil_temperature_5cm','soil_temperature_30cm'],
                Data_var=['T5','T30'],fmonth=fmonth, lmonth=lmonth, sim_idx=sim_idx, dataframe=True,l1=l1)

#plot_Tsoil(results.sel(date=(results['date.year']<2016)), site='Letto1', sim_idx=0)
#plot_Tsoil(results.sel(date=(results['date.year']>=2016)), site='Letto1', sim_idx=1)
#plot_Tsoil(results.sel(date=(results['date.year']>=2016)), site='Clear', sim_idx=2)

def plot_snow_runoff(results, sim_idx=0):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_timeseries_xr, plot_timeseries_df

    # Read observed WTD
    snow_depth = read_forcing("Lettosuo_meteo_2010_2018.csv", cols=['Snow_depth1','Snow_depth2','Snow_depth3'])

    # Read observed WTD
    Data = read_forcing("lettosuo_weir_data.csv", cols=['Runoff mm/h'])

    results['runoff'] = (results['soil_drainage'].copy() + results['soil_surface_runoff'].values) * 1e3 * 3600

    plt.figure()
    ax=plt.subplot(3,1,2)
    plot_timeseries_xr(results.sel(simulation=sim_idx), 'runoff', colors=['r'])
    plot_timeseries_df(Data, 'Runoff mm/h', colors='k', linestyle=':', xticks=False)
    plt.ylabel('[mm/h]')
    plt.subplot(3,1,3,sharex=ax)
    plot_wtd(results.sel(simulation=sim_idx), fig=False)
    plt.subplot(3,1,1,sharex=ax)
    plot_timeseries_df(snow_depth, ['Snow_depth1','Snow_depth2','Snow_depth3'], colors=['gray','gray','gray'])
    plot_timeseries_xr(results.sel(simulation=sim_idx), 'ffloor_snow_water_equivalent',unit_conversion={'unit':'mm', 'conversion':1e3}, xticks=False)


def plot_energy(results,treatment='control',fmonth=5, lmonth=9,sim_idx=0,norain=True,l1=True):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_fluxes

    results['ground_heat_flux'] = results['soil_heat_flux'].isel(soil=6).copy()

    Data = read_forcing("Lettosuo_EC_2010_2018.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[treatment]
    if treatment == 'partial':
        Data['LE'][
            (Data.index > '05-22-2018') & (Data.index < '06-09-2018')]=np.nan
#     period end (?)
    Data.index = Data.index - pd.Timedelta(hours=0.5)
    if treatment=='control':
        Data['NLWRAD'] = Data['NRAD'] - Data['NSWRAD']
        plot_fluxes(results, Data, norain=norain,
                res_var=['canopy_net_radiation','canopy_LWnet','canopy_SWnet', 'canopy_SH','canopy_LE','ground_heat_flux'],
                Data_var=['NRAD','NLWRAD','NSWRAD','SH','LE','GHF'],fmonth=fmonth, lmonth=lmonth, sim_idx=sim_idx,l1=l1)
    elif treatment=='partial':

        Data['NLWRAD'] = Data['NRAD'] - Data['NSWRAD']
        plot_fluxes(results, Data, norain=norain,
                res_var=['canopy_net_radiation','canopy_LWnet','canopy_SWnet', 'canopy_SH','canopy_LE','ground_heat_flux'],
                Data_var=['NRAD','NLWRAD','NSWRAD','SH','LE','GHF'],fmonth=fmonth, lmonth=lmonth, sim_idx=sim_idx,l1=l1)
    else:

        plot_fluxes(results, Data, norain=norain,
                res_var=['canopy_net_radiation', 'canopy_SH','canopy_LE','ground_heat_flux'],
                Data_var=['NRAD','SH','LE','GHF'],fmonth=fmonth, lmonth=lmonth, sim_idx=sim_idx,l1=l1)

def plot_fluxes_ebal(results,treatment='control',fmonth=5, lmonth=9,sim_idx=0,norain=True,l1=True):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_fluxes

    Data = read_forcing("Lettosuo_EC_2010_2018.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[treatment]
    # period end (?)
    Data.index = Data.index - pd.Timedelta(hours=0.5)
    Data['NEE'] *= 1.0 / 44.01e-3
    Data['GPP'] *= -1.0 / 44.01e-3
    Data['Reco'] *= 1.0 / 44.01e-3
    Data['GPP2'] = Data['Reco'] - Data['NEE']
    plot_fluxes(results, Data, norain=norain,
                res_var=['canopy_NEE','canopy_GPP','canopy_respiration','canopy_GPP','canopy_LE'],
                Data_var=['NEE','GPP','Reco','GPP2','LE'],fmonth=fmonth, lmonth=lmonth, sim_idx=sim_idx,l1=l1)

def plot_CO2(results,treatment='control',fmonth=5, lmonth=9,sim_idx=0,norain=True,l1=True):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_fluxes

    Data = read_forcing("Lettosuo_EC_2010_2018.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[treatment]
    # period end (?)
    Data.index = Data.index - pd.Timedelta(hours=0.5)
    Data['NEE'] *= 1.0 / 44.01e-3
    Data['GPP'] *= -1.0 / 44.01e-3
    Data['Reco'] *= 1.0 / 44.01e-3
    Data['GPP2'] = Data['Reco'] - Data['NEE']
    plot_fluxes(results, Data, res_var=['canopy_NEE','canopy_GPP','canopy_respiration'],
                Data_var=['NEE','GPP2','Reco'], norain=norain, fmonth=fmonth, lmonth=lmonth, sim_idx=sim_idx,l1=l1)

def plot_scatters(results, fyear=2010, lyear=2015, treatment='control',fmonth=5, lmonth=9, sim_idx=0, norain=True, legend=True, l1=True):

    from tools.iotools import read_forcing
    import matplotlib.dates
    from pyAPES_utilities.plotting import xarray_to_df, plot_xy, plot_timeseries_xr
    import seaborn as sns
    pal = sns.color_palette("hls", 6)

    Data = read_forcing("Lettosuo_EC_2010_2018.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[treatment]
    # period end (?)
    Data.index = Data.index - pd.Timedelta(hours=0.5)
    Data['NEE'] *= 1.0 / 44.01e-3
    Data['Reco'] *= 1.0 / 44.01e-3
    Data['GPP'] = Data['Reco'] - Data['NEE']
    Data['GPP'][Data['GPP']<0]=np.nan

    Data_var = ['NRAD', 'LE', 'SH', 'GHF','GPP','Reco']
    if treatment == 'partial':
        for var in Data_var:
            Data[var][
                (Data.index > '05-22-2018') & (Data.index < '06-09-2018')]=np.nan
    results['ground_heat_flux'] = results['soil_heat_flux'].isel(soil=6).copy()
    res_var=['canopy_net_radiation', 'canopy_LE', 'canopy_SH', 'ground_heat_flux', 'canopy_GPP','canopy_respiration']

#    if treatment == 'control':
#        Data_var.remove('GHF')
#        res_var.remove('ground_heat_flux')
#        pal.remove(pal[3])

    res_var.append('forcing_precipitation')

    df = xarray_to_df(results, set(res_var), sim_idx=sim_idx)
    Data = Data.merge(df, how='outer', left_index=True, right_index=True)

    dates = Data.index

    ix = Data['forcing_precipitation'].rolling(48, 1).sum()

    dryc = np.ones(len(dates))
    f = np.where(ix > 0.0)[0]  # wet canopy indices
    dryc[f] = 0.0

    months = Data.index.month
    year = Data.index.year

    years = range(fyear, lyear+1)
    N = len(years)
    M = len(Data_var)

    plt.figure(figsize=(N*2 + 0.5,M*2 + 0.5))
    for j in range(M):
        for i in range(N):
            if norain:
                ix = np.where((year == years[i]) & (months >= fmonth) & (months <= lmonth) & (dryc == 1) & np.isfinite(Data[Data_var[j]]))[0]
            else:
                ix = np.where((year == years[i]) & (months >= fmonth) & (months <= lmonth) & np.isfinite(Data[Data_var[j]]))[0]
            if i == 0:
                ax=plt.subplot(M, N, j*N+i+1)
                labels = {'x': 'Measured ' + Data_var[j], 'y': 'Modelled ' + Data_var[j]}
            else:
                labels = {'x': 'Measured ' + Data_var[j], 'y': ''}
                plt.subplot(M, N, j*N+i+1,sharex=ax, sharey=ax)
                plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
            if len(ix) > 0:
                plot_xy(Data[Data_var[j]][ix], Data[res_var[j]][ix], color=pal[j], axislabels=labels, l1=l1)
            if j == 0:
                plt.title(str(years[i]))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)

    plt.figure(figsize=(N*2 + 0.5,5))

    ax = plt.subplot(2, N, (1,N))
    plot_ET(results.sel(date=(results['date.year']>=fyear) & (results['date.year']<=lyear)), sim_idx=sim_idx, fmonth=fmonth, lmonth=lmonth, legend=legend, treatment=treatment)

    # Read observed WTD
    WTD = read_forcing("Lettosuo_WTD_pred.csv", cols='all')
    ax = plt.subplot(2, N, (N+1,2*N))

    plt.fill_between(WTD.index, WTD['control_max'].values, WTD['control_min'].values,
                     facecolor='k', alpha=0.3)
    plt.plot(WTD.index, WTD['control'].values,':k', linewidth=1.0)
    if treatment is not 'control':
        plot_timeseries_xr(results.isel(simulation=0), 'soil_ground_water_level', colors=['k'], xticks=True, legend=False)
        if treatment is 'partial':
            plt.fill_between(WTD.index, WTD['partial_max'].values, WTD['partial_min'].values,
                             facecolor='b', alpha=0.3)
            plt.plot(WTD.index, WTD['partial'].values,':b', linewidth=1.0)

        if treatment is 'clearcut':
            plt.fill_between(WTD.index, WTD['clearcut_max'].values, WTD['clearcut_min'].values,
                             facecolor='r', alpha=0.3)
            plt.plot(WTD.index, WTD['clearcut'].values,':r', linewidth=1.0)

    colors = {'control':  'k', 'partial': 'b', 'clearcut': 'r'}
    plot_timeseries_xr(results.isel(simulation=sim_idx), 'soil_ground_water_level', colors=[colors[treatment]], xticks=True, legend=False)
    plt.ylim([-1.0,0.0])
    plt.xlim(['01-01-'+str(fyear),'12-31-'+str(lyear)])

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(interval=6))
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%Y'))
    plt.ylabel('Water table level [m]')

    plt.tight_layout()

def plot_ET(results, sim_idx=0, fmonth=5, lmonth=9, legend=True, treatment='control'):

    try:
        ptnames=list(results['canopy_planttypes'].values)
    except:
        ptnames=['pine', 'spruce','decid','shrubs']
        print("no plant type names found")

    if sim_idx > 0:
        WB_ref={}
        WB_ref['year']=results['date.year'].sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').mean().values
        WB_ref['interception evaporation']=(results['canopy_evaporation'].isel(simulation=0).sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values + results['canopy_condensation'].isel(simulation=0).sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values)*1800*1000
        WB_ref['transpiration']=results['canopy_transpiration'].isel(simulation=0).sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values*1800*1000
        WB_ref['forest floor evaporation']=results['ffloor_evaporation'].isel(simulation=0).sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values*1800*1000
        WB_ref = pd.DataFrame.from_dict(WB_ref)

    WB={}
    WB['year']=results['date.year'].sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').mean().values
    WB['precipitation']=results['forcing_precipitation'].isel(simulation=sim_idx).sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values*1800*1000
    WB['interception evaporation']=(results['canopy_evaporation'].isel(simulation=sim_idx).sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values + results['canopy_condensation'].isel(simulation=sim_idx).sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values)*1800*1000
    WB['pine transpiration']=results['canopy_pt_transpiration'][:,sim_idx,ptnames.index('pine')].sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values*1800*1000
    WB['spruce transpiration']=results['canopy_pt_transpiration'][:,sim_idx,ptnames.index('spruce')].sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values*1800*1000
    WB['birch transpiration']=results['canopy_pt_transpiration'][:,sim_idx,ptnames.index('decid')].sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values*1800*1000
    WB['understory transpiration']=results['canopy_pt_transpiration'][:,sim_idx,ptnames.index('shrubs')].sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values*1800*1000
    WB['forest floor evaporation']=results['ffloor_evaporation'].isel(simulation=sim_idx).sel(date=(results['date.month']>=fmonth) & (results['date.month']<=lmonth)).groupby('date.year').sum().values*1800*1000

    WB = pd.DataFrame.from_dict(WB)

    # set width of bar
    barWidth = 0.15
    # Set position of bar on X axis
    r1 = np.arange(len(WB['year'])) - barWidth
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    # Make the plot
    plt.bar(r1, WB['interception evaporation'], color=pal[3], width=barWidth, edgecolor='white', label='Interception E')
    plt.bar(r2, WB['pine transpiration'], color=vege[3], width=barWidth, edgecolor='white', label='Pine T')
    plt.bar(r2, WB['spruce transpiration'], bottom=WB['pine transpiration'], color=vege[2], edgecolor='white', width=barWidth, label='Spruce T')
    plt.bar(r2, WB['birch transpiration'], bottom=WB['spruce transpiration']+WB['pine transpiration'], color=vege[1], edgecolor='white', width=barWidth, label='Birch T')
    plt.bar(r2, WB['understory transpiration'], bottom=WB['birch transpiration']+WB['spruce transpiration']+WB['pine transpiration'], color=vege[0], edgecolor='white', width=barWidth, label='Understory T')
    plt.bar(r3, WB['forest floor evaporation'], color=pal[4], width=barWidth, edgecolor='white', label='Forest floor E')

    plt.xticks([r for r in range(len(WB['year']))], [int(WB['year'][i]) for i in range(len(WB['year']))])
    plt.ylabel('ET components [mm]')
    if legend:
        plt.legend(loc='upper center', ncol=3, frameon=False, borderpad=0.0, borderaxespad=-3)

    if treatment is not 'control':
        plt.bar(r1, WB_ref['interception evaporation'], fill=False, width=barWidth, edgecolor='black', linestyle=':')
        plt.bar(r2, WB_ref['transpiration'], fill=False, width=barWidth, edgecolor='black', linestyle=':')
        plt.bar(r3, WB_ref['forest floor evaporation'], fill=False, width=barWidth, edgecolor='black', linestyle=':')

    plt.xlim([-0.5, len(WB['year']) - 0.5])
#    plt.tight_layout(rect=(0, 0, 0.95, 1))
    return WB

def plot_WUE(results,treatment='control',fmonth=5, lmonth=9,sim_idx=0,norain=True):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_fluxes

    Data = read_forcing("Lettosuo_EC_2010_2018.csv", cols='all')

    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[treatment]
    # period end (?)
    Data.index = Data.index - pd.Timedelta(hours=0.5)
    Data['NEE'] *= 1.0 / 44.01e-3
    Data['GPP'] *= -1.0 / 44.01e-3
    Data['Reco'] *= 1.0 / 44.01e-3
    Data['GPP2'] = Data['Reco'] - Data['NEE']

    Forc = read_forcing("Lettosuo_forcing_2010_2018.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)

    plt.figure()
    plt.plot(Data['GPP2'], Data['LE'],'.r')

def data_analysis(results, fmonth=5, lmonth=9,norain=True):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_fluxes, xarray_to_df
    from canopy.micromet import e_sat

    Data = read_forcing("Lettosuo_EC_2010_2018.csv", cols='all',
                        start_time="2010-01-01", end_time="2019-01-01")

    Data['partial_LE'][
        (Data.index > '05-22-2018') & (Data.index < '06-09-2018')]=np.nan
#     period end (?)
    Data.index = Data.index - pd.Timedelta(hours=0.5)

    Forc = read_forcing("Lettosuo_forcing_2010_2018.csv", cols='all',
                    start_time="2010-01-01", end_time="2019-01-01")
    # vapor pressure
    esat, s = e_sat(Forc['Tair'].values)
    Forc['VPD'] = 1e-3 * (esat - Forc['H2O'].values * Forc['P'].values)

    Data = Data.merge(Forc, how='outer', left_index=True, right_index=True)

    dryc = np.ones(len(Data.index))
#    if norain:
#        ix = Forc['Prec'].rolling(48, 1).sum()
#        f = np.where(ix > EPS)[0]  # wet canopy indices
#        dryc[f] = 0.0

    Data['RH'] = Data['H2O'] * Data['P'].values / esat * 100
    f = np.where(Data['RH'] > 70)[0]  # wet canopy indices
    dryc[f] = 0.0
    f = np.where(Data['diffPar'] + Data['dirPar'] < 100)[0]  # wet canopy indices
    dryc[f] = 0.0

    months = Data.index.month

    if norain:
        ix = np.where((months >= fmonth) & (months <= lmonth) & (dryc == 1))[0]
    else:
        ix = np.where((months >= fmonth) & (months <= lmonth))[0]

    df = xarray_to_df(results, ['canopy_LE'], sim_idx=1)
    Data = Data.merge(df, how='outer', left_index=True, right_index=True)
    Data['mod_Gs'] = Data['canopy_LE'] / Data['VPD']

    plt.figure()
    ax=plt.subplot(1,3,1)
    plt.scatter(Data['VPD'][ix],Data['control_LE'][ix],alpha=0.1)
    plt.subplot(1,3,2,sharey=ax)
    plt.scatter(Data['VPD'][ix],Data['partial_LE'][ix],color='r',alpha=0.1)
    plt.subplot(1,3,3,sharey=ax)
    plt.scatter(Data['VPD'][ix],Data['clearcut_LE'][ix],color='g',alpha=0.1)

    Data['control_Gs'] = Data['control_LE']/Data['VPD']
    Data['partial_Gs'] = Data['partial_LE']/Data['VPD']
    Data['clearcut_Gs'] = Data['clearcut_LE']/Data['VPD']
    plt.figure()
    ax=plt.subplot(1,3,1)
    plt.scatter(Data['VPD'][ix],Data['control_Gs'][ix],alpha=0.1)
    plt.scatter(Data['VPD'][ix],Data['mod_Gs'][ix],color='k',alpha=0.1)
    plt.subplot(1,3,2,sharey=ax)
    plt.scatter(Data['VPD'][ix],Data['partial_Gs'][ix],color='r',alpha=0.1)
    plt.subplot(1,3,3,sharey=ax)
    plt.scatter(Data['VPD'][ix],Data['clearcut_Gs'][ix],color='g',alpha=0.1)

    plt.figure()
    plt.scatter(Data['control_Gs'][ix],Data['mod_Gs'][ix],alpha=0.1)

#    Data_daily=Data.resample('D',how=lambda x: x.values.mean())
#
#    Data_daily.plot(subplots=True,kind='line',marker='o',markersize=1)

def plot_ET_WTD(results, fyear=2010, lyear=2015, treatment='control',fmonth=5, lmonth=9, sim_idx=0, norain=True, legend=True, l1=True):

    from tools.iotools import read_forcing
    import matplotlib.dates
    from pyAPES_utilities.plotting import xarray_to_df, plot_xy, plot_timeseries_xr

    years = range(fyear, lyear+1)
    N = len(years)

    plt.figure(figsize=(N*2 + 0.5,5))

    ax = plt.subplot(2, N, (1,N))
    plot_ET(results.sel(date=(results['date.year']>=fyear) & (results['date.year']<=lyear)), treatment=treatment, sim_idx=sim_idx, fmonth=fmonth, lmonth=lmonth, legend=legend)

    # Read observed WTD
    WTD = read_forcing("Lettosuo_WTD_pred.csv", cols='all')
    ax = plt.subplot(2, N, (N+1,2*N))

    plt.fill_between(WTD.index, WTD['control_max'].values, WTD['control_min'].values,
                     facecolor='k', alpha=0.3)
    plt.plot(WTD.index, WTD['control'].values,':k', linewidth=1.0)
    if treatment is not 'control':
        plot_timeseries_xr(results.isel(simulation=0), 'soil_ground_water_level', colors=['k'], xticks=True, legend=False)
        if treatment is 'partial':
            plt.fill_between(WTD.index, WTD['partial_max'].values, WTD['partial_min'].values,
                             facecolor='b', alpha=0.3)
            plt.plot(WTD.index, WTD['partial'].values,':b', linewidth=1.0)

        if treatment is 'clearcut':
            plt.fill_between(WTD.index, WTD['clearcut_max'].values, WTD['clearcut_min'].values,
                             facecolor='r', alpha=0.3)
            plt.plot(WTD.index, WTD['clearcut'].values,':r', linewidth=1.0)

    colors = {'control':  'k', 'partial': 'b', 'clearcut': 'r'}
    plot_timeseries_xr(results.isel(simulation=sim_idx), 'soil_ground_water_level', colors=[colors[treatment]], xticks=True, legend=False)
    plt.ylim([-1.0,0.0])
    plt.xlim(['01-01-'+str(fyear),'12-31-'+str(lyear)])

    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator(interval=6))
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%Y'))
    plt.ylabel('Water table level [m]')

    plt.tight_layout()

def plot_daily(results,treatment='control',fmonth=5, lmonth=9,sim_idx=0, lim=0.6, l1=False):

    from tools.iotools import read_forcing
    from pyAPES_utilities.plotting import plot_xy, xarray_to_df, plot_timeseries_df

    variables=['NRAD','LE','SH','GHF']
    res_var=['canopy_net_radiation','canopy_LE','canopy_SH','ground_heat_flux']

    results['ground_heat_flux'] = results['soil_heat_flux'].isel(soil=6).copy()

    Data = read_forcing("Lettosuo_EC_2010_2018_gapped.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[treatment]
    if treatment == 'partial':
        Data['LE'][
            (Data.index > '05-22-2018') & (Data.index < '06-09-2018')]=np.nan
    # period end
    Data.index = Data.index - pd.Timedelta(hours=0.5)


    df = xarray_to_df(results, set(res_var), sim_idx=sim_idx)
    Data = Data.merge(df, how='outer', left_index=True, right_index=True)

    Data_daily=Data.resample('D',how=lambda x: x.values.mean())

    for i in range(len(variables)):
        Data_daily[variables[i] + '_filtered'] = np.where(
            Data_daily[variables[i]+'gap'] > lim, np.nan, Data_daily[variables[i]])
        Data_daily[variables[i] + '_err'] = Data_daily[variables[i] + '_filtered'] - Data_daily[res_var[i]]

    plt.figure()
    for i in range(len(variables)):
        if i==0:
            ax=plt.subplot(len(variables),1,i+1)
        else:
            plt.subplot(len(variables),1,i+1,sharex=ax)
        plot_timeseries_df(Data_daily,[res_var[i],variables[i],variables[i] +'_filtered'], colors=pal)

    months = Data_daily.index.month

    ix = np.where((months >= fmonth) & (months <= lmonth))[0]

    plt.figure(figsize=(8.5, 2.5))
    for i in range(len(variables)):
        plt.subplot(1,len(variables),i+1)
        plot_xy(Data_daily[variables[i] +'_filtered'][ix], Data_daily[res_var[i]][ix], color=pal[i], alpha=0.3, l1=l1)
        plt.title(variables[i])
    plt.tight_layout()

    year = Data_daily.index.year

    fyear = min(year)
    lyear = max(year)
    years = range(fyear, lyear+1)
    N = len(years)
    M = len(variables)

    plt.figure(figsize=(N*2 + 0.5,M*2 + 0.5))
    for j in range(M):
        for i in range(N):
            ix = np.where((year == years[i]) & (months >= fmonth) & (months <= lmonth))[0]
            if i == 0:
                ax=plt.subplot(M, N, j*N+i+1)
                labels = {'x': 'Measured ' + variables[j], 'y': 'Modelled ' + variables[j]}
            else:
                labels = {'x': 'Measured ' + variables[j], 'y': ''}
                plt.subplot(M, N, j*N+i+1,sharex=ax, sharey=ax)
                plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
            if len(ix) > 0:
                plot_xy(Data_daily[variables[j] +'_filtered'][ix], Data_daily[res_var[j]][ix], color=pal[j], axislabels=labels, alpha=0.6, l1=l1)
            if j == 0:
                plt.title(str(years[i]))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)

    plt.figure(figsize=(10,4*1.7 + 0.5))
    for i in range(len(variables)):
        if i==0:
            ax=plt.subplot(len(variables),1,i+1)
        else:
            plt.subplot(len(variables),1,i+1,sharex=ax)
        plt.plot(Data_daily.index.dayofyear, Data_daily[variables[i] +'_err'],linestyle='',marker='o', color=pal[i], alpha=0.2)
        plt.plot([0, 365], [0, 0], 'k--', linewidth=1)
        plt.title('Measured ' + variables[i] + ' - Modelled ' + variables[i])
        if i+1 < len(variables):
            plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
        else:
            plt.xlabel('Day of year')
    plt.xlim([fmonth*30.5, (lmonth + 1)*30.5])
    plt.tight_layout()