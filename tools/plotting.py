# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:09:06 2018

@author: L1656
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from timeseries_tools import diurnal_cycle, yearly_cumulative
from iotools import read_forcing, xarray_to_df
from canopy.constants import LATENT_HEAT, MOLAR_MASS_H2O, MOLAR_MASS_CO2, PAR_TO_UMOL

import seaborn as sns
pal = sns.color_palette("hls", 6)

prop_cycle = plt.rcParams['axes.prop_cycle']
default = prop_cycle.by_key()['color']

def plot_results(results):

    # Read snow

    # Read weir
    weir = read_forcing("Lettosuo_weir.csv", cols=['runf'])

    # Read gwl
    gwl_meas = read_forcing("lettosuo_WTD_pred.csv", cols=['part','clear','ctrl'])
#    gwl_meas = read_forcing("Lettosuo_gwl.csv", cols=['WT_E','WT_N','WT_W','WT_S'], na_values=-999)

    plt.figure(figsize=(10,9))
    ax = plt.subplot(711)
    plot_timeseries_xr(results, 'forcing_air_temperature', colors=pal, xticks=False)
    plt.subplot(712, sharex=ax)
    plot_timeseries_xr(results, ['forcing_precipitation', 'canopy_throughfall'],
                       colors=pal[3:], xticks=False, unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
    plt.subplot(713, sharex=ax)
    plot_timeseries_xr(results, ['ffloor_evaporation_bryo', 'ffloor_evaporation_soil', 'canopy_transpiration', 'canopy_evaporation',
                       'soil_drainage', 'soil_surface_runoff'],
                       colors=pal, cum=True, stack=True, xticks=False,
                       unit_conversion={'unit':'mm', 'conversion':1e3},)
    plt.subplot(714, sharex=ax)
    plot_timeseries_xr(results, ['canopy_evaporation', 'canopy_transpiration', 'ffloor_evaporation_bryo', 'ffloor_evaporation_soil'],
                       colors=[pal[3]] + [pal[2]] + [pal[0]] + [pal[1]], xticks=False,
                       unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
    plt.subplot(715, sharex=ax)
    plot_timeseries_df(weir, 'runf', unit_conversion = {'unit':'mm h-1', 'conversion':3600},
                       labels='measured runoff', colors=['k'], xticks=False)
    plot_timeseries_xr(results, 'soil_drainage', colors=pal, xticks=False,
                      unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
    plt.subplot(716, sharex=ax)
    plot_timeseries_xr(results, 'ffloor_snow_water_equivalent', colors=['gray'], xticks=False, stack=True,
                       unit_conversion={'unit':'mm', 'conversion':1e3})
    plt.subplot(717, sharex=ax)
    plot_timeseries_df(gwl_meas, ['part','clear','ctrl'], colors=pal[1:], xticks=True, limits=False)
    plot_timeseries_xr(results, 'soil_ground_water_level', colors=pal[0:], xticks=True)
    plt.ylim(-1.0, 0.0)

    plt.tight_layout(rect=(0, 0, 0.8, 1))

def plot_fluxes(results, treatment='control-N', sim_idx=0):
    Data = read_forcing("Lettosuo_EC.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[treatment]
    Data['ET'] = Data.LE / LATENT_HEAT * MOLAR_MASS_H2O * 3600  # W m-2 / J mol-1 * kg mol-1 * s h-1 = kg m-2 h-1 = mm h-1
    Data.GPP = -Data.GPP

    variables=['canopy_NEE','canopy_GPP','canopy_respiration','canopy_transpiration','canopy_evaporation','forcing_precipitation','ffloor_evaporation']
    df = xarray_to_df(results, variables, sim_idx=sim_idx)
    Data = Data.merge(df, how='outer', left_index=True, right_index=True)
    Data['ET_mod'] = (Data.canopy_transpiration + Data.ffloor_evaporation) * 1e3 * 3600
    Data.canopy_GPP = Data.canopy_GPP * MOLAR_MASS_CO2  # umol m-2 s-1 * kg mol-1 = mg m-2 s-1
    Data.canopy_respiration = Data.canopy_respiration * MOLAR_MASS_CO2

    dates = Data.index

    ix = pd.rolling_mean(Data.forcing_precipitation.values, 48, 1)
    dryc = np.ones(len(dates))
    f = np.where(ix > 0)[0]  # wet canopy indices
    dryc[f] = 0.0
    months = Data.index.month
    fmonth = 4
    lmonth = 9
    Data.ET[Data.gapped == 1] = Data.ET[Data.gapped == 1] * np.nan
    ixET = np.where((months >= fmonth) & (months <= lmonth) & (dryc == 1) & np.isfinite(Data.ET))[0]
    Data.GPP[Data.gapped == 1] = Data.GPP[Data.gapped == 1] * np.nan
    ixGPP = np.where((months >= fmonth) & (months <= lmonth) & np.isfinite(Data.GPP))[0]
    Data.Reco[Data.gapped == 1] = Data.Reco[Data.gapped == 1] * np.nan
    ixReco = np.where((months >= fmonth) & (months <= lmonth) & np.isfinite(Data.Reco))[0]
    labels=['Modelled', 'Measured']

    plt.figure(figsize=(10,6))
    plt.subplot(341)
    plot_xy(Data.GPP[ixGPP], Data.canopy_GPP[ixGPP], color=pal[0], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(345)
    plot_xy(Data.Reco[ixReco], Data.canopy_respiration[ixReco], color=pal[1], axislabels={'x': '', 'y': 'Modelled'})

    plt.subplot(349)
    plot_xy(Data.ET[ixET], Data.ET_mod[ixET], color=pal[2], axislabels={'x': 'Measured', 'y': 'Modelled'})

    ax = plt.subplot(3,4,(2,3))
    plot_timeseries_df(Data, ['canopy_GPP', 'GPP'], colors=[pal[0],'k'], xticks=False,
                       labels=labels, marker=[None, '.'])
    plt.title('GPP [mg CO2 m-2 s-1]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(3,4,(6,7), sharex=ax)
    plot_timeseries_df(Data, ['canopy_respiration','Reco'], colors=[pal[1],'k'], xticks=False,
                       labels=labels, marker=[None, '.'])
    plt.title('Reco [mg CO2 m-2 s-1]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    plt.subplot(3,4,(10,11), sharex=ax)
    plot_timeseries_df(Data, ['ET_mod','ET'], colors=[pal[2],'k'], xticks=True,
                       labels=labels, marker=[None, '.'])
    plt.title('ET [mm h-1]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)

    ax =plt.subplot(344)
    plot_diurnal(Data.GPP[ixGPP], color='k', legend=False)
    plot_diurnal(Data.canopy_GPP[ixGPP], color=pal[0], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')

    plt.subplot(348, sharex=ax)
    plot_diurnal(Data.Reco[ixReco], color='k', legend=False)
    plot_diurnal(Data.canopy_respiration[ixReco], color=pal[1], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')

    plt.subplot(3,4,12, sharex=ax)
    plot_diurnal(Data.ET[ixET], color='k', legend=False)
    plot_diurnal(Data.ET_mod[ixET], color=pal[2], legend=False)

    plt.tight_layout(rect=(0, 0, 0.88, 1), pad=0.5)

def plot_pt_results(results, variable):
    if variable == 'canopy_pt_transpiration':
        unit1={'unit':'mm h-1', 'conversion':1e3*3600}
        unit2={'unit':'mm', 'conversion':1e3}
    else:
        unit1={'unit':'mg CO2 m-2 s-1', 'conversion': 44.01 * 1e-3}
        unit2={'unit':'kg CO2 m-2', 'conversion': 44.01 * 1e-9}
    plt.figure(figsize=(10,4))
    ax = plt.subplot(211)
    plot_timeseries_pt(results, variable, sim_idx=0,
                       unit_conversion=unit1,
                       xticks=False, stack=False, cum=False)
    plt.subplot(212, sharex=ax)
    plot_timeseries_pt(results, variable, sim_idx=0,
                       unit_conversion=unit2,
                       xticks=True, stack=True, cum=True)
    plt.title('')
    plt.tight_layout(rect=(0, 0, 0.9, 1))

def plot_timeseries_pt(results, variable, sim_idx=0, unit_conversion={'unit':None, 'conversion':1.0},
                    xticks=True, stack=True, cum=True, legend=True, limits=True, colors=None):
    """
    Plot results by plant type.
    Args:
        results (xarray): xarray dataset
        variables (str): variable name to plot
        sim_idx (int): index of simulation in xarray dataset
        unit_converrsion (dict):
            'unit' (str): unit of plotted variable (if different from unit in dataset)
            'conversion' (float): conversion needed to get to 'unit'
        xticks (boolean): False for no xticklabels
        stack (boolean): stack plotted timeseries (fills areas)
        cum (boolean): plot yearly cumulative timeseries
    """
    species=['Pine','Spruce','Decid']
    labels=[]
    if colors == None:
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
    sub_planttypes = (len(results.planttype) - 1) / 3
    if sub_planttypes > 1:
        for sp in species:
            labels += [sp + '_' + str(k) for k in range(sub_planttypes)]
    else:
        labels = species
    labels.append('Shrubs')
    plot_timeseries_xr([results.isel(planttype=i) for i in range(len(results.planttype))],
                       variable, colors=colors, xticks=xticks,
                       stack=stack, cum=cum, sim_idx=sim_idx,
                       unit_conversion=unit_conversion, labels=labels, legend=legend, limits=limits)

def plot_timeseries_xr(results, variables, unit_conversion = {'unit':None, 'conversion':1.0},
                       colors=default, xticks=True, stack=False, cum=False, labels=None, limits=True, legend=True):
    """
    Plot timeseries from xarray results.
    Args:
        results (xarray or list of xarray): xarray dataset or list of xarray datasets
        variables (str or list of str): variable names to plot
        unit_converrsion (dict):
            'unit' (str): unit of plotted variable (if different from unit in dataset)
            'conversion' (float): conversion needed to get to 'unit'
        colors (list): list of color codes
        xticks (boolean): False for no xticklabels
        stack (boolean): stack plotted timeseries (fills areas)
        cum (boolean): plot yearly cumulative timeseries
        labels (list of str): labels corresponding to results[0], results[1],...
    """
    if type(results) != list:
        if 'simulation' in results.dims:
            results = [results.sel(simulation=i) for i in results.simulation.values]
        else:
            results = [results]
    else:
        results = [result.isel(simulation=0) for result in results]
    if type(variables) != list:
        variables = [variables]

    values_all=[]
    labels=[]
    for result in results:
        if cum:
            values = yearly_cumulative(result, variables)
            for k in range(len(values)):
                values_all.append(values[k])
        else:
            for var in variables:
                values_all.append(result[var].values)
        for var in variables:
            labels.append(result[var].units.split('[')[0])
        title = ''

    unit = '[' + results[0][variables[0]].units.split('[')[-1]
    if cum:
        unit = unit.split(' ')[-2]+']'
    if unit_conversion['unit'] != None:
        unit = '[' + unit_conversion['unit'] + ']'
        values_all = [val * unit_conversion['conversion'] for val in values_all]

    if len(values_all) > len(colors):
        colors = colors * (len(values_all) / len(colors) + 1)

    if stack:
        plt.stackplot(results[0].date.values, values_all, labels=labels, colors=colors)
        ymax = max(sum(values_all))
    else:
        for i in range(len(values_all)):
            plt.plot(results[0].date.values, values_all[i], color=colors[i], linewidth=1, label=labels[i])
        ymax = max([max(val) for val in values_all])
    plt.title(title)
    plt.xlim([results[0].date.values[0], results[0].date.values[-1]])
    if limits:
        plt.ylim(min([min(val) for val in values_all]), ymax)
    plt.ylabel(unit)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=xticks)
    plt.xlabel('')
    if legend:
        plt.legend(bbox_to_anchor=(1.01,0.5), loc="center left", frameon=False, borderpad=0.0, fontsize=8)

def plot_timeseries_df(data, variables, unit_conversion = {'unit':None, 'conversion':1.0},
                       labels=None, marker=None, colors=default, xticks=True, stack=False, cum=False, linestyle='-', limits=True,legend=True):
    """
    Plot timeseries from dataframe data.
    Args:
        data (DataFrame): Dataframe of datasets with datetime index
        variables (str or list of str): variable names to plot
        unit_converrsion (dict):
            'unit' (str): unit of plotted variable (if different from unit in dataset)
            'conversion' (float): conversion needed to get to 'unit'
        labels (str or list of str): labels corresponding to variables (if None uses variable names)
        colors (list): list of color codes
        xticks (boolean): False for no xticklabels
        stack (boolean): stack plotted timeseries (fills areas)
        cum (boolean): plot yearly cumulative timeseries
        limits (boolean): set axis limits based on plotted data
    """
    if type(variables) != list:
        variables = [variables]
    values_all=[]
    if cum:
        values = yearly_cumulative(data, variables)
        for k in range(len(values)):
            values_all.append(values[k])
    else:
        for var in variables:
            values_all.append(data[var].values)
    if labels == None:
        labels = variables
    if type(labels) != list:
        labels = [labels]

    if unit_conversion['unit'] != None:
        unit = '[' + unit_conversion['unit'] + ']'
        values_all = [val * unit_conversion['conversion'] for val in values_all]
    else:
        unit = ''

    if len(values_all) > len(colors):
        colors = colors * (len(values_all) / len(colors) + 1)

    if stack:
        plt.stackplot(data.index, values_all, labels=labels, colors=colors)
        ymax = max(sum(values_all))
    else:
        for i in range(len(values_all)):
            if marker is None:
                markerstyle = None
            else:
                markerstyle = marker[i]
                if marker[i] is not None:
                    linestyle = 'None'
            plt.plot(data.index, values_all[i], color=colors[i], linewidth=1, label=labels[i], marker=markerstyle, markersize=1, linestyle=linestyle)
        ymax = max([np.nanmax(val) for val in values_all])

    if limits:
        plt.xlim([data.index[0], data.index[-1]])
        plt.ylim(min([np.nanmin(val) for val in values_all]), ymax)

    plt.ylabel(unit)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=xticks)
    plt.xlabel('')
    if legend:
        plt.legend(bbox_to_anchor=(1.01,0.5), loc="center left", frameon=False, borderpad=0.0, fontsize=8)

def plot_columns(data, col_index=None, slope=None):
    """
    Plots specified columns from dataframe as timeseries and correlation matrix.
    Args:
        data (DataFrame): dataframe including data to plot
        col_index (list of int): indices of columns to plot as list
            (if None plot all colums in data)
    """
    col_names=[]
    if col_index == None:
        col_index = range(len(data.columns))
    for i in col_index:
        col_names.append(data.columns[i])
    data[col_names].plot(kind='line',marker='o',markersize=1)
    plt.legend()
    axes = pd.plotting.scatter_matrix(data[col_names], figsize=(10, 10), alpha=.2)
    corr = data[col_names].corr().as_matrix()
    lim = [data[col_names].min().min(), data[col_names].max().max()]
    for i in range(len(col_index)):
        for j in range(len(col_index)):
            if i != j:
                idx = np.isfinite(data[col_names[i]]) & np.isfinite(data[col_names[j]])
                if idx.sum() != 0.0:
                    if slope==None:
                        p = np.polyfit(data[col_names[j]][idx], data[col_names[i]][idx], 1)
                        R2 = corr[i,j]**2
                    else:
                        p = [slope, np.mean(data[col_names[i]][idx] - slope* data[col_names[j]][idx])]
                        SStot = np.sum((data[col_names[i]][idx] - np.mean(data[col_names[i]][idx]))**2)
                        SSres = np.sum((data[col_names[i]][idx] - (p[0]*data[col_names[j]][idx] + p[1]))**2)
                        R2 = 1 - SSres/SStot
                    axes[i, j].annotate("y = %.2fx + %.2f \n R2 = %.3f" % (p[0], p[1], R2), (0.4, 0.9), xycoords='axes fraction', ha='center', va='center')
                    axes[i, j].plot(lim, [p[0]*lim[0] + p[1], p[0]*lim[1] + p[1]], 'r', linewidth=1)
                axes[i, j].plot(lim, lim, 'k--', linewidth=1)
            axes[i, j].set_ylim(lim)
            axes[i, j].set_xlim(lim)

def plot_lad_profiles(filename="letto2016_partial.txt", normed=False, quantiles = [1.0]):
    """
    Plots stand leaf area density profiles from given file.
    Args:
        filename (str): name of file with dbh series
        normed (boolean): normalized profiles
        quantiles (list): cumulative frequency limits for grouping trees, e.g. [0.5, 1.0]
    """
    from parameters.utilities import model_trees

    z = np.linspace(0, 30.0, 100)
    lad_p, lad_s, lad_d, _, _, _, lai_p, lai_s, lai_d = model_trees(z, quantiles,
        normed=normed, dbhfile="parameters/runkolukusarjat/" + filename, plot=True)
    if normed == False:
        lad = z * 0.0
        for k in range(len(quantiles)):
            lad += lad_p[:,k] + lad_s[:,k] +lad_d[:,k]
        lai_tot = sum(lai_p) + sum(lai_s) + sum(lai_d)
        plt.plot(lad, z,':k', label='total, %.2f m$^2$m$^{-2}$' % lai_tot)
    plt.legend(frameon=False, borderpad=0.0, labelspacing=0.1)

def plot_xy(x, y, color=default[0], title='', axislabels={'x':'', 'y':''}):
    """
    Plot x,y scatter with linear regression line, info of relationship and 1:1 line.
    Args:
        x,y (array): arrays for x and y data
        color (str or tuple): color code
        title (str): title of plot
        axislabels (dict)
            'x' (str): x-axis label
            'y' (str): y-axis label
    """
    plt.scatter(x, y, marker='o', color=color, alpha=.1)
    idx = np.isfinite(x) & np.isfinite(y)
    p = np.polyfit(x[idx], y[idx], 1)
    corr = np.corrcoef(x[idx], y[idx])
    plt.annotate("y = %.2fx + %.2f \nR2 = %.2f" % (p[0], p[1], corr[1,0]**2), (0.45, 0.85), xycoords='axes fraction', ha='center', va='center', fontsize=9)
    lim = [min(min(y[np.isfinite(y)]),
               min(x[np.isfinite(x)])),
           max(max(y[np.isfinite(y)]),
               max(x[np.isfinite(x)]))]
    plt.plot(lim, [p[0]*lim[0] + p[1], p[0]*lim[1] + p[1]], 'r', linewidth=1)
    plt.plot(lim, lim, 'k--', linewidth=1)
    plt.ylim(lim)
    plt.xlim(lim)
    plt.title(title)
    plt.xlabel(axislabels['x'])
    plt.ylabel(axislabels['y'])

def plot_diurnal(var, quantiles=False, color=default[0], title='', ylabel='', label='median', legend=True):
    """
    Plot diurnal cycle (hourly) for variable
    Args:
        var (Dataframe): dataframe with datetime index and variable to plot
        quantiles (boolean): plot quantiles as ranges
        color (str or tuple): color code
        title (str): title of plot
        ylabel (str): y-axis label
    """
    
    var_diurnal = diurnal_cycle(var)[var.name]
    if quantiles:
        plt.fill_between(var_diurnal['hour'], var_diurnal['5th'], var_diurnal['95th'],
                         color=color, alpha=0.2, label='5...95%')
        plt.fill_between(var_diurnal['hour'], var_diurnal['25th'], var_diurnal['75th'],
                         color=color, alpha=0.4, label='25...75%')
    plt.plot(var_diurnal['hour'], var_diurnal['median'],'o-', markersize=3, color=color,label=label)
    print 'Daily sum :' + str(sum(var_diurnal['median'])) + ' mm/d'
    plt.title(title)
    plt.xticks(np.linspace(0,24,4))
    plt.xlabel('Time [h]')
    plt.xlim(0,24)
    plt.ylabel(ylabel)
    if legend:
        plt.legend()
        plt.legend(frameon=False, borderpad=0.0,loc="center right")

def plot_efficiencies(results, treatment='control', sim_idx=0):
    PAR = read_forcing("Lettosuo_forcing_2010_2018.csv",cols=['diffPar','dirPar'],
                       start_time=results.date[0].values, end_time=results.date[-1].values)
    Data = read_forcing("Lettosuo_EC.csv", cols='all',
                        start_time=results.date[0].values, end_time=results.date[-1].values)
    Data.columns = Data.columns.str.split('_', expand=True)
    Data = Data[treatment]
    Data['ET'] = Data.LE / LATENT_HEAT * 1e3 # [mmol m-2 s-1]
    Data.GPP = -Data.GPP / MOLAR_MASS_CO2  # [umol m-2 s-1]
    Data['PAR'] = (PAR['diffPar'] + PAR['dirPar']) * PAR_TO_UMOL  # [umol m-2 s-1]

    variables=['canopy_GPP','canopy_transpiration','canopy_evaporation','forcing_precipitation','ffloor_evaporation']
    df = xarray_to_df(results, variables, sim_idx=sim_idx)
    Data = Data.merge(df, how='outer', left_index=True, right_index=True)
    Data['ET_mod'] = (Data.canopy_transpiration + Data.ffloor_evaporation) / MOLAR_MASS_H2O * 1e3 * 1e3 # [mmol m-2 s-1]
    Data.canopy_GPP = Data.canopy_GPP  # [umol m-2 s-1]

    dates = Data.index

    ix = pd.rolling_mean(Data.forcing_precipitation.values, 48, 1)
    dryc = np.ones(len(dates))
    f = np.where(ix > 0)[0]  # wet canopy indices
    dryc[f] = 0.0
    months = Data.index.month
    fmonth = 4
    lmonth = 9
    hours = Data.index.hour
    fhour = 0
    lhour = 24
    Data.ET[Data.gapped == 1] = Data.ET[Data.gapped == 1] * np.nan

    Data.GPP[Data.gapped == 1] = Data.GPP[Data.gapped == 1] * np.nan

    labels=['Modelled', 'Measured']

    Data['LUE'] = Data.GPP / Data.PAR
    Data['WUE'] = Data.GPP / Data.ET

    Data['LUE_mod'] = Data.canopy_GPP / Data.PAR
    Data['WUE_mod'] = Data.canopy_GPP / Data.ET_mod
    ixET = np.where((months >= fmonth) & (months <= lmonth) &
                    (hours >= fhour) & (hours <= lhour) &
                    (dryc == 1) & np.isfinite(Data.WUE))[0]
    ixGPP = np.where((months >= fmonth) & (months <= lmonth) &
                    (hours >= fhour) & (hours <= lhour) &
                    np.isfinite(Data.LUE))[0]

    plt.figure(figsize=(10,4.5))
    plt.subplot(241)
    plot_xy(Data.LUE[ixGPP], Data.LUE_mod[ixGPP], color=pal[0], axislabels={'x': '', 'y': 'Modelled'})
    plt.ylim(0, 0.05)
    plt.xlim(0, 0.05)

    plt.subplot(245)
    plot_xy(Data.WUE[ixET], Data.WUE_mod[ixET], color=pal[2], axislabels={'x': 'Measured', 'y': 'Modelled'})
    plt.ylim(0, 20)
    plt.xlim(0, 20)
    ax = plt.subplot(2,4,(2,3))
    plot_timeseries_df(Data, ['LUE_mod', 'LUE'], colors=[pal[0],'k'], xticks=False,
                       labels=labels, marker=[None, '.'], limits=False)
    plt.title('Ligth use efficiency [umol/umol]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)
    plt.ylim(0, 0.05)

    plt.subplot(2,4,(6,7), sharex=ax)
    plot_timeseries_df(Data, ['WUE_mod','WUE'], colors=[pal[2],'k'], xticks=True,
                       labels=labels, marker=[None, '.'], limits=False)
    plt.title('Water use efficiency [umol/mmol]', fontsize=10)
    plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)
    plt.ylim(0, 20)

    ax =plt.subplot(244)
    plot_diurnal(Data.LUE[ixGPP], color='k', legend=False)
    plot_diurnal(Data.LUE_mod[ixGPP], color=pal[0], legend=False)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
    plt.xlabel('')

    plt.subplot(248, sharex=ax)
    plot_diurnal(Data.WUE[ixET], color='k', legend=False)
    plot_diurnal(Data.WUE_mod[ixET], color=pal[2], legend=False)

    plt.tight_layout(rect=(0, 0, 0.88, 1), pad=0.5)