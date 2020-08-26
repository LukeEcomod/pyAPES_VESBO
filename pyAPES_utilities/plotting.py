# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:09:06 2018

Note:
    migrated to python3
    - absolute imports
    - range() is not wrapped in list as consumed in for-each-loop

@author: Kersti Haahti
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import minimize

from pyAPES_utilities.timeseries_tools import diurnal_cycle, yearly_cumulative

EPS = np.finfo(float).eps

prop_cycle = plt.rcParams['axes.prop_cycle']
default = prop_cycle.by_key()['color']
pal=default

def plot_fluxes(results, Data,
                res_var=['canopy_Rnet','canopy_SH','canopy_LE'],
                Data_var=['NRAD','SH','LE'],
                sim_idx=0, fmonth=4, lmonth=9, norain=True, dataframe=False, l1=False,
                save_to_path=None
               ):

    N=len(Data_var)

    if norain:
        res_var.append('forcing_precipitation')

    if dataframe is False:
        df = xarray_to_df(results, set(res_var), sim_idx=sim_idx)
        Data = Data.merge(df, how='outer', left_index=True, right_index=True)
        res_dates=results.date.values
    else:
        Data = Data.merge(results[res_var], how='outer', left_index=True, right_index=True)
        res_dates=results.index

    dates = Data.index

    dryc = np.ones(len(dates))
    if norain:
        ix = Data['forcing_precipitation'].rolling(48, 1).sum()
        f = np.where(ix > EPS)[0]  # wet canopy indices
        dryc[f] = 0.0

    months = Data.index.month

    labels = {'x': '', 'y': 'Modelled'}

    fig = plt.figure(figsize=(10,N*1.7 + 0.5))
    for i in range(N):
        if norain:
            ix = np.where((months >= fmonth) & (months <= lmonth) & (dryc == 1) & np.isfinite(Data[Data_var[i]]))[0]
        else:
            ix = np.where((months >= fmonth) & (months <= lmonth) & np.isfinite(Data[Data_var[i]]))[0]
        plt.subplot(N, 4, i*4+1)
        if i + 1 == N:
            labels = {'x': 'Measured', 'y': 'Modelled'}
        plot_xy(Data[Data_var[i]][ix], Data[res_var[i]][ix], color=pal[i], axislabels=labels,l1=l1)
        if i == 0:
            ax1 = plt.subplot(N,4,(i*4+2,i*4+3))
            plot_timeseries_df(Data, [res_var[i],Data_var[i]], colors=[pal[i],'k'], xticks=False,
                       labels=['Modelled', 'Measured'], marker=[None, '.'])
        else:
            ax = plt.subplot(N,4,(i*4+2,i*4+3), sharex=ax1)
            plot_timeseries_df(Data, [res_var[i],Data_var[i]], colors=[pal[i],'k'], xticks=True,
                       labels=['Modelled', 'Measured'], marker=[None, '.'])
            plt.xlim([res_dates[0], res_dates[-1]])
        try:
            plt.title(results[res_var[i]].units, fontsize=10)
        except:
            plt.title(res_var[i], fontsize=10)
        plt.legend(bbox_to_anchor=(1.6,0.5), loc="center left", frameon=False, borderpad=0.0)
        if i + 1 < N:
            plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
            plt.xlabel('')
        if i == 0:
            ax2 = plt.subplot(N,4,i*4+4)
        else:
            plt.subplot(N,4,i*4+4, sharex=ax2)
        plot_diurnal(Data[Data_var[i]][ix], color='k', legend=False)
        plot_diurnal(Data[res_var[i]][ix], color=pal[i], legend=False)
        if i + 1 < N:
            plt.setp(plt.gca().axes.get_xticklabels(), visible=False)
            plt.xlabel('')

    plt.tight_layout(rect=(0, 0, 0.88, 1), pad=0.5)

    if save_to_path:
        plt.savefig(save_to_path)

    return fig

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
            try:
                labels.append(result[var].units.split('[')[0])
            except:
                labels.append(var)
        title = ''

    try:
        unit = '[' + results[0][variables[0]].units.split('[')[-1]
    except:
        unit='unit'
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
            plt.plot(results[0].date.values, values_all[i], color=colors[i], linewidth=1.5, label=labels[i])
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
            plt.plot(data.index, values_all[i], color=colors[i], linewidth=1.5, label=labels[i], marker=markerstyle, markersize=1, linestyle=linestyle)
        ymax = max([np.nanmax(val) for val in values_all])

    if limits:
        plt.xlim([data.index[0], data.index[-1]])
        plt.ylim(min([np.nanmin(val) for val in values_all]), ymax)

    plt.ylabel(unit)
    plt.setp(plt.gca().axes.get_xticklabels(), visible=xticks)
    plt.xlabel('')
    if legend:
        plt.legend(bbox_to_anchor=(1.01,0.5), loc="center left", frameon=False, borderpad=0.0, fontsize=8)

def plot_columns(data, col_index=None, slope=None, plot_timeseries=True):
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
    if plot_timeseries:
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

def plot_lad_profiles(filename="letto2016_partial.txt", normed=False, quantiles = [1.0],
                      subplot=1, subplots=1, biomass_function='marklund'):
    """
    Plots stand leaf area density profiles from given file.
    Args:
        filename (str): name of file with dbh series
        normed (boolean): normalized profiles
        quantiles (list): cumulative frequency limits for grouping trees, e.g. [0.5, 1.0]
    """
    from .parameter_utilities import model_trees

    z = np.linspace(0, 30.0, 100)
    lad_p, lad_s, lad_d, _, _, _, lai_p, lai_s, lai_d = model_trees(z, quantiles,
        normed=False, dbhfile="pyAPES_utilities/runkolukusarjat/" + filename, plot=False, biomass_function=biomass_function)

    lad = z * 0.0
    for k in range(len(quantiles)):
        lad += lad_p[:,k] + lad_s[:,k] +lad_d[:,k]
    lai_tot = sum(lai_p) + sum(lai_s) + sum(lai_d)

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
#    plt.figure(figsize=(3.2,4.5))
    ax=plt.subplot(1,subplots,subplot)
    if normed:
        for k in range(len(quantiles)):
            plt.plot(lad_p[:, k]/lai_tot,z,color=colors[0], label=r'Pine, $%.2f\times \mathrm{LAI_{tot}}$' % (lai_p[k]/lai_tot))#(lai_p[k]/lai_tot))#,lad_g,z)
            plt.plot(lad_s[:, k]/lai_tot,z,color=colors[1], label=r'Spruce, $%.2f\times \mathrm{LAI_{tot}}$' % (lai_s[k]/lai_tot))
            plt.plot(lad_d[:, k]/lai_tot,z,color=colors[2], label=r'Birch, $%.2f\times \mathrm{LAI_{tot}}$' % (lai_d[k]/lai_tot))
        plt.title("  ")#dbhfile.split("/")[-1])
        plt.ylabel('Height [m]')
        plt.xlabel(r'Normalized leaf area density [m$^2$m$^{-3}$]')
        ax.set_xticks([0.0,0.05,0.1,0.15])
        plt.plot(lad/lai_tot, z,':k', label='Total')
    else:
        for k in range(len(quantiles)):
            plt.plot(lad_p[:, k],z,color=colors[0], label='Pine, %.2f m$^2$m$^{-2}$' % lai_p[k])#,lad_g,z)
            plt.plot(lad_s[:, k],z,color=colors[1], label='Spruce, %.2f m$^2$m$^{-2}$' % lai_s[k])
            plt.plot(lad_d[:, k],z,color=colors[2], label='Birch, %.2f m$^2$m$^{-2}$' % lai_d[k])
        plt.title("  ")#dbhfile.split("/")[-1])
        plt.ylabel('Height [m]')
        plt.xlabel('Leaf area density [m$^2$m$^{-3}$]')
        plt.plot(lad, z,':k', label='Total, %.2f m$^2$m$^{-2}$' % lai_tot)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend(frameon=False, borderpad=0.0, labelspacing=0.1, loc="upper right",bbox_to_anchor=(1.1,1.1))
    plt.tight_layout()

def plot_xy(x, y, color=default[0], title='', axislabels={'x':'', 'y':''},slope=None,alpha=0.05,l1=False):
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
    plt.scatter(x, y, marker='o', color=color, alpha=alpha)
    idx = np.isfinite(x) & np.isfinite(y)

    if l1:
        p = l1_fit(x[idx], y[idx])
    else:
        if slope==None:
            p = np.polyfit(x[idx], y[idx], 1)
        else:
            p = [slope, np.mean(slope * y[idx] - x[idx])]

    residuals = y[idx] - (p[0]*x[idx] + p[1])
    R2 = 1 - sum(residuals**2)/sum((y[idx]-np.mean(y[idx]))**2)
    MAE = np.mean(np.abs(residuals))
#    plt.annotate("y = %.2fx + %.2f\nR$^2$ = %.2f\nMAE = %.1f" % (p[0], p[1], R2, MAE), (0.45, 0.85), xycoords='axes fraction', ha='center', va='center', fontsize=9)
    plt.annotate("y = %.2fx + %.2f\nR$^2$ = %.2f" % (p[0], p[1], R2), (0.45, 0.85), xycoords='axes fraction', ha='center', va='center', fontsize=9)

    lim = [min(min(y[idx]),
               min(x[idx]))-1000,
           max(max(y[idx]),
               max(x[idx]))+1000]
    lim2 = [min(min(y[idx]),
               min(x[idx])),
           max(max(y[idx]),
               max(x[idx]))]
    add = (lim2[1] - lim2[0]) * 0.1
    lim2[0] = lim2[0] - add
    lim2[1] = lim2[1] + add
    plt.plot(lim, [p[0]*lim[0] + p[1], p[0]*lim[1] + p[1]], 'r', linewidth=1)
    plt.plot(lim, lim, 'k--', linewidth=1)
    plt.ylim(lim2)
    plt.xlim(lim2)
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
    print('Daily sum :' + str(sum(var_diurnal['median'])) + ' mm/d')
    plt.title(title)
    plt.xticks(np.linspace(0,24,4))
    plt.xlabel('Time [h]')
    plt.xlim(0,24)
    plt.ylabel(ylabel)
    if legend:
        plt.legend()
        plt.legend(frameon=False, borderpad=0.0,loc="center right")

def xarray_to_df(results, variables, sim_idx=0):
    # note: works only for 1d variables
    series = []
    for var in variables:
        #print(var)
        series.append(results[var].isel(simulation=sim_idx).to_pandas())
    df = pd.concat(series, axis=1)
    #print(np.shape(df), len(variables), df.columns)
    df.columns = variables
    df.index = df.index.round('30T')

    return df

def l1_fit(U, v):
    """
    Find a least absolute error solution (m, k) to U * m + k = v + e.
    Minimize sum of absolute values of vector e (the residuals).
    Returned result is a dictionary with fit parameters result["m"] and result["k"]
    and other information.
    source:
    https://github.com/flatironinstitute/least_absolute_regression/blob/master/lae_regression/lae_regression/least_abs_err_regression.py
    """

    def fit(x, params):
        y = params[0] * x + params[1]
        return y

    def cost_function(params, x, y):
        return np.sum(np.abs(y - fit(x, params)))

    output = minimize(cost_function, np.array([1,1]), args=(U, v))

    return output.x
