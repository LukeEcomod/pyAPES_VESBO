# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 10:19:19 2018

@author: L1656
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def yearly_cumulative(results, variables):
    if type(results) is pd.DataFrame:
        years = results.index.year
    else: # xarray
        years = results.date.dt.year.values

    yearly_cum = np.empty([len(variables), len(years)])

    for k in range (0, len(variables)):
        yearly_cum[k,:] = np.cumsum(results[variables[k]].values)

    for t in range (years[0], years[-1]):
        ix = np.where(years > t)
        for i in range(0, len(variables)):
            yearly_cum[i,ix] = yearly_cum[i,ix] - yearly_cum[i,ix[0][0]]
    return yearly_cum

def diurnal_cycle(data, ap='hour'):
    """
    computes ensemble diurnal cycle of flux or environmental data
    Args:
        data - pd.DataFrame or pd.Series. index = pd.datetime, i.e. '1996-01-01 00:30:00'
        ap - averaging period. ap='hour' or 'minute'
    Returns:
        res - dict: keys == data.columns
                    values == pd.dataframe where
                                columns = ['hour', 'minu', 'N', 'mean', 'std', 'se',
                                           'median', '5th', '25th', '75th', '95th'
                                          ]
    NOTE:
        seeks for unique hours and minutes in data, ensembles them and returns statistics.
        Nodata == np.NaN are omited when statistics are computed.
    Samuli Launiainen, Luke Jan 7th, 2018
    """
    
    if isinstance(data, pd.Series):
        data = data.to_frame()
    
    if isinstance(data, pd.DataFrame):
        r, c = np.shape(data)  # rows, cols
        hr = data.index.hour
        mn = data.index.minute
        hour = np.unique(hr)
        minu = np.unique(mn)
        cols = data.columns
        

    else:
        print('diurnal_cycle: data must be pd.DataFrame or pd.Series')

    print '********** computing diurnal cycles *********'
    res = {}
    for k in range(0, c):
        if ap.lower() == 'hour':
            N = len(hour)             
            x = np.ones((N, 11))*np.NaN
            x[:, 0] = hour
            x[:, 1] = 0.0
            
            n = 0
            for t in hour:
                y = data.iloc[:, k]  # column k
                f = np.where((hr == t) & (np.isfinite(y)))[0]
                
                x[n, 2] = len(f)  # no of observations
                x[n, 3] = np.mean(y[f]) 
                x[n, 4] = np.std(y[f])
                x[n, 5] = x[n, 3] / x[n, 2]  # s.e.
                x[n, 6:] = np.percentile(y[f], [50.0, 5.0, 25.0, 75.0, 95.0])
                n += 1

            res[cols[k]] = pd.DataFrame(x, columns=['hour', 'minu', 'N', 'mean', 'std', 'se',
                                                    'median', '5th', '25th', '75th', '95th'])

        if ap.lower() == 'minute':
            N = len(hour) * len(minu)
            x = np.ones((N, 11))*np.NaN

            n = 0
            for t in hour:
                for p in minu:
                    # print(k, t, p)
                    y = data.iloc[:, k]  # column k
                    f = np.where((hr == t) & (mn == p) & (np.isfinite(y)))[0]
                    # print f
                    x[n, 0] = t
                    x[n, 1] = p
                    x[n, 2] = len(f)  # no of observations
                    x[n, 3] = np.mean(y[f])
                    x[n, 4] = np.std(y[f])
                    x[n, 5] = x[n, 3] / x[n, 2]  # s.e.
                    x[n, 6:] = np.percentile(y[f], [50.0, 5.0, 25.0, 75.0, 95.0])
                    n += 1
            res[cols[k]] =  pd.DataFrame(x, columns=['hour', 'minu', 'N', 'mean', 'std', 'se',
                                                     'median', '5th', '25th', '75th', '95th'])

    return res

""" --- Gap filling scripts---"""

def fill_gaps(df, res_col_name, description, fill_nan = None, plot=False):
    """
    Fill gaps with other available data
    Args:
        df (DataFrame): columns in priority order to be used
        res_col_name (string): name of column in dataframe returned by function
    Returns:
        df (Dataframe): resulting dataframe with
            'res_col_name': data collected from input
            'res_col_name_flag': flags defining data source
        info (string): flag definitions and their fequency
    """
    col_names = list(df.columns.values)
    info = "\n" + res_col_name + ": " + description
    flag = res_col_name + "_flag"
    NN = float(len(df))
    for i in range(len(col_names)):
        col_name = col_names[i]
        if i == 0:
            df[res_col_name] = df[col_name]
            df[flag] = np.where(df[res_col_name].notnull(), i, len(col_names))
            info += "\n  flag %s (%.2f" % (i, df[res_col_name].notnull().sum()/NN * 100) + "%): " + col_name
        else:
            df[flag][df[res_col_name].isnull() & df[col_name].notnull()] = i
            df[res_col_name] = df[res_col_name].fillna(df[col_name])
            info += "\n  flag %s (%.2f" % (i, sum(df[flag]==i)/NN * 100) + "%): " + col_name
    if fill_nan == 'linear':
        df = df.interpolate()
        message = "linearly interpolated"
        # fill nans in beginning and end
        if df[res_col_name].isnull().sum() > 0:
            df[flag][df[res_col_name].isnull()] = len(col_names) + 1
            df = df.fillna(method='bfill')
            df = df.fillna(method='ffill')
            message += "\n  flag %s (%.2f" % (i, sum(df[flag]==len(col_names) + 1)/NN * 100) + "%): filled with nearest"
    elif type(fill_nan) == float:
        df = df.fillna(fill_nan)
        message = "filled with " + res_col_name + " = "+ str(fill_nan)
    else:
        message = "still missing!!!"
    info += "\n  flag %s (%.2f" % (len(col_names), sum(df[flag]==len(col_names))/NN * 100) + "%): " + message
    if plot:
        df[[res_col_name, flag]].plot(subplots=True)
        plt.title(info, y=-0.6,  fontsize=9)
    return df[[res_col_name, flag]], info

