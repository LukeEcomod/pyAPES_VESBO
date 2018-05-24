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

def fill_gaps(df, res_col_name, plot=False):
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
    info = "\n" + res_col_name + " flags:"
    flag = res_col_name + "_flag"
    NN = float(len(df))
    for i in range(len(col_names)):
        col_name = col_names[i]
        if i == 0:
            df[res_col_name] = df[col_name]
            df[flag] = np.where(df[res_col_name].notnull(), i+1, 0)
            info += "\nFlag %s (%.2f" % (i+1, df[res_col_name].notnull().sum()/NN * 100) + "%): " + col_name
        else:
            df[flag] = np.where(df[res_col_name].isnull() & df[col_name].notnull(), i+1, df[flag])
            df[res_col_name] = df[res_col_name].fillna(df[col_name])
            info += "\nFlag %s (%.2f" % (i+1, sum(df[flag]==i+1)/NN * 100) + "%): " + col_name
    info += "\nFlag %s (%.2f" % (0, df[res_col_name].isnull().sum()/NN * 100) + "%): still missing!!!"
    # HOW TO FILL REST!??
    if plot:
        df[[res_col_name, flag]].plot(subplots=True)
        plt.title(info, y=-0.6,  fontsize=9)
    return df[[res_col_name, flag]], info

def continuous_prec(Prec_data):

    """ --- continuous prec timeseries for jokioinen --- """
    # FMI open data jokioinen + FMI jokioinen gauge 1
    Prec_data['jokioinen_prec'] = Prec_data['jokioinen_meteo: Precipitation amount'].fillna(Prec_data['jokioinen_prec1: Prec [mm h-1]'])
    Prec_data['jokioinen_prec_flag']=np.where(Prec_data['jokioinen_prec'].isnull(), 10.0, 0.0)
    # + FMI jokioinen gauge 2
    Prec_data['jokioinen_prec'] = Prec_data['jokioinen_prec'].fillna(Prec_data['jokioinen_prec2: Prec [mm h-1]'])
    Prec_data['jokioinen_prec_flag']=np.where(Prec_data['jokioinen_prec'].isnull(), 20.0, Prec_data['jokioinen_prec_flag'])
    # + FMI open data somero
    Prec_data['jokioinen_prec'] = Prec_data['jokioinen_prec'].fillna(Prec_data['somero_prec: Precipitation amount'])
    Prec_data['jokioinen_prec_flag']=np.where(Prec_data['jokioinen_prec'].isnull(), 30.0, Prec_data['jokioinen_prec_flag'])
    # fill rest with zero
    Prec_data['jokioinen_prec'] = Prec_data['jokioinen_prec'].fillna(0)
    # plot
    #plot_columns(Prec_data[['jokioinen_prec_flag','jokioinen_prec']])
    
    """ --- continuous prec timeseries for lettosuo --- """
    Prec_data['lettosuo_prec_flag']=np.where(Prec_data['Letto1_metsanpohja: avg(Rain (mm))'].isnull(), 10.0, 0.0)
    # consider prec data unrealiable when Tair < 2C
    Prec_data['lettosuo_prec_flag']=np.where(Prec_data['Letto1_metsanpohja: avg(Temp (C))'] < 2.0,
             20.0, Prec_data['lettosuo_prec_flag'])
    Prec_data['lettosuo_prec']=np.where(Prec_data['Letto1_metsanpohja: avg(Temp (C))'] < 2.0,
             np.nan, Prec_data['Letto1_metsanpohja: avg(Rain (mm))'])
    # if rolling 3 day mean prec less than 10% of jokionen rolling mean, remove
    Prec_data_daily = pd.rolling_sum(Prec_data.fillna(0), 3 * 48, 1)
    Prec_data_daily['lettosuo_prec'] = np.where(Prec_data['lettosuo_prec'].isnull(), np.nan, Prec_data_daily['lettosuo_prec'])
    Prec_data['lettosuo_prec_flag']=np.where(Prec_data_daily['lettosuo_prec'] < 0.1 * Prec_data_daily['jokioinen_prec'],
             30.0, Prec_data['lettosuo_prec_flag'])
    Prec_data['lettosuo_prec']=np.where(Prec_data_daily['lettosuo_prec'] < 0.1 * Prec_data_daily['jokioinen_prec'],
             np.nan, Prec_data['lettosuo_prec'])
    # fill gaps with jokioinen data
    Prec_data['lettosuo_prec'] = Prec_data['lettosuo_prec'].fillna(Prec_data['jokioinen_prec'])
    # plot
    #plot_columns(Prec_data[['lettosuo_prec_flag','lettosuo_prec']])

    return Prec_data