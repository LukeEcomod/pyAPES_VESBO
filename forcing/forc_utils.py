# -*- coding: utf-8 -*-
"""
UTILITIES FOR READING AND EDITING FORCING DATA
"""

import pandas as pd
import numpy as np
from datetime import datetime
import sys
from canopy.radiation import solar_angles
from canopy.evapotranspiration import e_sat

def read_forcing(forc_filename, start_time, end_time, 
                 MLM=False, loc=None, cols=None, dt=None, na_values=-999):
    """
    Reads forcing data
    Args:
        forc_filename (str): forcing file name with comma separator
        start_time (str): stranting time [yyyy-mm-dd]
        end_time (str): ending time [yyyy-mm-dd]
        MLM (boolean): True for multilayer computation - needed when reading inputs
        loc (dict): site location - needed when reading inputs
            'lat' (float): latitude
            'lon' (float): longitude
        cols (list): header names to read from file, if None reads input data
        dt (float): time step [s], if given checks that dt in file is equal to this
    Returns:
        Forc (dataframe): dataframe with datetime as index and cols read from file
    """

    # filepath
    forc_fp = "forcing/" + forc_filename
    dat = pd.read_csv(forc_fp, sep=',', header='infer', na_values=-999)

    # set to dataframe index
    dat.index = pd.to_datetime({'year': dat['yyyy'],
                                'month': dat['mo'],
                                'day': dat['dd'],
                                'hour': dat['hh'],
                                'minute': dat['mm']})

    dat = dat[(dat.index >= start_time) & (dat.index <= end_time)]

    # read inputs if cols is not defined
    if cols is None:
        if MLM is False:  # read inputs needed by simple capony description
            cols = ['Prec',  # Precipitation [mm/dt]
                    'Tair',  # Air temperature [degC]
                    'RH',  # Relative humidity [%]  (to compute VPD)
                    'Rg']  # Global radiation [W/m2]
            if 'P' in dat:  cols.append('P')  # Athmospheric pressure [kPa]
            if 'CO2' in dat: cols.append('CO2')  # CO2 mixing ration [ppm]
            if 'U' in dat: cols.append('U')  # Wind speed 10 min avg. [m/s]
            if 'Par' in dat:  cols.append('Par')  # photosynthetically active radiation [W/m2]
        else:  # read inputs needed by multilayer canopy description
            cols = ['Prec',  # Precipitation [mm/dt]
                    'Tair',  # Air temperature [degC]
                    'RH',  # Relative humidity [%]  (to compute VPD)
                    'P',  # Athmospheric pressure [kPa]
                    'U',  # Wind speed 10 min avg. [m/s]
                    'Ustar',  # friction velocity [m/s]
                    'H2O',  # H2O mixing ration [ppth]
                    'CO2',  # CO2 mixing ration [ppm]
                    'dirPar',  # direct Par [W/m2]
                    'diffPar',  # diffuse Par [W/m2]
                    'dirNir',  # direct Nir [W/m2]
                    'diffNir',  # diffuse Nir [W/m2]
                    'LWin',  # longwave radiation [w/m2]
                    'LWout',
                    'Wh',  # top layer soil moisture [m3/m3]   ----> FROM SOIL WATER MODEL
                    'Wa',  # deeper layer soil moisture [m3/m3]   ----> FROM SOIL WATER MODEL
                    'Tsa']  # soil temperature [degC]      ----> FROM SOIL HEAT MODEL

    # Forc dataframe from specified columns
    Forc = dat[cols]

    # Check time step if specified
    if dt is not None:
        if len(set(Forc.index[1:]-Forc.index[:-1])) > 1:
            sys.exit("Forcing file does not have constant time step")
        if (Forc.index[1] - Forc.index[0]).total_seconds() != dt:
            sys.exit("Forcing file time step differs from dt given in general parameters")

    if cols[0] == 'Prec':
        # Calculate forcings needed by model
        Forc = edit_forcing(Forc, MLM, loc)

    return Forc

def edit_forcing(Forc, MLM, loc):
    """
    Calculate forcings needed by model from given meteorology
    Args:
        Forc (dataframe)
        MLM (boolean): True for multilayer computation
        loc (dict): site location
            'lat' (float): latitude
            'lon' (float): longitude
    Returns:
        Forc (dataframe) with edited and new attributes:
            'doy': day of year [days]
            'Prec': converts units to [m s-1]
            'Tdaily': daily temperature as rolling mean [degC]
            'vpd': vapor pressure deficit [kPa]
            'Par': fotosynthetically active radiation [W m-2]
    """

    # day of year
    Forc.loc[:,'doy'] = Forc.index.dayofyear
    jday = Forc.index.dayofyear + Forc.index.hour / 24.0 + Forc.index.minute / 1440.0

    # zenith angle
    Forc.loc[:,'Zen'], _, _, _, _, _ = solar_angles(loc['lat'], loc['lon'], jday, timezone=+2.0)

    # default values for P, CO2, and U if not in inputs
    if 'P' not in Forc:  Forc.loc[:,'P'] = 101.3  # [kPa]
    if 'CO2' not in Forc: Forc.loc[:,'CO2'] = 380.0  # [ppm]
    if 'U' not in Forc: Forc.loc[:,'U'] = 2.0  # [m/s]

    # precipitaion unit from [mm/dt] to [m/s]
    dt = (Forc.index[1] - Forc.index[0]).total_seconds()
    Forc.loc[:,'Prec'] = Forc.loc[:,'Prec'] * 1e-3 / dt
    # atm. pressure [kPa] -> [Pa]
    Forc.loc[:,'P'] = Forc.loc[:,'P'] * 1e3

    # daily temperature
    Forc.loc[:,'Tdaily'] = pd.rolling_mean(Forc['Tair'], int((24*3600)/dt), 1)

    # vapor pressure deficit [kPa]
    esat, _, _ = e_sat(Forc.loc[:,'Tair'],Forc.loc[:,'P'])  # [Pa]
    Forc.loc[:,'vpd'] = (1 - Forc.loc[:,'RH'] / 100.0) * esat * 1e-3  # [kPa]

    if MLM is False:
        # estimate Par [Wm-2] if not in inputs
        if 'Par' not in Forc: Forc.loc[:,'Par'] = 0.5 * Forc.loc[:,'Rg']

    else:
        # H2O mixing ratio [ppth] -> [mol/mol]
        Forc.loc[:,'H2O'] = 1e-3 * Forc.loc[:,'H2O']
        # approximate net radiation [W/m2]
        Forc.loc[:,'Rnet'] = (1.0 - 0.08) * (Forc.loc[:,'dirPar'] + 
                                             Forc.loc[:,'diffPar'] + 
                                             Forc.loc[:,'dirNir'] + 
                                             Forc.loc[:,'diffNir'] + 
                                             Forc.loc[:,'LWin'] - 
                                             Forc.loc[:,'LWout'])

    return Forc

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
        print '********** computing diurnal cycles *********'
    else:
        print('diurnal_cycle: data must be pd.DataFrame or pd.Series')

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