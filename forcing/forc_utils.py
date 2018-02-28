# -*- coding: utf-8 -*-
"""
UTILITIES FOR READING AND EDITING FORCING DATA
"""

import pandas as pd
import numpy as np
from datetime import datetime
import sys

def read_forcing(forc_filename, start_time, end_time, cols=None, dt=None):
    """
    Reads forcing data
    Args:
        forc_filename (str): forcing file name with comma separator
        start_time (str): stranting time [yyyy-mm-dd]
        end_time (str): ending time [yyyy-mm-dd]
        cols (list): header names to read from file, if none reads input data
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
        cols = [
                'Prec',  # Precipitation [mm/dt]
                'Tair',  # Air temperature [degC]
                'U',  # Wind speed 10 min avg. [m/s]
                'RH',  # Relative humidity [%]
                'Rg'  # Global radiation [W/m2]
                ]

    # Forc dataframe from necessary variables
    Forc = dat[cols]

    # Check time step if specified
    if dt is not None:
        if len(set(Forc.index[1:]-Forc.index[:-1])) > 1:
            sys.exit("Forcing file does not have constant time step")
        if (Forc.index[1] - Forc.index[0]).total_seconds() != dt:
            sys.exit("Forcing file time step differs from dt given in general parameters")

    if cols[0] == 'Prec':
        # Calculate forcings needed by model
        Forc = edit_forcing(Forc)

    return Forc

def edit_forcing(Forc):
    """
    Calculate forcings needed by model from given meteorology
    Args:
        Forc (dataframe)
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

    # precipitaion unit from [mm/dt] to [m/s]
    dt = (Forc.index[1] - Forc.index[0]).total_seconds()
    Forc.loc[:,'Prec'] = Forc.loc[:,'Prec'] / 1000.0 / dt

    # daily temperature
    Forc.loc[:,'Tdaily'] = pd.rolling_mean(Forc['Tair'], int((24*3600)/dt), 1)

    # vapor pressure deficit [kPa]
    esat = 0.6112 * np.exp((17.67 * Forc.loc[:,'Tair']) / (Forc.loc[:,'Tair'] + 273.16 - 29.66))  # [kPa]
    Forc.loc[:,'vpd'] = (1 - Forc.loc[:,'RH'] / 100.0) * esat

    # photosynthetically active radiation [Wm-2]
    Forc.loc[:,'Par'] = 0.5 * Forc.loc[:,'Rg']  # Should coefficient be model parameter?

    return Forc