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
        cols = ['doy',
                'Prec',
                'P',
                'Tair',
                'Tdaily',
                'U',
                'Ustar',
                'H2O',
                'CO2',
                'Zen',
                'Rg',
                'LWin',
                'LWout',
                'diffPar',
                'dirPar',
                'diffNir',
                'dirNir']

    # Forc dataframe from specified columns
    Forc = dat[cols]

    # Check time step if specified
    if dt is not None:
        if len(set(Forc.index[1:]-Forc.index[:-1])) > 1:
            sys.exit("Forcing file does not have constant time step")
        if (Forc.index[1] - Forc.index[0]).total_seconds() != dt:
            sys.exit("Forcing file time step differs from dt given in general parameters")

    if cols[0] == 'doy':
        # Calculate forcings needed by model
        Forc['Par'] = Forc['diffPar'] + Forc['dirPar']
        #Forc = edit_forcing(Forc, MLM, loc)

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

    if MLM is False:
        # estimate Par [Wm-2] if not in inputs
        if 'Par' not in Forc: Forc.loc[:,'Par'] = 0.5 * Forc.loc[:,'Rg']
        # H2O mixing ratio from RH
        esat, _, _ = e_sat(Forc['Tair'])
        Forc.loc[:,'H2O'] = (Forc.loc[:,'RH'] / 100.0) * esat / Forc.loc[:,'P']
    else:
        # H2O mixing ratio [ppth] -> [mol/mol]
        Forc.loc[:,'H2O'] = 1e-3 * Forc.loc[:,'H2O']
        # approximate net radiation [W/m2]
        Forc.loc[:,'Rnet'] = (1.0 - 0.08) * (Forc.loc[:,'dirPar'] \
                + Forc.loc[:,'diffPar'] + Forc.loc[:,'dirNir'] \
                + Forc.loc[:,'diffNir']) + Forc.loc[:,'LWin'] - Forc.loc[:,'LWout']

    return Forc