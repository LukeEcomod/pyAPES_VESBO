# -*- coding: utf-8 -*-
"""
UTILITIES FOR READING AND EDITING FORCING DATA
"""

import pandas as pd
import numpy as np
from datetime import datetime

def read_forcing(forc_filename, start_time, end_time):
    """
    Reads forcing data
    Args:
        forc_filename
        start_time
        end_time
    Returns:
        Forc
        Nsteps
    """

    # filepath
    forc_fp = "forcing/" + forc_filename
    dat = pd.read_csv(forc_fp, sep=',', header='infer', na_values=-999)

    # parse first columns to datetime
    N = len(dat)
    t = []
    for k in range(N):
        t.append(
            datetime(dat['yyyy'].iloc[k],
                     dat['mo'].iloc[k],
                     dat['dd'].iloc[k],
                     dat['hh'].iloc[k],
                     dat['mm'].iloc[k]))

    # set to dataframe index
    dat.index = t
    dat = dat[(dat.index >= start_time) & (dat.index <= end_time)]

    # Forc dataframe from necessary variables
    Forc = dat[['Prec', 'Tair', 'U', 'RH', 'Rg', 'SnowD']]

    # Calculate forcings needed by model
    Forc = edit_forcing(Forc)

    # number of timesteps to simulate
    Nsteps = len(Forc)
    del dat, t, N

    return Forc, Nsteps

def edit_forcing(Forc):
    """
    Calculate forcings needed by model from given meteorology
    Args:
        Forc
    Returns:
        Forc
            doy
            Prec in m/s
            vpd
            Par
    """

    # day of year
    Forc['doy'] = Forc.index.dayofyear

    # precipitaion unit from [mm/dt] to [m/s]
    dt = (Forc.index[1]-Forc.index[0]).total_seconds()
    Forc['Prec'] = Forc['Prec'] / 1000.0 / dt

    # vapor pressure deficit [kPa]
    esat = 0.6112 * np.exp((17.67 * Forc['Tair']) / (Forc['Tair'] + 273.16 - 29.66))  # [kPa]
    Forc['vpd'] = (1 - Forc['RH'] / 100.0) * esat

    # photosynthetically active radiation [Wm-2]
    Forc['Par'] = 0.5 * Forc['Rg']  # Should coefficient be model parameter?

    return Forc