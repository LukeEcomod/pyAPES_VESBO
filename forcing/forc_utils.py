# -*- coding: utf-8 -*-
"""
UTILITIES FOR READING AND EDITING FORCING DATA
"""

import pandas as pd
import numpy as np
from datetime import datetime

def read_forcing(forc_filename, start_time, end_time):  # cols=None jos on määritelty lukee inputit muuten määritellyt columnit??
    """
    Reads forcing data
    Args:
        forc_filename
        start_time
        end_time
    Returns:
        Forc
        Nsteps   # tarviiko?
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

    # Forc dataframe from necessary variables
    Forc = dat[['Prec', 'Tair', 'U', 'RH', 'Rg', 'SnowD']]  # yksiköt!

    # Calculate forcings needed by model
    Forc = edit_forcing(Forc)

    # number of timesteps to simulate
    Nsteps = len(Forc)

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
    Forc.loc[:,'doy'] = Forc.index.dayofyear

    # precipitaion unit from [mm/dt] to [m/s]
    dt = (Forc.index[1] - Forc.index[0]).total_seconds()
    Forc.loc[:,'Prec'] = Forc.loc[:,'Prec'] / 1000.0 / dt

    # vapor pressure deficit [kPa]
    esat = 0.6112 * np.exp((17.67 * Forc.loc[:,'Tair']) / (Forc.loc[:,'Tair'] + 273.16 - 29.66))  # [kPa]
    Forc.loc[:,'vpd'] = (1 - Forc.loc[:,'RH'] / 100.0) * esat

    # photosynthetically active radiation [Wm-2]
    Forc.loc[:,'Par'] = 0.5 * Forc.loc[:,'Rg']  # Should coefficient be model parameter?

    return Forc