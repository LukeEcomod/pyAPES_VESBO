# -*- coding: utf-8 -*-
"""
I/O tools for handling CCFPeat simulations.

Created on Fri Jun 08 10:32:44 2018

@author: Kersti Haahti
"""

import sys
import pandas as pd
import datetime

def read_forcing(forc_filename, start_time=None, end_time=None,
                 cols=None, dt=None, na_values='NaN'):
    """
    Reads forcing or other data from csv file to dataframe
    Args:
        forc_filename (str): forcing file name with comma separator
        start_time (str): starting time [yyyy-mm-dd], if None first date in file used
        end_time (str): ending time [yyyy-mm-dd], if None last date in file used
        cols (list): header names to read from file, if None reads variables need as input data
        dt (float): time step [s], if given checks that dt in file is equal to this
        na_values (str/float): nan value representation in file
    Returns:
        Forc (dataframe): dataframe with datetime as index and cols read from file
    """

    # filepath
    forc_fp = "forcing/" + forc_filename
    dat = pd.read_csv(forc_fp, sep=',', header='infer', na_values=na_values)

    # set to dataframe index
    dat.index = pd.to_datetime({'year': dat['yyyy'],
                                'month': dat['mo'],
                                'day': dat['dd'],
                                'hour': dat['hh'],
                                'minute': dat['mm']})

    if start_time == None:
        start_time = dat.index[0]
    if end_time == None:
        end_time = dat.index[-1]
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
        Forc.loc[:,'Par'] = Forc['diffPar'].values + Forc['dirPar'].values

    return Forc

def read_results(outputfiles):
    """
    Opens simulation results netcdf4 dataset in xarray
    (or multiple in list of xarrays)
    Args:
        outputfiles (str or list of str):
            outputfilenameor list of outputfilenames
    Returns:
        results (xarray or list of xarrays):
            simulation results from given outputfile(s)
    """

    import xarray as xr

    direc = "results/"

    if type(outputfiles) != list:
        outputfiles = [outputfiles]

    results = []
    for outputfile in outputfiles:
        fp = direc + outputfile
        result = xr.open_dataset(fp)
        result.coords['simulation'] = result.simulation.values
        result.coords['soil'] = result.soil_z.values
        result.coords['canopy'] = result.canopy_z.values
        result.coords['planttype'] = ['pine','spruce','decid','shrubs']
        results.append(result)

    if len(results) == 1:
        return results[0]
    else:
        return results

def save_df_to_csv(df, fn, readme='', fp="forcing/", timezone = +2):
    """
    Save dataframe with datetime index to csv file with corresponding readme.txt.
    Args:
        df (DataFrame): data to save
        fn (str): filename of saved file (.csv not included)
        readme (str): readme corresponding to df to save as txt file
        fp (str): filepath, forcing folder used as default
        timezone (float): time zone in refernce to UTC, default UTC + 2
    """

    # add datetime as columns
    df.insert(0, 'yyyy', df.index.year.values)
    df.insert(1, 'mo', df.index.month.values)
    df.insert(2, 'dd', df.index.day.values)
    df.insert(3, 'hh', df.index.hour.values)
    df.insert(4, 'mm', df.index.minute.values)

    df.to_csv(path_or_buf=fp + fn + ".csv", sep=',', na_rep='NaN', index=False)
    Readme = "Readme for " + fn + ".csv"
    Readme += "\n\nKersti Haahti, Luke " + str(datetime.datetime.now().date())
    Readme += "\n\nyyyy, mo, dd, hh, mm: datetime [UTC + %.1f]" % timezone
    Readme += readme
    outF = open(fp + fn + "_readme.txt", "w")
    print >>outF, Readme
    outF.close()

def xarray_to_df(results, variables, sim_idx=0):
    series = []
    for var in variables:
        series.append(results[var].isel(simulation=sim_idx).to_pandas())
    df = pd.concat(series, axis=1)
    df.columns=variables
    return df