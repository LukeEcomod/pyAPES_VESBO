# -*- coding: utf-8 -*-
"""
I/O tools for handling CCFPeat simulations.

Created on Fri Jun 08 10:32:44 2018

Note:
    migrated to python3
    - print >>out to print(file=out)

@author: Kersti Haahti
"""

import sys
import os
import pandas as pd
import numpy as np
import datetime
import json


def write_ncf(nsim=None, results=None, ncf=None):
    """ Writes model simultaion results in netCDF4-file

    Args:
        index (int): model loop index
        results (dict): calculation results from group
        ncf (object): netCDF4-file handle
    """

    keys = results.keys()
    variables = ncf.variables.keys()

    for key in keys:

        if key in variables and key != 'time':
            length = np.asarray(results[key]).ndim

            if length > 1:
                ncf[key][:, nsim, :] = results[key]
            elif key == 'soil_z' or key == 'canopy_z':
                if nsim == 0:
                    ncf[key][:] = results[key]
            else:
                ncf[key][:, nsim] = results[key]


def initialize_netcdf(variables,
                      sim,
                      soil_nodes,
                      canopy_nodes,
                      plant_nodes,
                      forcing,
                      filepath='results/',
                      filename='climoss.nc',
                      description='Simulation results'):
    """ Climoss netCDF4 format output file initialization

    Args:
        variables (list): list of variables to be saved in netCDF4
        sim (int): number of simulations
        soil_nodes (int): number of soil calculation nodes
        canopy_nodes (int): number of canopy calculation nodes
        forcing: forcing data (pd.dataframe)
        filepath: path for saving results
        filename: filename
    """
    from netCDF4 import Dataset, date2num
    from datetime import datetime

    # dimensions
    date_dimension = None
    simulation_dimension = sim
    soil_dimension = soil_nodes
    canopy_dimension = canopy_nodes
    ptypes_dimension = plant_nodes

    pyAPES_folder = os.getcwd()
    filepath = os.path.join(pyAPES_folder, filepath)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ff = os.path.join(filepath, filename)

    # create dataset and dimensions
    ncf = Dataset(ff, 'w')
    ncf.description = description
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'pyAPES_beta2018'

    ncf.createDimension('date', date_dimension)
    ncf.createDimension('simulation', simulation_dimension)
    ncf.createDimension('soil', soil_dimension)
    ncf.createDimension('canopy', canopy_dimension)
    ncf.createDimension('planttype', ptypes_dimension)

    time = ncf.createVariable('date', 'f8', ('date',))
    time.units = 'days since 0001-01-01 00:00:00.0'
    time.calendar = 'standard'
#    tvec = [k.to_datetime() for k in forcing.index] is depricated
    tvec = [pd.to_datetime(k) for k in forcing.index]
    time[:] = date2num(tvec, units=time.units, calendar=time.calendar)

    for var in variables:

        var_name = var[0]
        var_unit = var[1]
        var_dim = var[2]

        variable = ncf.createVariable(
                var_name, 'f4', var_dim)

        variable.units = var_unit

    return ncf, ff


def read_forcing(forc_filename, start_time=None, end_time=None,
                 cols=None, dt=None, na_values='NaN'):
    """
    Reads forcing or other data from csv file to dataframe
    Args:
        forc_filename (str): forcing file name with comma separator
        start_time (str): starting time [yyyy-mm-dd], if None first date in
            file used
        end_time (str): ending time [yyyy-mm-dd], if None last date
            in file used
        cols (list): header names to read from file, if None reads
            variables need as input data
        dt (float): time step [s], if given checks
            that dt in file is equal to this
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
    elif cols == 'all':
        cols = [col for col in dat]
    # Forc dataframe from specified columns
    Forc = dat[cols].copy()

    # Check time step if specified
    if dt is not None:
        if len(set(Forc.index[1:]-Forc.index[:-1])) > 1:
            sys.exit("Forcing file does not have constant time step")
        if (Forc.index[1] - Forc.index[0]).total_seconds() != dt:
            sys.exit("Forcing file time step differs from dt given in general parameters")

    if Forc.columns[0] == 'doy':
        Forc.loc[:,'Par'] = Forc['diffPar'].values + Forc['dirPar'].values

    return Forc

def read_results(outputfiles):
    """
    Opens simulation results netcdf4 dataset in xarray
    (or multiple in list of xarrays)
    Args:
        outputfiles (str or list of str):
            outputfilename or list of outputfilenames
    Returns:
        results (xarray or list of xarrays):
            simulation results from given outputfile(s)
    """

    import xarray as xr

    if type(outputfiles) != list:
        outputfiles = [outputfiles]

    results = []
    for outputfile in outputfiles:
        fp = outputfile
        result = xr.open_dataset(fp)
        result.coords['simulation'] = result.simulation.values
        result.coords['soil'] = result.soil_z.values
        result.coords['canopy'] = result.canopy_z.values
#        result.coords['planttype'] = ['pine','spruce','decid','shrubs']
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
    print(Readme, file=outF)
    outF.close()


def xarray_to_df(results, variables, sim_idx=0):
    series = []
    for var in variables:
        series.append(results[var].isel(simulation=sim_idx).to_pandas())
    df = pd.concat(series, axis=1)
    df.columns = variables
    return df


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def jsonify(params, file_name='parameter_space.json'):
    """ Dumps simulation parameters into json format
    """

    os.path.join(file_name)

    with open(file_name, 'w+') as fp:
        json.dump(params, fp, cls=NumpyEncoder)


def open_json(file_path):
    """ Opens a json file
    """
    os.path.join(file_path)
    with open(file_path) as json_data:
        data = json.load(json_data)

    return data

