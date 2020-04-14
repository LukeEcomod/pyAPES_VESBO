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
            # if key == 'canopy_planttypes':
            #     print(key, length, type(results[key]), results[key], np.shape(ncf[key]))
            if length > 1:
                ncf[key][:, nsim, :] = results[key]
            elif key == 'soil_z' or key == 'canopy_z' or \
                 key == 'canopy_planttypes' or key == 'ffloor_groundtypes':
                if nsim == 0:
                    ncf[key][:] = results[key]
            else:
                ncf[key][:, nsim] = results[key]


def initialize_netcdf(variables,
                      sim,
                      soil_nodes,
                      canopy_nodes,
                      planttypes,
                      groundtypes,
                      time_index,
                      filepath='results/',
                      filename='climoss.nc',
                      description='Simulation results'):
    """ Climoss netCDF4 format output file initialization

    Args:
        variables (list): list of variables to be saved in netCDF4
        sim (int): number of simulations
        soil_nodes (int): number of soil calculation nodes
        canopy_nodes (int): number of canopy calculation nodes
        time_index: time_index of default forcing data (pd.DataSeries)
        filepath: path for saving results
        filename: filename
    """
    from netCDF4 import Dataset, date2num
    from datetime import datetime

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

    ncf.createDimension('date', None)
    ncf.createDimension('simulation', sim)
    ncf.createDimension('soil', soil_nodes)
    ncf.createDimension('canopy', canopy_nodes)
    ncf.createDimension('planttype', planttypes)
    ncf.createDimension('groundtype', groundtypes)
    
    time = ncf.createVariable('date', 'f8', ('date',))
    time.units = 'days since 0001-01-01 00:00:00.0'
    time.calendar = 'standard'
#    tvec = [k.to_datetime() for k in forcing.index] is depricated
    tvec = [pd.to_datetime(k) for k in time_index]
    time[:] = date2num(tvec, units=time.units, calendar=time.calendar)

    for var in variables:

        var_name = var[0]
        var_unit = var[1]
        var_dim = var[2]

        if var_name == 'canopy_planttypes' or var_name == 'ffloor_groundtypes':
            variable = ncf.createVariable(
                var_name, 'S10', var_dim)
        else:
            variable = ncf.createVariable(
                    var_name, 'f4', var_dim)

        variable.units = var_unit

    return ncf, ff


def read_forcing(forc_filename, start_time=None, end_time=None,
                 cols=None, dt=1800.0, na_values='NaN', sep=';'):
    """
    Reads forcing or other data from csv file to dataframe
    Args:
        forc_filename (str): forcing file name
        start_time (str): starting time [yyyy-mm-dd], if None first date in
            file used
        end_time (str): ending time [yyyy-mm-dd], if None last date
            in file used
        cols (list): header names to read from file, if None reads
            variables need as input data
        dt (float): time step [s], if given checks
            that dt in file is equal to this
        na_values (str/float): nan value representation in file
        sep (str): field separator
    Returns:
        Forc (dataframe): dataframe with datetime as index and cols read from file
    """

    # filepath
    forc_fp = "forcing/" + forc_filename
    dat = pd.read_csv(forc_fp, header='infer', na_values=na_values, sep=sep)

    # set to dataframe index
    tvec = pd.to_datetime(dat[['year', 'month', 'day', 'hour', 'minute']])
    tvec = pd.DatetimeIndex(tvec)
    dat.index = tvec

    # select time period    
    if start_time == None:
        start_time = dat.index[0]
    if end_time == None:
        end_time = dat.index[-1]
    
    dat = dat[(dat.index >= start_time) & (dat.index <= end_time)]

    # convert: H2O mmol / mol --> mol / mol; Prec kg m-2 in dt --> kg m-2 s-1
    dat['H2O'] = 1e-3 * dat['H2O']
    dat['Prec'] = dat['Prec'] / dt
    
    # if cols is not defined; return these
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
                'LWin',
                'diffPar',
                'dirPar',
                'diffNir',
                'dirNir']
        # these for phenology model initialization
        if 'X' in dat:
            cols.append('X')
        if 'DDsum' in dat:
            cols.append('DDsum')
    # or return all columns
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

    return Forc

def read_data(ffile, start_time=None, end_time=None, na_values='NaN', sep=';'):
    dat = pd.read_csv(ffile, header='infer', na_values=na_values, sep=sep)
    # set to dataframe index
    tvec = pd.to_datetime(dat[['year', 'month', 'day', 'hour', 'minute']])
    tvec = pd.DatetimeIndex(tvec)
    dat.index = tvec

    # select time period    
    if start_time == None:
        start_time = dat.index[0]
    if end_time == None:
        end_time = dat.index[-1]
    
    dat = dat[(dat.index >= start_time) & (dat.index <= end_time)]

    return dat

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

