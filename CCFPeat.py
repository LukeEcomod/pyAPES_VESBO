# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 15:09:03 2018

@author: L1656
"""
import os
import numpy as np
import pandas as pd
from forcing.forc_utils import read_forcing
from canopy.canopy_model import CanopyModel
from soilprofile.soil_model import SoilModel

from copy import deepcopy

def driver(create_ncf=False, LAI_sensitivity=False):
    """
    """

    # Import general parameters
    from parameters.general_parameters import gpara
    # Import canopy model parameters
    from parameters.canopy_parameters import cpara
    # Import soil model parameters
    from parameters.soil_parameters import spara

    if LAI_sensitivity:
        from parameters.sensitivity_sampling import LAIcombinations
        Nsim = len(LAIcombinations)
    else:
        Nsim = 1

    # Read forcing
    forcing = read_forcing(gpara['forc_filename'],
                        gpara['start_time'],
                        gpara['end_time'],
                        MLM=cpara['ctr']['multilayer_model']['ON'],
                        loc=cpara['loc'],
                        dt=gpara['dt'])

    tasks = []

    for k in range(Nsim):
        if LAI_sensitivity:
            for pt in range(3):
                cpara['plant_types'][pt].update({'LAImax': [LAIcombinations[k][pt]]})
        tasks.append(Model(gpara, cpara, spara, forcing, nsim=k))

    if create_ncf:
        import time
        timestr = time.strftime('%Y%m%d%H%M')
        filename = timestr + '_CCFPeat_results.nc'

        ncf, _ = initialize_netcdf(
                gpara['variables'],
                Nsim,
                tasks[k].Nsoil_nodes,
                tasks[k].Ncanopy_nodes,
                forcing,
                filename=filename)

        for task in tasks:
            print('Running simulation number: {}' .format(task.Nsim))
            running_time = time.time()
            results = task._run()
            print('Running time %.2f seconds' % (time.time() - running_time))
            _write_ncf(nsim=task.Nsim, results=results, ncf=ncf)

            del results

        output_file = "results/" + filename
        print('Ready! Results are in: ' + output_file)
        ncf.close()

    else:
        results = {task.Nsim: task._run() for task in tasks}
        output_file = results

    return output_file

class Model():
    """
    """
    def __init__(self,
                 gen_para,
                 canopy_para,
                 soil_para,
                 forcing,
                 nsim=0):

        self.dt = gen_para['dt']

        self.Nsteps = len(forcing)
        self.forcing = forcing
        self.Nsim = nsim

        self.Nsoil_nodes = len(soil_para['z'])
        if canopy_para['ctr']['multilayer_model']['ON']:
            self.Ncanopy_nodes = canopy_para['grid']['Nlayers']
        else:
            self.Ncanopy_nodes = 1

        # create canopy model instance
        self.canopy_model = CanopyModel(canopy_para)

        # create soil model instance
        self.soil_model = SoilModel(soil_para['z'], soil_para)

        self.results = _create_results(gen_para['variables'],
                                       self.Nsteps,
                                       self.Nsoil_nodes,
                                       self.Ncanopy_nodes)

    def _run(self):
        """ Runs atmosphere-canopy-soil--continuum model"""
        k_steps=np.arange(0, self.Nsteps, int(self.Nsteps/10))
        for k in range(0, self.Nsteps):
            if k in k_steps[:-1]:
                s = str(np.where(k_steps==k)[0][0]*10) + '%'
                print '{0}..\r'.format(s),

            # Soil moisture forcing for canopy model
            Rew = 1.0
            beta = 1.0

            """ Canopy, moss and Snow """
            # run daily loop (phenology and seasonal LAI)
            if self.forcing['doy'].iloc[k] != self.forcing['doy'].iloc[k-1] or k == 0:
                self.canopy_model._run_daily(
                        self.forcing['doy'].iloc[k],
                        self.forcing['Tdaily'].iloc[k])

            # run timestep loop
            canopy_flux, canopy_state = self.canopy_model._run_timestep(
                    dt=self.dt,
                    forcing=self.forcing.iloc[k],
                    beta=beta,
                    Rew=Rew)

            """ Water and Heat in soil """
            # potential infiltration and evaporation from ground surface
            ubc_w = {'Prec': canopy_flux['potential_infiltration'] / self.dt,
                     'Evap': 0.0}
            # transpiration sink  --- SOIL MODELIN SISÄÄN?!?
            rootsink = np.zeros(self.soil_model.Nlayers)
            rootsink[0] = canopy_flux['transpiration'] / self.dt / self.soil_model.dz[0]  # ekasta layerista, ei väliä tasapainolaskennassa..
            # temperature above soil surface
            ubc_T = {'type': 'flux', 'value': None}

            # run soil water and heat flow
            soil_flux, soil_state = self.soil_model._run(
                    self.dt,
                    ubc_w, ubc_T,
                    water_sink=rootsink)

            forcing_state = {
                    'wind_speed': self.forcing['U'].iloc[k],
                    'air_temperature': self.forcing['Tair'].iloc[k],
                    'precipitation': self.forcing['Prec'].iloc[k] * self.dt,
                    'h2o': self.forcing['H2O'].iloc[k],
                    'co2': self.forcing['CO2'].iloc[k]}

            canopy_state.update(canopy_flux)
            soil_state.update(soil_flux)

            self.results = _append_results('canopy', k, canopy_state, self.results)
            self.results = _append_results('soil', k, soil_state, self.results)
            self.results = _append_results('forcing', k, forcing_state, self.results)

        print '100%'
        self.results = _append_results('canopy', None, {'z': self.canopy_model.z}, self.results)
        self.results = _append_results('soil', None, {'z': self.soil_model.z}, self.results)
        return self.results

def _create_results(variables, Nstep, Nsoil_nodes, Ncanopy_nodes):
    """
    Creates temporary results dictionary to accumulate simulation results
    """

    results = {}

    for var in variables:

        var_name = var[0]
        dimensions = var[2]

        if 'canopy' in dimensions:
            if 'date' in dimensions:
                var_shape = [Nstep, Ncanopy_nodes]
            else:
                var_shape = [Ncanopy_nodes]

        elif 'soil' in dimensions:
            if 'date' in dimensions:
                var_shape = [Nstep, Nsoil_nodes]
            else:
                var_shape = [Nsoil_nodes]

        else:
            var_shape = [Nstep]

        results[var_name] = np.full(var_shape, np.NAN)

    return results


def _append_results(group, step, step_results, results):
    """
    Adds results from each simulation steps to temporary results dictionary
    """

    results_keys = results.keys()
    step_results_keys = step_results.keys()

    for key in step_results_keys:
        variable = group + '_' + key

        if variable in results_keys:
            if variable == 'z':
                results[variable] = step_results[key]
            else:
                results[variable][step] = step_results[key]

    return results

def _write_ncf(nsim=None, results=None, ncf=None):
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

    print("Writing results of simulation number: {} is finished".format(nsim))

def initialize_netcdf(variables, sim, soil_nodes, canopy_nodes, forcing, filepath='results', filename='climoss.nc'):
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

#    filepath = os.path.join(climoss_path, filepath)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ff = os.path.join(filepath, filename)

    # create dataset and dimensions
    ncf = Dataset(ff, 'w')
    ncf.description = 'CliMoss model run results.'
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'APES_Jan2016'

    ncf.createDimension('date', date_dimension)
    ncf.createDimension('simulation', simulation_dimension)
    ncf.createDimension('soil', soil_dimension)
    ncf.createDimension('canopy', canopy_dimension)

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

#    print("netCDF4 path: " + ff)
    return ncf, ff