# -*- coding: utf-8 -*-
"""
.. module: pyAPES
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Model framework for soil-bottom layer-canopy-atmosphere interactions.
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Oct 02 09:04:05 2018

References:
Launiainen, S., Katul, G.G., Lauren, A. and Kolari, P., 2015. Coupling boreal
forest CO2, H2O and energy flows by a vertically structured forest canopy – 
Soil model with separate bryophyte layer. Ecological modelling, 312, pp.385-405.

To call model and read results:
    from tools.iotools import read_results
    outputfile = driver(create_ncf=True, dbhfile="letto2014.txt")
    results = read_results(outputfile)
"""

import os
import numpy as np
import pandas as pd
from tools.iotools import read_forcing
from canopy.canopy import CanopyModel
from soilprofile.soilprofile import SoilModel
from parameters.canopy import get_cpara

from copy import deepcopy

def driver(create_ncf=False, dbhfile="letto2014.txt"):
    """
    """

    # Import general parameters
    from parameters.general import gpara
    # Import canopy model parameters
    cpara = get_cpara(dbhfile)
    # Import soil model parameters
    from parameters.soil import spara

    Nsim = 1

    # Read forcing
    forcing = read_forcing(gpara['forc_filename'],
                           gpara['start_time'],
                           gpara['end_time'],
                           dt=gpara['dt'])

    tasks = []

    for k in range(Nsim):
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
                tasks[k].Nplant_types,
                forcing,
                filename=filename,
                description=dbhfile)

        for task in tasks:
            print('Running simulation number: {}' .format(task.Nsim))
            running_time = time.time()
            results = task.run()
            print('Running time %.2f seconds' % (time.time() - running_time))
            _write_ncf(nsim=task.Nsim, results=results, ncf=ncf)

            del results

        output_file = "results/" + filename
        print('Ready! Results are in: ' + output_file)
        ncf.close()

    else:
        results = {task.Nsim: task.run() for task in tasks}
        output_file = results

    return output_file

class Model(object):
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
        self.Ncanopy_nodes = canopy_para['grid']['Nlayers']

        # create soil model instance
        self.soil_model = SoilModel(soil_para['z'], soil_para)

        # create canopy model instance
        self.canopy_model = CanopyModel(canopy_para, self.soil_model.grid['dz'])

        self.Nplant_types = len(self.canopy_model.Ptypes)

        self.results = _create_results(gen_para['variables'],
                                       self.Nsteps,
                                       self.Nsoil_nodes,
                                       self.Ncanopy_nodes,
                                       self.Nplant_types)

    def run(self):
        """ Runs atmosphere-canopy-soil--continuum model"""
        k_steps=np.arange(0, self.Nsteps, int(self.Nsteps/10))
        for k in range(0, self.Nsteps):
            # print str(k)
            if k in k_steps[:-1]:
                s = str(np.where(k_steps==k)[0][0]*10) + '%'
                print '{0}..\r'.format(s),

# Check Tsoil
            # Soil moisture forcing for canopy model
            Tsoil = self.forcing['Tair'].iloc[k]  # should come from soil model!
            Wsoil = self.soil_model.Wliq[0]  # certain depth?!

            """ Canopy, moss and Snow """
            # run daily loop (phenology and seasonal LAI)
            if self.forcing['doy'].iloc[k] != self.forcing['doy'].iloc[k-1] or k == 0:
                self.canopy_model.run_daily(
                        self.forcing['doy'].iloc[k],
                        self.forcing['Tdaily'].iloc[k])
            # properties of first soil node
            canopy_forcing = {'depth': self.soil_model.grid['z'][0],
                              'soil_temperature': self.soil_model.T[0],
                              'soil_water_potential': self.soil_model.h[0],
                              'soil_volumetric_water': self.soil_model.Wliq[0],
                              'soil_pond_storage': self.soil_model.pond,
                              'soil_hydraulic_conductivity': self.soil_model.Kv[0],
                              'soil_thermal_conductivity': self.soil_model.Lambda[0],
                              'wind_speed': self.forcing['U'].iloc[k],
                              'air_temperature': self.forcing['Tair'].iloc[k],
                              'precipitation': self.forcing['Prec'].iloc[k],
                              'h2o': self.forcing['H2O'].iloc[k],
                              'co2': self.forcing['CO2'].iloc[k],
                              'PAR': {'direct': self.forcing['dirPar'].iloc[k],
                                      'diffuse': self.forcing['diffPar'].iloc[k]},
                              'NIR': {'direct': self.forcing['dirNir'].iloc[k],
                                      'diffuse': self.forcing['diffNir'].iloc[k]},
                              'lw_in': self.forcing['LWin'].iloc[k],
                              'friction_velocity': self.forcing['Ustar'].iloc[k],
                              'air_pressure': self.forcing['P'].iloc[k],
                              'zenith_angle': self.forcing['Zen'].iloc[k]
                              }

            # run timestep loop
            canopy_flux, canopy_state, ffloor_flux, ffloor_state = self.canopy_model.run_timestep(
                    dt=self.dt,
                    forcing=canopy_forcing)

            """ Water and Heat in soil """
            # potential infiltration and evaporation from ground surface
            ubc_w = {'Prec': ffloor_flux['potential_infiltration'],
                     'Evap': ffloor_flux['soil_evaporation']}

            # transpiration sink
            rootsink = np.zeros(self.soil_model.Nlayers)
            rootsink[0:len(self.canopy_model.rad)] = self.canopy_model.rad * canopy_flux['transpiration']
            rootsink = rootsink / self.soil_model.grid['dz']

            ubc_T = {'type': 'flux', 'value': -ffloor_flux['ground_heat_flux']}

            # run soil water and heat flow
            soil_flux, soil_state = self.soil_model._run(
                    self.dt,
                    ubc_w, ubc_T,
                    water_sink=rootsink)

            forcing_state = {
                    'wind_speed': self.forcing['U'].iloc[k],
                    'air_temperature': self.forcing['Tair'].iloc[k],
                    'precipitation': self.forcing['Prec'].iloc[k],
                    'h2o': self.forcing['H2O'].iloc[k],
                    'co2': self.forcing['CO2'].iloc[k]}

            canopy_state.update(canopy_flux)
            ffloor_state.update(ffloor_flux)
            soil_state.update(soil_flux)

            self.results = _append_results('forcing', k, forcing_state, self.results)
            self.results = _append_results('canopy', k, canopy_state, self.results)
            self.results = _append_results('ffloor', k, ffloor_state, self.results)
            self.results = _append_results('soil', k, soil_state, self.results)

        print '100%'
        self.results = _append_results('canopy', None, {'z': self.canopy_model.z}, self.results)
        self.results = _append_results('soil', None, {'z': self.soil_model.grid['z']}, self.results)

        return self.results

def _create_results(variables, Nstep, Nsoil_nodes, Ncanopy_nodes, Nplant_types):
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

        elif 'planttype' in dimensions:
            if 'date' in dimensions:
                var_shape = [Nstep, Nplant_types]
            else:
                var_shape = [Nplant_types]


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

#        print(variable)

        if variable in results_keys:
            if variable == 'z':
                results[variable] = step_results[key]
            else:
#                print(type(step_results[key]))
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

def initialize_netcdf(variables, sim, soil_nodes, canopy_nodes, plant_nodes, forcing, filepath='results', filename='climoss.nc',
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

#    filepath = os.path.join(climoss_path, filepath)

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

#    print("netCDF4 path: " + ff)
    return ncf, ff