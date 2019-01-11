# -*- coding: utf-8 -*-
"""
.. module: soil
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Note:
    migrated to python3
    - absolute imports
    - dict.keys are not wrapped in list()

Soilprofile

Created on Tue Oct 02 09:04:05 2018
"""

import numpy as np
import pandas as pd
from datetime import datetime
from matplotlib import pyplot as plt
from .constants import EPS

from .water import Water
from .heat import Heat

class Soil(object):

    def __init__(self, p):
        r""" Initializes soilprofile for 1D solutions of water and heat balance

        Args:
            p (dict):
                'grid' (dict):
                    'zh' (list/array?): bottom elevation of layers with different parametrizations [m], soil surface = 0.0
                    'dz' (list/array): thickness of computational layers [m]
                'profile_properties' (dict): all parameters given as list/array of len(zh)
                    'pF' (dict): water retention parameters (van Genuchten)
                        'ThetaS' (list/array): saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                        'ThetaR' (list/array): residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                        'alpha' (list/array): air entry suction [cm\ :sup:`-1`]
                        'n' (list/array): pore size distribution [-]
                    'saturated_conductivity_vertical' (list/array): [m s-1]
                    'saturated_conductivity_horizontal' (list/array): [m s-1]
                    'dry_heat_capacity'  (list/array or None): [J m-3 (total volume) K-1]
                        ! if None, estimated from organic/mineral composition
                    'solid_composition' (dict): volume fractions of solid volume [-]
                        'organic' (list/array)
                        'sand' (list/array)
                        'silt' (list/array)
                        'clay' (list/array)
                    'freezing_curve' (list/array): freezing curve parameter [-]
                    'bedrock' (dict): 
                        'heat_capacity' (float): [J m-3 (total volume) K-1]
                        'thermal_conductivity' (float): [W m-1 K-1]
                'water_model':
                    'solve': True/False
                    'type': 'Richards',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
                    'pond_storage_max': 0.01,  #  maximum pond depth [m]
                    'initial_condition' (dict): initial conditions
                        'ground_water_level' (float): groundwater depth [m]
                        'pond_storage' (float): pond depth at surface [m]
                    'lower_boundary (dict)': lower boundary condition
                                  'type' (str): 'impermeable', 'flux', 'free_drain' or 'head'
                                  'value' (float or None): give for 'head' [m] and 'flux' [m s-1]
                                  'depth' (float): depth of impermeable boundary (=bedrock) [m]
                    ('drainage_equation' (dict): drainage equation and drainage parameters
                                      'type' (str): 'hooghoudt' (no other options yet)
                                      'depth' (float): drain depth [m]
                                      'spacing' (float): drain spacing [m]
                                      'width' (float):  drain width [m])
                'heat_model':
                    'solve': True/False
                    'initial_condition':
                        'temperature' (array/float): initial temperature [degC]
                    'lower_boundary':
                        'type' (str): 'flux' or 'temperature'
                        'value' (float): value of flux [W m-2] or temperature [degC]

        Returns:
            self (object)
        """

        self.solve_heat = p['heat_model']['solve']
        self.solve_water = p['water_model']['solve']

        # grid
        dz = np.array(p['grid']['dz'])
        z = dz / 2 - np.cumsum(dz)  # depth of nodes from soil surface
        N = len(dz)
        self.Nlayers = N

        dzu = np.zeros(N)  # distances between grid points i-1 and i
        dzu[1:] = z[:-1] - z[1:]
        dzu[0] = dz[0] / 2.0  # from soil surface to first node

        dzl = np.zeros(N)  # distances between grid points i and i+1
        dzl[:-1] = z[:-1] - z[1:]
        dzl[-1] = dz[-1] / 2.0 #  from last node to bottom surface

        self.grid = {'z': z,
                     'dz': dz,
                     'dzu': dzu,
                     'dzl': dzl}

        self.ones = np.ones(N)

        profile_properties = form_profile(z, p['grid']['zh'], p['soil_properties'],
                                          p['water_model']['lower_boundary'])

        # initialize water and heat intance
        self.water = Water(self.grid, profile_properties, p['water_model'])

        self.heat = Heat(self.grid, profile_properties, p['heat_model'],
                                   self._fill(self.water.Wtot, 0.0))

    def run(self, dt, forcing, water_sink=None, heat_sink=None, lbc_water=None, lbc_heat=None, state=None):
        r""" Runs soil model for one timestep.
        Args:
            dt: time step [s]
            forcing (dict):
                'potential_infiltration': [m s-1]
                'potential_evaporation': [m s-1]
                'atmospheric_pressure_head': [m]
                'ground_heat_flux' (float): heat flux from soil surface [W m-2]
                    OR 'temperature' (float): soil surface temperature [degC]
                'state_water': if water model is not solved this allows state can be changed over time
                    'ground_water_level' OR 'volumetric_water_content'
                'state_heat': if heat model is not solved this allows state can be changed over time
                    'temperature'
            water_sink (array or None): water sink from layers, e.g. root sink [m s-1]
                ! array length can be only root zone or whole soil profile
                ! if None set to zero
            heat_sink (array or None): heat sink from layers [??????]
            lbc_water (dict or None): allows boundary to be changed in time
                'type': 'impermeable', 'flux', 'free_drain', 'head'
                'value': give for 'head' [m] and 'flux' [m s-1]
            lbc_heat (dict or None): allows boundary to be changed in time
                '':
                '':

        Returns:
            fluxes (dict):[m s-1]
                'infiltration'
                'evaporation'
                'drainage'
                'transpiration'
                'surface_runoff'
                'water_closure'
            state
                'temperature' [degC]
                'water_potential' [m]
                'volumetric_water_content' [m3m-3]
                'ice_content' [m3m-3]
                'pond_storage' [m]
                'ground_water_level' [m]
                'hydraulic_conductivity' [ms-1]
                'thermal_conductivity [Wm-1s-1]
        """
        fluxes = {}

        if self.solve_water:
            water_fluxes = self.water.run(dt,
                                          forcing,
                                          water_sink=water_sink,
                                          lower_boundary=lbc_water)
            fluxes.update(water_fluxes)

        elif 'state_water' in forcing:
            self.water.update_state(forcing['state_water'])

        state = {'water_potential': self._fill(self.water.h),
                 'volumetric_water_content': self._fill(self.water.Wtot, 0.0),
                 'hydraulic_conductivity': self._fill(self.water.Kv),
                 'pond_storage': self.water.h_pond,
                 'ground_water_level': self.water.gwl}

        if self.solve_heat:
            heat_fluxes = self.heat.run(dt,
                                        forcing,
                                        state['volumetric_water_content'],
                                        heat_sink=heat_sink,
                                        lower_boundary=lbc_heat)

            fluxes.update(heat_fluxes)
        else:
            if'state_heat' in forcing:
                self.heat.update_state(state=forcing['state_heat'],
                                       Wtot=state['volumetric_water_content'])
            else:
                self.heat.update_state(Wtot=state['volumetric_water_content'])

        state.update({'volumetric_ice_content': self.heat.Wice,
                      'temperature': self.heat.T,
                      'thermal_conductivity': self.heat.thermal_conductivity})

        return fluxes, state

    def _fill(self, x, value=np.nan):
        r""" Fills arrays to entire soil profile.
        Args:
            x (array): array to change to len(z)
            value (float): value used to fill array
        Returns:
            x (array): input x extended to len(z)
        """
        if len(x) < len(self.ones):
            x_filled = self.ones * value
            x_filled[self.water.ix] = x
            return x_filled
        else:
            return x

def form_profile(z, zh, p, lbc_water):
    r""" Forms dictonary of soil profile parameters (all computational nodes)
    Args:
        z (array): node elevations, soil surface = 0.0 [m]
        zh (list): bottom elevation of layers with different parametrizations [m]
        p (dict): all parameters given as list/array of len(zh)
            'pF' (dict): water retention parameters (van Genuchten)
                'ThetaS' (list/array): saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                'ThetaR' (list/array): residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                'alpha' (list/array): air entry suction [cm\ :sup:`-1`]
                'n' (list/array): pore size distribution [-]
            'saturated_conductivity_vertical' (list/array): [m s-1]
            'saturated_conductivity_horizontal' (list/array): [m s-1]
            'dry_heat_capacity'  (list/array or None): [J m-3 (total volume) K-1]
                ! if None, estimated from organic/mineral composition
            'solid_composition' (dict): fractions of solid volume [-]
                'organic' (list/array)
                'sand' (list/array)
                'silt' (list/array)
                'clay' (list/array)
            'freezing_curve' (list/array): freezing curve parameter [-]
            'bedrock' (dict): 
                'heat_capacity' (float): [J m-3 (total volume) K-1]
                'thermal_conductivity' (float): [W m-1 K-1]
        lbc_water (dict):
            'type' (str): 'impermeable', 'flux', 'free_drain' or 'head'
            'depth' (float): depth of impermeable boundary (=bedrock) [m]
    Returns:
        prop (dict): all parameters given as arrays of len(z)
            'pF' (dict): water retention parameters (van Genuchten)
                'ThetaS' (array): saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                'ThetaR' (array): residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                'alpha' (array): air entry suction [cm\ :sup:`-1`]
                'n' (array): pore size distribution [-]
            'saturated_conductivity_vertical' (array): [m s-1]
            'saturated_conductivity_horizontal' (array): [m s-1]
            'porosity' (array): soil porosity (=ThetaS) [m\ :sup:`3` m\ :sup:`-3`\ ]
            'dry_heat_capacity' (array): [J m-3 (total volume) K-1]
                ! if nan, estimated from organic/mineral composition
            'solid_composition' (dict): fractions of solid volume [-]
                'organic' (array)
                'sand' (array)
                'silt' (array)
                'clay' (array)
            'freezing_curve' (array): freezing curve parameter [-]
            'bedrock_thermal_conductivity' (array): nan above bedrock depth[W m-1 K-1]
    """
    if p['solid_heat_capacity'] is None:
        p['solid_heat_capacity'] = [np.nan] * len(zh)
    N = len(z)
    prop = {}
    # horizons to computational nodes
    ix = np.zeros(N)
    for depth in zh:
        ix += np.where(z < depth, 1, 0)

    for key in p.keys():
        if (key == 'pF' or key == 'solid_composition'):
            prop.update({key: {}})

            for subkey in p[key].keys():
                pp = np.array([p[key][subkey][int(ix[i])] for i in range(N)])
                prop[key].update({subkey: pp})

        elif key != 'bedrock':
            #print(key)
            pp = np.array([p[key][int(ix[i])] for i in range(N)])
            prop.update({key: pp})

    prop.update({'porosity': prop['pF']['ThetaS']})
    prop.update({'bedrock_thermal_conductivity': np.ones(N) * np.nan})

    if lbc_water['type'] == 'impermeable' and z[-1] < lbc_water['depth']:
        ixx = np.where(z < lbc_water['depth'])[0]
        prop['porosity'][ixx] = EPS
        prop['solid_heat_capacity'][ixx] = p['bedrock']['solid_heat_capacity']
        prop['bedrock_thermal_conductivity'][ixx] = p['bedrock']['thermal_conductivity']

    return prop