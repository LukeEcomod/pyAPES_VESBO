#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module: bryophyte
    :synopsis: APES-model component
.. moduleauthor:: Antti-Jussi Kieloaho

Created on Fri Oct 26 08:37:46 2018


Note:
    - set soil_hydraulic_conductivity to ZERO

@author: ajkieloaho

"""

from copy import deepcopy

import numpy as np

from .heat_and_water import heat_and_water_exchange, evaporation_through
from .heat_and_water import convert_hydraulic_parameters
from .carbon import soil_respiration

from canopy.constants import WATER_DENSITY, MOLAR_MASS_H2O, MOLAR_MASS_C, EPS


class Litter(object):
    """
    """

    def __init__(self, properties, initial_conditions=None):
        """
        """
        dry_mass = properties['bulk_density'] * properties['height']
        properties['dry_mass'] = dry_mass

        residual_water_content = (properties['min_water_content']
                                  / WATER_DENSITY
                                  * properties['bulk_density'])

        field_capacity = (properties['max_water_content']
                          / WATER_DENSITY
                          * properties['bulk_density'])

        saturated_water_content = properties['porosity']

        if 'water_retention' not in properties:
            water_retention = {}

            water_retention['theta_r'] = residual_water_content
            water_retention['theta_s'] = saturated_water_content
            water_retention['field_capacity'] = field_capacity

            fraction = (
                    (saturated_water_content - field_capacity)
                    / (saturated_water_content - residual_water_content))

            water_retention['fraction'] = fraction
            water_retention['saturated_conductivity'] = properties['saturated_conductivity']

            water_retention['alpha'] = properties['alpha']
            water_retention['n'] = properties['n']

            if 'pore_connectivity' in properties:
                water_retention['pore_connectivity'] = properties['pore_connectivity']

            if 'compressibility' in properties:
                water_retention['compressability'] = properties['compressability']

            properties['water_retention'] = water_retention

        else:
            water_retention = properties['water_retention']
            water_retention['field_capacity'] = field_capacity

        self.porosity = properties['water_retention']['theta_s']

        self.properties = properties

        # set initial conditions
        if initial_conditions is not None:
            #: [:math:`^{\circ}`\ C]
            self.temperature = initial_conditions['temperature']

            #: [g g\ :sup:`-1`\ ]
            if initial_conditions['water_content'] <= properties['max_water_content']:
                self.water_content = initial_conditions['water_content']

            else:
                self.water_content = properties['max_water_content']

        else:
            #: [:math:`^{\circ}`\ C]
            self.temperature = 10.

            #: [g g\ :sup:`-1`\ ]
            self.water_content = (
                    properties['max_water_content']
                    + properties['min_water_content']) / 2.0

        self.water_storage = self.water_content * properties['dry_mass']

        #: [m\ :sup:`3` m\ :sup:`-3`\ ]
        self.volumetric_water = (
            self.water_content / WATER_DENSITY * properties['bulk_density'])

        self.volumetric_air = self.porosity - self.volumetric_water

        #: [m]
        self.water_potential = convert_hydraulic_parameters(
                self.volumetric_water,
                self.properties['water_retention'],
                'volumetric_water')

        #: [kg C m-2], 'free carbon pool'
        self.carbon_pool = 0.0
        self.coverage = properties['ground_coverage']

        self.old_carbon_pool = self.carbon_pool
        self.old_water_content = self.water_content
        self.old_water_storage = self.water_storage
        self.old_volumetric_water = self.volumetric_water
        self.old_volumetric_air = self.volumetric_air
        self.old_water_potential = self.water_potential
        self.old_temperature = self.temperature

    def update(self):
        """ Updates old states to states after iteration.
        """

        self.old_carbon_pool = self.carbon_pool
        self.old_water_content = self.water_content
        self.old_water_storage = self.water_storage
        self.old_volumetric_water = self.volumetric_water
        self.old_volumetric_air = self.volumetric_air
        self.old_water_potential = self.water_potential
        self.old_temperature = self.temperature

    def restore(self):
        """ Restores new states back to states before iteration.
        """

        self.carbon_pool = self.old_carbon_pool
        self.water_content = self.old_water_content
        self.water_storage = self.old_water_storage
        self.volumetric_water = self.old_volumetric_water
        self.volumetric_air = self.old_volumetric_air
        self.water_potential = self.old_water_potential
        self.temperature = self.old_temperature

    def run(self, dt, forcing, solver=False):
        r""" Calculates one timestep and updates states of Litter instance.
        Args:
            dt: timestep [s]
            forcing (dict): states of microclimate
                'throughfall': [mm s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ]
                'lwdn': [W m\ :sup:`-2`\ ]
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'precipitation_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_depth': [m]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [m]
                'soil_hydraulic_conductivity': [m s\ :sup:`-1`\ ]
                'soil_thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ]
                'nsteps' number of steps in odesolver

        Returns:
            fluxes (dict)
            states (dict)
        """

        # calculate moss energy and water balance and new state
        litter_forcing = deepcopy(forcing)

        # No capillar connection assumed
        litter_forcing['soil_hydraulic_conductivity'] = 0.0
        fluxes, states = heat_and_water_exchange(self.properties,
                                                 self.old_temperature,
                                                 self.old_water_content,
                                                 dt,
                                                 litter_forcing)

        # update state variables
        self.temperature = states['temperature']
        self.water_content = states['water_content']
        self.water_storage = states['water_storage']
        self.volumetric_water = states['volumetric_water']
        self.volumetric_air = self.porosity - states['volumetric_water']
        self.water_potential = states['water_potential']


        # calculate respiration
        respiration = soil_respiration(self.properties,
                                self.temperature,
                                self.volumetric_water,
                                self.volumetric_air
                                )

        fluxes.update({'respiration_rate': respiration})

        self.carbon_pool = self.old_carbon_pool + 1e3 * MOLAR_MASS_C * 1e-6 * respiration

        states['carbon_pool'] = self.carbon_pool

        # compute soil evaporation through moss layer

        # [mol mol-1] -> [Pa]
        h2o = forcing['h2o'] * forcing['air_pressure']

        # [mol m-2 s-1]
        soil_evaporation = evaporation_through(
            self.properties,
            self.volumetric_water,  # old
            self.temperature,  # old
            forcing['air_temperature'],
            h2o,
            forcing['wind_speed'],
            forcing['air_pressure'],
            forcing['soil_temperature'],
            forcing['soil_water_potential'],
            forcing['soil_hydraulic_conductivity'],
            forcing['depth'])

        # unit conversion: 1000 kg m-2 s-1 = mm s-1

        soil_evaporation = {key: value * MOLAR_MASS_H2O for key, value in soil_evaporation.items()}

        fluxes.update(soil_evaporation)

        return fluxes, states

# EOF
