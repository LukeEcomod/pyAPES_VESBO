#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 08:54:17 2018

@author: ajkieloaho
"""

import numpy as np

EPS = np.finfo(float).eps  # machine epsilon

class DegreeDaySnow(object):
    """ Represents snow cover-soil-atmosphere interactions with simple
     zero-dimensional degree-day model
    """

    def __init__(self, properties):
        """
        snowpack = {
        'kmelt': 2.31e-8,  # Melting coefficient [m degC-1 s-1] (=2.0 mm/C/d)
        'kfreeze': 5.79e-9,  # Freezing  coefficient [m degC-1 s-1] (=0.5 mm/C/d)
        'retention': 0.2,  # max fraction of liquid water in snow [-]
        'Tmelt': 0.0,  # temperature when melting starts [degC]
        'swe_ini': 0.0,  # initial snow water equivalent [m],
        'optical_properties': {
                'albedo': {'PAR': 0.8, 'NIR': 0.8}
                'emissivity': 0.97,
                'albedo_PAR': 0.8,
                'albedo_NIR': 0.8,
                }
        }

        """

        #self.properties = properties

        # melting and freezing coefficients [kg m-2 s-1]
        self.kmelt = properties['kmelt']
        self.kfreeze = properties['kfreeze']

        # max fraction of liquid water in snow [-]
        self.retention = properties['retention']
        self.Tmelt = properties['Tmelt']

        self.optical_properties = properties['optical_properties']

        # state variables:
        self.temperature = properties['initial_conditions']['temperature']
        self.swe = properties['initial_conditions']['snow_water_equivalent']  # [kg m-2]
        self.ice = properties['initial_conditions']['snow_water_equivalent'] # ice content
        self.liq = 0.0  # liquid water storage in snowpack [kg m-2]

        # temporary storage of iteration results
        self.iteration_state = None

    def update(self):
        """ Updates new states to the snowpack.
        """
        self.temperature = self.iteration_state['temperature']
        self.ice = self.iteration_state['ice']
        self.liq = self.iteration_state['liq']
        self.swe = self.iteration_state['swe']

    def run(self, dt, forcing):
        """ Calculates one timestep and updates states of SnowpackModel instance

        Args:
            dt: timestep [s]
            forcing (dict):
                'air_temperature': [degC]
                'precipitation_rain': [kg m-2 s-1]
                'precipitation_snow': [kg m-2 s-1]

        Returns:
            fluxes (dict):
                'throughfall': [kg m-2 s-1]
                'water_closure': [kg m-2 s-1]
            states (dict):
                'snow_water_equivalent': [kg m-2 s-1]
        """

        """ --- initial conditions for calculating mass balance error --"""


        """ --- melting and freezing in snopack --- """
        if forcing['air_temperature'] >= self.Tmelt:
            # [m]
            melt = np.minimum(self.ice,
                              self.kmelt * dt * (forcing['air_temperature'] - self.Tmelt))
            freeze = 0.0

        else:
            melt = 0.0
            freeze = np.minimum(self.liq,
                                self.kfreeze * dt * (self.Tmelt - forcing['air_temperature']))

        """ --- update state of snowpack and compute potential infiltration --- """
        ice = np.maximum(0.0,
                         self.ice + forcing['precipitation_snow'] * dt + freeze - melt)

        liq = np.maximum(0.0,
                         self.liq + forcing['precipitation_rain'] * dt - freeze + melt)

        pot_inf = np.maximum(0.0, liq - ice * self.retention)

        # liquid water and ice in snow, and snow water equivalent [m]
        liq = np.maximum(0.0, liq - pot_inf)
        ice = ice
        swe = liq + ice

        # mass-balance error [kg m-2]
        water_closure = ((swe - self.swe)
                         - (forcing['precipitation_rain'] * dt + forcing['precipitation_snow'] * dt - pot_inf))

        # store iteration state
        self.iteration_state =  {'temperature': forcing['air_temperature'],
                                 'swe': swe,
                                 'ice': ice,
                                 'liq': liq}

        fluxes = {'potential_infiltration': pot_inf / dt,
                  'water_closure': water_closure / dt
                 }

        states = {'snow_water_equivalent': swe,
                  'temperature': forcing['air_temperature']
                 }

        return fluxes, states

# EOF

