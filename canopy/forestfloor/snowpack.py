#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 08:54:17 2018

@author: ajkieloaho
"""

import numpy as np

EPS = np.finfo(float).eps  # machine epsilon

class Snowpack(object):
    """ Represents snow cover-soil-atmosphere interactions
    """

    def __init__(self, properties, initital_conditions):
        """
        """

        self.properties = properties

        # melting and freezing coefficients [m/s]
        self.kmelt = properties['kmelt']
        self.kfreeze = properties['kfreeze']
        # max fraction of liquid water in snow [-]
        self.reten = properties['retention']
        self.Tmelt = properties['Tmelt']

        # state variables:
        self.swe_ice = properties['swe_ini']  # water equivalent of ice in snowpack [m]
        self.swe_liq = 0.0  # liquid water storage in snowpack [m]
        self.swe = self.swe_liq + self.swe_ice  # snow water equivalent [m]

        self.old_ice = self.swe_ice
        self.old_liq = self.swe_liq
        self.old_swe = self.swe


    def update(self):
        """
        """

        self.old_swe_ice = self.swe_ice
        self.old_swe_liq = self.swe_liq
        self.old_swe =  self.swe

    def restore(self):
        """
        """

    def run(self, dt, forcing):
        """ Calculates one timestep and updates states of SnowpackModel instance

        Args:
            dt: timestep [s]
            T: air temperature [degC]
            Trfall_rain: throughfall as rainfall [m]
            Trfall_snow: throughfall as snowfall [m]

        Returns:
            Potinf: potential infiltration [m]
            MBE: mass balance error [m]
        """

        """ --- initial conditions for calculating mass balance error --"""
        self.swe = self.old_swe  # initial state m

        """ --- melting and freezing in snopack --- """
        if forcing['air_temperature'] >= self.Tmelt:
            # [m]
            Melt = np.minimum(self.old_swe_ice,
                              self.kmelt * dt * (forcing['air_temperature'] - self.Tmelt))
            Freeze = 0.0

        else:
            # [m]
            Melt = 0.0
            Freeze = np.minimum(self.old_swe_liq,
                                self.kfreeze * dt * (self.Tmelt - forcing['air_temperature']))  # m

        """ --- update state of snowpack and compute potential infiltration --- """
        swe_ice = np.maximum(0.0,
                             self.oldswe_ice + Trfall_snow + Freeze - Melt)

        swe_liq = np.maximum(0.0,
                             self.oldswe_liq + Trfall_rain - Freeze + Melt)

        # potential infiltration [m]
        pot_inf = np.maximum(0.0,
                            swe_liq - swe_ice * self.reten)

        # liquid water and ice in snow, and snow water equivalent [m]
        self.swe_liq = np.maximum(0.0,
                                  swe_liq - PotInf)
        self.swe_ice = swe_ice
        self.swe = self.swe_liq + self.swe_ice

        # mass-balance error [m]
        water_closure = ((self.swe - self.old_swe)
                         - (Trfall_rain + Trfall_snow - pot_inf))


        return {'potential_infiltration': pot_inf,
                'water_closure': water_closure}

# EOF

