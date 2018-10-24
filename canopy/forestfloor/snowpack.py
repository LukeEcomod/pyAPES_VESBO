#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 08:54:17 2018

Note:
    migrated to python3
    - no changes

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
        self.retention = properties['retention']
        self.Tmelt = properties['Tmelt']

        # state variables:
        self.ice = properties['swe_ini']  # water equivalent of ice in snowpack [m]
        self.liq = 0.0  # liquid water storage in snowpack [m]
        self.swe = self.liq + self.ice  # snow water equivalent [m]

        self.old_ice = self.ice
        self.old_liq = self.liq
        self.old_swe = self.swe

    def snowcover(self):
        """ Returns true if there is snow cover else false
        """
        return self.swe > 0.0


    def update(self):
        """ Updates new states to the snowpack.
        """

        self.old_ice = self.ice
        self.old_liq = self.liq
        self.old_swe = self.swe

    def restore(self):
        """ Restores new states back to states before iteration.
        """

        self.ice = self.old_ice
        self.liq = self.old_liq
        self.swe = self.old_swe


    def run(self, dt, forcing):
        """ Calculates one timestep and updates states of SnowpackModel instance

        Args:
            dt: timestep [s]
            T: air temperature [degC]
            Trfall_rain: throughfall as rainfall [m s-1]
            Trfall_snow: throughfall as snowfall [m s-1]

        Returns:
            Potinf: potential infiltration [m s-1]
            MBE: mass balance error [m s-1]
        """

        """ --- initial conditions for calculating mass balance error --"""
        self.swe = self.old_swe  # initial state m

        """ --- melting and freezing in snopack --- """
        if forcing['air_temperature'] >= self.Tmelt:
            # [m]
            melt = np.minimum(self.old_ice,
                              self.kmelt * dt * (forcing['air_temperature'] - self.Tmelt))
            freeze = 0.0

        else:
            # [m]
            melt = 0.0
            freeze = np.minimum(self.old_liq,
                                self.kfreeze * dt * (self.Tmelt - forcing['air_temperature']))  # m

        """ --- update state of snowpack and compute potential infiltration --- """
        ice = np.maximum(0.0,
                             self.old_ice + forcing['throughfall_snow'] * dt + freeze - melt)

        liq = np.maximum(0.0,
                             self.old_liq + forcing['throughfall_rain'] * dt - freeze + melt)

        # potential infiltration [m]
        pot_inf = np.maximum(0.0,
                            liq - ice * self.retention)

        # liquid water and ice in snow, and snow water equivalent [m]
        self.liq = np.maximum(0.0,
                                  liq - pot_inf)
        self.ice = ice
        self.swe = self.liq + self.ice

        # mass-balance error [m]
        water_closure = ((self.swe - self.old_swe)
                         - (forcing['throughfall_rain'] * dt + forcing['throughfall_snow'] * dt - pot_inf))


        fluxes = {'potential_infiltration': pot_inf / dt,
                  'water_closure': water_closure / dt}

        states = {'snow_water_equivalent': self.swe}

        return fluxes, states

# EOF

