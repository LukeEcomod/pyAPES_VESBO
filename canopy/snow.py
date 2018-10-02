# -*- coding: utf-8 -*-
"""
.. module: snow
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes temperature-based snow accumulation and melt.

Created on Thu Mar 01 13:21:29 2018

Note: implement energy balance approach!
"""

import numpy as np
eps = np.finfo(float).eps  # machine epsilon

class Snowpack(object):
    r"""Describes temperature-based snow accumulation and melt.
    """
    def __init__(self, p):
        """
        Args:
            p (dict):
               'kmelt':  melting coefficient [m degC-1 s-1]
               'kfreeze': freezing  coefficient [m degC-1 s-1]
               'retention': max fraction of liquid water in snow [-]
               'Tmelt': temperature when melting starts [degC]
               'swe_ini':  initial snow water equivalent [m]
        Returns:
            self (object)
                .swe_ice: water equivalent of ice in snowpack [m]
                .swe_liq: liquid water storage in snowpack [m]
                .swe: snow water equivalent [m]
        """

        # parameters:
        # melting and freezing coefficients [m/s]
        self.kmelt = p['kmelt']
        self.kfreeze = p['kfreeze']
        # max fraction of liquid water in snow [-]
        self.reten = p['retention']
        self.Tmelt = p['Tmelt']

        # state variables:
        self.swe_ice = p['swe_ini']  # water equivalent of ice in snowpack [m]
        self.swe_liq = 0.0  # liquid water storage in snowpack [m]
        self.swe = self.swe_liq + self.swe_ice  # snow water equivalent [m]

        self.update()

    def _run(self, dt, T, Trfall_rain, Trfall_snow):
        """
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
        self.swe = self.oldswe  # initial state m

        """ --- melting and freezing in snopack --- """
        if T >= self.Tmelt:
            Melt = np.minimum(self.oldswe_ice, self.kmelt * dt * (T - self.Tmelt))  # m
            Freeze = 0.0
        else:
            Melt = 0.0
            Freeze = np.minimum(self.oldswe_liq, self.kfreeze * dt * (self.Tmelt - T))  # m

        """ --- update state of snowpack and compute potential infiltration --- """
        swei = np.maximum(0.0, self.oldswe_ice + Trfall_snow + Freeze - Melt)
        swel = np.maximum(0.0, self.oldswe_liq + Trfall_rain - Freeze + Melt)
        # potential infiltration [m]
        PotInf = np.maximum(0.0, swel - swei * self.reten)
        # liquid water and ice in snow, and snow water equivalent [m]
        self.swe_liq = np.maximum(0.0, swel - PotInf)
        self.swe_ice = swei
        self.swe = self.swe_liq + self.swe_ice

        # mass-balance error [m]
        MBE = (self.swe - self.oldswe) - (Trfall_rain + Trfall_snow - PotInf)

        return PotInf, MBE

    def update(self):
        """Update snowpack state to old state
        """
        self.oldswe_ice = self.swe_ice
        self.oldswe_liq = self.swe_liq
        self.oldswe =  self.swe