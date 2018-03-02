# -*- coding: utf-8 -*-
"""
Functions and classes for snowpack descriptions

@author: L1656
"""

import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from evapotranspiration import penman_monteith

class Snowpack():
    """
    Computes snowpack state and potential infiltration to soil
    Now simple degreeday model
    ! Add option for energy balance based description of two-layer snowpack?
    """
    def __init__(self, p, cf):
        """
        Args:
            dt: timestep [s]
            p parameter dict
               'kmelt':  melting coefficient [m degC-1 s-1]
               'kfreeze':  freezing  coefficient [m degC-1 s-1]
               'retention':  max fraction of liquid water in snow [-]
               'Tmelt': temperature when melting starts [degC]
               'swe_ini':  initial snow water equivalent [m]
            cf: canopy closure [-]
        Returns:
            Snowpack -object
        """

        # parameters:
        # melting and freezing coefficients [m/s]
        self.kmelt = p['kmelt'] - 1.64 / (24 * 3600 * 1000) * cf  # Kuusisto E, 'Lumi Suomessa'
        self.kfreeze = p['kfreeze']
        # max fraction of liquid water in snow [-]
        self.reten = p['retention']
        self.Tmelt = p['Tmelt']

        # state variables:
        self.swe_ice = p['swe_ini']  # water equivalent of ice in snowpack [m]
        self.swe_liq = 0.0  # liquid water storage in snowpack [m]
        self.swe = self.swe_liq + self.swe_ice  # snow water equivalent [m]

    def _run(self, dt, T, Trfall_rain, Trfall_snow):
        """
        Args:
            dt: timestep [s]
            T: air temperature [degC]
            Trfall_rain: throughfall as rainfall [m]
            Trfall_snow: throughfall as snowfall [m]
        Returns:
            updates sate self.swe_ice, self.swe_liq
            swe: snow water equivalent [m]
            Potinf: potential infiltration [m]
            MBE: mass balance error
        """

        """ --- initial conditions for calculating mass balance error --"""
        swe = self.swe  # initial state m

        """ --- melting and freezing in snopack --- """
        if T >= self.Tmelt:
            Melt = np.minimum(self.swe_ice, self.kmelt * dt * (T - self.Tmelt))  # m
            Freeze = 0.0
        else:
            Melt = 0.0
            Freeze = np.minimum(self.swe_liq, self.kfreeze * dt * (self.Tmelt - T))  # m

        """ --- update state of snowpack and compute potential infiltration --- """
        swei = np.maximum(0.0, self.swe_ice + Trfall_snow + Freeze - Melt)
        swel = np.maximum(0.0, self.swe_liq + Trfall_rain - Freeze + Melt)
        # potential infiltration [m]
        PotInf = np.maximum(0.0, swel - swei * self.reten)
        # liquid water and ice in snow, and snow water equivalent [m]
        self.swe_liq = np.maximum(0.0, swel - PotInf)
        self.swe_ice = swei
        self.swe = self.swe_liq + self.swe_ice

        # mass-balance error mm
        MBE = (self.swe - swe) - (Trfall_rain + Trfall_snow - PotInf)

        return PotInf, MBE
