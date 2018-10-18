#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 15:31:23 2018

@author: ajkieloaho
"""

import numpy as np
from canopy.constants import EPS, MOLAR_MASS_H2O
from heat_and_water import saturation_vapor_pressure

def water_exchage(self, dt, forcing, properties):
    """

    Args:
        dt: [s]
        forcing (dict):
            'throughfall_rain': [mm]
            'windspeed': [m s-1]
            'air_temperature': [degC]
            'h2o': [mol mol-1]
            'air_pressure': [Pa]
    """

    U
    T
    H2O

    es, _ = saturation_vapor_pressure(forcing['air_temperature'])

    # [mol mol-1]
    D = np.maximum(0.0, es / forcing['air_pressure'] - forcing['h2o'])

    # initial water content
    Wo = self.Wold

    # interception and throughfall rate, new storage
    interception = np.maximum(0.0,
                              np.minimum(forcing['throughfall_rain'],
                                         self.Wmax - Wo))

    # [mm]
    potential_infiltration = forcing['throughfall_rain'] - interception

    # intermediate storage
    # [mm]
    W = Wo + interception

    # evaporation from moss layer: actual conductance is boundary layer x
    # correction for internal resistance
    grel = np.minimum(0.1285 * W / self.Mdry - 0.1285, 1.0)
    gb = grel * self._boundary_layer_conductance(U)

    erate = gb * D  # mol m-2 s-1
    # rate = 1.26*eq_evap(Rn, T, units='mol')  # unrestricted rate
    Evap = np.minimum(erate * MOLAR_MASS_H2O * dt, W - self.Wmin)  # mm
    self.W = W - Evap  # mm

    Mbe = (self.W - Wo) - (forcing['throughfall_rain'] - Evap - potential_infiltration)
    # print('Mbe', Mbe)

    return Evap/dt, Trfall, Mbe