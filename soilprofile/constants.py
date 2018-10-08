#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
.. module: constants
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Constants used in soilprofile module.

Created on Wed Oct 03 15:29:40 2018
"""

import numpy as np

EPS = np.finfo(float).eps  # machine epsilon

#: [m s\ :sup: `-2`\ ], standard gravity
GRAVITY = 9.81
#: [K], zero degrees celsius in Kelvin
DEG_TO_KELVIN = 273.15
#: [J kg\ :sup:`-1`\ ], latent heat of freezing
LATENT_HEAT_FREEZING = 333700.0
#: [\ :math:`^{\circ}`\ C], freezing point of water
FREEZING_POINT_H2O = 0.0

#: [kg m\ :sup:`-3`\ ], densities
ICE_DENSITY = 917.0
WATER_DENSITY = 1.0e3
MINERAL_DENSITY = 2.65e3
ORGANIC_DENSITY = 1.30e3
AIR_DENSITY = 1.25

#: [J m\ :sup:`-3`\ K \ :sup:`-1`\], thermal condutivities
K_WATER = 0.57
K_ICE = 2.2
K_AIR = 0.025
K_QUARTZ = 8.8
K_MINERAL = 2.9
K_ORG = 0.25
K_SAND = 7.7  # Tian et al. 2007
K_SILT = 2.74  # Tian et al. 2007
K_CLAY = 1.93  # Tian et al. 2007

#: [J kg\ :sup:`-1`\ K \ :sup:`-1`\], specific heats capacities
CP_AIR = 1297.0
CP_WATER = 4180.0
CP_ICE = 2100.0
CP_ORGANIC = 1920.0
CP_MINERAL = 870.0
CP_GRANITE = 820.0
CP_QUARTZ = 800.0

