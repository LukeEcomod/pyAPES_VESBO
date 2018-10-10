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

#: [J kg\ :sup:`-1`\ ], latent heat of freezing
LATENT_HEAT_FREEZING = 333700.0
#: [\ :math:`^{\circ}`\ C], freezing point of water
FREEZING_POINT_H2O = 0.0
#: [kg m\ :sup:`-3`\ ], densities
ICE_DENSITY = 917.0

#: [J m\ :sup:`-3`\ K \ :sup:`-1`\], thermal condutivities
K_WATER = 0.57
K_ICE = 2.2
K_AIR = 0.025
K_ORG = 0.25
K_SAND = 7.7  # Tian et al. 2007
K_SILT = 2.74  # Tian et al. 2007
K_CLAY = 1.93  # Tian et al. 2007

#: volumetric heat capacieties  [J m\ :sup:`-3`\ K \ :sup:`-1`\]
CV_AIR = 1297.0  # air at 101kPa
CV_WATER = 4.18e6  # water
CV_ICE = 1.93e6  # ice
CV_ORGANIC = 2.50e6  # dry organic matter
CV_MINERAL = 2.31e6  # soil minerals

