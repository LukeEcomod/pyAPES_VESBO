# -*- coding: utf-8 -*-
"""
.. module: rootzone
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes roots of planttype.
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Jul 24 10:49:59 2018

Note:
    migrated to python3
    - nothing changed

References: --- SHOULD WE USE THIS??? Do we need to solve hroot???
Volpe, V., Marani, M., Albertson, J.D. and Katul, G., 2013. Root controls
on water redistribution and carbon uptake in the soil–plant system under
current and future climate. Advances in Water resources, 60, pp.110-120.
"""

import numpy as np
import matplotlib.pyplot as plt

from canopy.constants import EPS

class RootUptake(object):
    r""" Describes roots of planttype.
    """
    def __init__(self, p, dz_soil, LAImax):
        r""" Initializes rootuptake object.

        Args:
            p (dict):
                'root_depth': depth of rooting zone [m]
                'beta': shape parameter for root distribution model
                'RAI_LAI_multiplier': multiplier for total fine root area index (RAI = 2*LAImax)
                'fine_radius': fine root radius [m]
                'radial_K': maximum bulk root membrane conductance in radial direction [s-1]
            dz_soil (array): thickness of soilprofile layers from top to bottom [m]
            LAImax (float): maximum leaf area index [m2 m-2]
        Returns:
            self (object)
        """
        # parameters
        self.root_depth = p['root_depth']
        self.fine_radius = p['fine_radius']  # fine root radius [m]
        self.root_cond = p['root_cond']  # [s]
        self.RAI = p['RAI_LAI_multiplier']*LAImax  # total fine root area index (m2/m2)
        self.rad = self.RAI * RootDistribution(p['beta'], dz_soil, p['root_depth'])
        self.ix = np.where(np.isfinite(self.rad))
        self.dz = dz_soil[self.ix]

        # state variables
        self.h_root = 0.0

    def wateruptake(self, transpiration, h_soil, kh_soil):
        r""" Root wateruptake based on root water potential (Volpe et al 2013)

        Agrs:

        Returns:

        """

        # conductance from soil to root-soil interface [s-1]
        alpha = np.sqrt(self.root_depth/self.RAI) / np.sqrt(2.0 * self.fine_radius)
        ks = alpha * kh_soil[self.ix] * self.rad

        # conductance from soil-root interface to base of xylem [s-1]
        kr = self.rad * self.dz / self.root_cond

        # soil to xylem conductance [s-1]
        g_sr = ks * kr / (ks + kr + EPS)

        # assume total root uptake equals transpiration rate and solve uniform root pressure [m]
        self.h_root = -(transpiration - sum(g_sr * h_soil[self.ix])) / sum(g_sr)

        # root uptake [m s-1]
        rootsink = g_sr * (h_soil[self.ix] - self.h_root)

        return rootsink

def RootDistribution(beta, dz, root_depth):
    r""" Computes normalized root area density distribution with depth.
    (sum(dz*R) = 1)

    Args:
        beta: shape parameter for root distribution model
        dz (array):  thickness soil layers from top to bottom [m]
        root_depth: depth of rooting zone [m]
    Returns:
        R (array): normalized root area density distribution with depth,
            extends only to depth of rooting zone

    Reference:
        Gale and Grigal, 1987 Can. J. For.Res., 17, 829 - 834.
    """
    z = np.concatenate([[0.0], np.cumsum(dz)])
    z = np.concatenate([z[z < root_depth], [root_depth]])
    d = abs(z * 100.0)  # depth in cm

    Y = 1.0 - beta**d  # cumulative distribution (Gale & Grigal 1987)
    R = Y[1:] - Y[:-1]  # root area density distribution
# TESTI, SET FIRST LAYER WITH NO ROOTS
    R[0] = 0.0

    # addjust distribution to match soil profile depth
    R = R / sum(R) / dz[:len(R)]

    return R