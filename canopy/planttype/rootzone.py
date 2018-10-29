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
import matplotlib as plt

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
        self.RAI = p['RAI_LAI_multiplier']*LAImax  # total fine root area index (m3/m2)
        self.rad = self.RAI * RootDistribution(p['beta'], dz_soil, p['root_depth'])
        self.fine_radius = p['fine_radius']  # fine root radius [m]
        self.radial_K = p['radial_K']  # maximum bulk root membrane conductance in radial direction [s-1]

        # state variables

def RootDistribution(beta, dz, root_depth):
    r""" Computes normalized root area density distribution with depth.

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
    d = abs(z * 100.0 )  # depth in cm

    Y = 1.0 - beta**d  # cumulative distribution (Gale & Grigal 1987)
    R = Y[1:] - Y[:-1]  # root area density distribution

    # addjust distribution to match soil profile depth
    R = R / sum(R)

    return R