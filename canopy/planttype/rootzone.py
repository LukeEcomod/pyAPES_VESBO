# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 10:49:59 2018

@author: L1656
"""

import numpy as np
import matplotlib as plt

class RootUptake():
    """
    Do we need to solve hroot???
    """
    def __init__(self, p, dz_soil, LAImax):
        """
        Args:
            p - parameter dict
        Returns:
            Roots -instance
        """
        # parameters
        self.RAI = p['RAI_LAI_multiplier']*LAImax  # total fine root area index (m3/m2)
        self.rad = self.RAI * RootDistribution(p['beta'], dz_soil, p['root_depth'])
        self.fine_radius = p['fine_radius']  # fine root radius [m]
        self.radial_K = p['radial_K']  # maximum bulk root membrane conductance in radial direction [s-1]

        # state variables
        

def RootDistribution(beta, dz, root_depth):
    """
    Returns cumulative (Y) and root area density (R) distribution
    with depth. sum(Y)=1.
    Uses model of Gale and Grigal, 1987 Can. J. For.Res., 17, 829 - 834.
    Args:
        beta: dimensionless factor
        dz: soil layers from top to bottom [m]
        root_depth: depth of rooting zone [m]
    """
    z = np.concatenate([[0.0], np.cumsum(dz)])
    z = np.concatenate([z[z < root_depth], [root_depth]])
    d = abs(z * 100.0 )  # depth in cm

    Y = 1.0 - beta**d  # cumulative distribution (Gale & Grigal 1987)
    R = Y[1:] - Y[:-1]  # root area density distribution

    # addjust distribution to match soil profile depth
    R = R / sum(R)

    return R