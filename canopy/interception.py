# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 13:21:29 2018

@author: L1656
"""

import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from evapotranspiration import penman_monteith

class Interception():
    """
    Big-leaf interception model for rain and snowfall
    ****** Kersti! ********
    We need to change this to be multi-layer; see APES Matlab-codes. 
    But for simplicity let's first start with assumption that Tleaf == Tair and neglect
    leaf energy budget.
    """
    def __init__(self, p, LAI):
        """
        Args:
            dt: timestep [s]
            p - parameter dict
                'wmax':  maximum interception storage capacity for rain [m per unit of LAI]
                'wmaxsnow': maximum interception storage capacity for snow [m per unit of LAI]
                'Tmin': temperature below which all is snow [degC]
                'Tmax': temperature above which all is water [degC]
                'w_ini': initial canopy storage [m]
        Returns:
            Interception -object
        """
        # parameters:
        # maximum storage capacities [m per unit of LAI]
        self.wmax = p['wmax']  # for rainfall
        self.wmaxsnow = p['wmaxsnow']  # for snowfall
        # quality of precipitation [degC]
        self.Tmin = p['Tmin']
        self.Tmax = p['Tmax']
        self.Tmelt = p['Tmelt']

        # state variables
        self.W = np.minimum(p['w_ini'], p['wmax'] * LAI) # interception storage [m]

    def _run(self, dt, LAI, cf, T, Prec, AE, VPD, Ra=25.0, U=2.0):
        """
        Args:
            dt: timestep [s]
            LAI: leaf area index [m\ :sup:`2`\ m\ :sup:`-2`\]
            cf: canopy closure [-]
            T: air temperature [degC]
            Prec: precipitation rate [m s\ :sup:`-1`\]
            AE: available energy at canopy level  [W m\ :sup:`-2`\]
            VPD: vapor pressure deficit [kPa]
            Ra: canopy aerodynamic resistance [s m\ :sup:`-1`\]
            U: wind speed  at hc [m/s] ----------------------------------------------- ???
        Returns:
            updates sate self.W
            Trfall_rain: throughfall as rainfall [m]
            Trfall_snow: throughfall as snowfall [m]
            Interc: interception [m]
            Evap: evaporation from canopy store [m]
            MBE: mass balance error
        """
        # interception storage capacities [m]
        Wmax = self.wmax * LAI
        Wmaxsnow = self.wmaxsnow * LAI

        Prec = Prec * dt  # [m/s] -> [m]

        """ 'potential' evaporation and sublimation rates """
        if (Prec == 0) & (T <= self.Tmin):  # sublimation case
            Ce = 0.01 * ((self.W + eps) / Wmaxsnow)**(-0.4)  # exposure coeff [-]
            Sh = (1.79 + 3.0 * U**0.5)  # Sherwood numbner [-]
            gi = Sh * self.W * 1000 * Ce / 7.68 + eps # [m/s]
            erate = penman_monteith(AE, 1e3*VPD, T, gi, 1./Ra,  units='m', type='sublimation') * dt
        elif (Prec == 0) & (T > self.Tmin):  # evaporation case
            gs = 1e6
            erate = penman_monteith(AE, 1e3*VPD, T, gs, 1./Ra, units='m', type='evaporation') * dt
        else:  # negelect evaporation during precipitation events
            erate = 0.0
    
        """ --- state of precipitation --- """
        # fraction as water [-]
        if T >= self.Tmax:
            fW = 1.0
        elif T <= self.Tmin:
            fW = 0.0
        else:
            fW = (T - self.Tmin) / (self.Tmax - self.Tmin)

        """ --- canopy water storage change --- """
        W = self.W  # initiate

        # snow unloading from canopy, ensures also that seasonal
        # LAI development does not mess up computations
        if T >= self.Tmax and W > Wmax:
            Unload = W - Wmax
        else:
            Unload = 0.0
        # update canopy storage [m]
        W = W - Unload

        # interception of rain or snow: asymptotic approach of saturation.
        # Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
        # Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
        if T < self.Tmin:
            Interc = (Wmaxsnow - W) * (1.0 - np.exp(-(cf / Wmaxsnow) * Prec))
        else:  # above Tmin, interception capacity equals that of liquid precip
            Interc = np.maximum(0.0, (Wmax - W)) * (1.0 - np.exp(-(cf / Wmax) * Prec))
        # update canopy storage [m]
        W = W + Interc

        # throughfall to field layer or snowpack
        Trfall = Prec + Unload - Interc
        Trfall_rain = fW * Trfall
        Trfall_snow = (1 - fW) * Trfall

        # evaporate from canopy and update storage [m]
        Evap = np.minimum(erate, W)
        W = W - Evap

        # mass-balance error [m] ! self.W is old storage
        MBE = (W - self.W) - (Prec - Evap - (Trfall_rain + Trfall_snow))

        # update state variables
        self.W = W

        return Trfall_rain, Trfall_snow, Interc, Evap, MBE
