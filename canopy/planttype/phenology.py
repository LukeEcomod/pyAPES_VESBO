# -*- coding: utf-8 -*-
"""
.. module: phenology
    :synopsis: APES-model component
.. moduleauthor:: Samuli Launiainen & Kersti Haahti

Note:
    migrated to python3
    - nothing changed

Describes plant phenology, seasonal cycle of photosynthetic 
capacity and leaf-area development.

Created on Mon May 15 13:43:44 2017
"""

import numpy as np
from canopy.constants import DEG_TO_RAD

class Photo_cycle(object):
    r"""Seasonal cycle of photosynthetic machinery.

    References:
        Kolari et al. 2007 Tellus.
    """
    def __init__(self, p):
        r""" Initializes photo cycle model.

        Args:
            p (dict):
                'Xo': initial delayed temperature [degC]
                'fmin': minimum photocapacity [-]
                'Tbase': base temperature [degC]
                'tau': time constant [days]
                'smax': threshold for full acclimation [degC]
        Returns:
            self (object)
        """
        self.tau = p['tau']  # time constant (days)
        self.Tbase = p['Tbase']  # base temperature (degC)
        self.Smax = p['smax']  # threshold for full acclimation (degC)
        self.fmin = p['fmin']  # minimum photocapacity (-)

        # state variables
        self.X = p['Xo']  # initial delayed temperature (degC)
        self.f = 1.0  # relative photocapacity

    def run(self, T, out=False):
        r""" Computes new stage of temperature acclimation and phenology modifier.

        Args:
            T (float): mean daily air temperature [degC]
            out (bool): if true returns phenology modifier [0...1]

        Note: Call once per day
        """
        self.X = self.X + 1.0 / self.tau * (T - self.X)  # degC

        S = np.maximum(self.X - self.Tbase, 0.0)
        self.f = np.maximum(self.fmin, 
                            np.minimum(S / (self.Smax - self.Tbase), 1.0))

        if out:
            return self.f

class LAI_cycle(object):
    r"""Dercribes seasonal cycle of leaf-area index (LAI)
    """
    def __init__(self, p, loc):
        r""" Initializes LAI cycle model.

        Args:
            'laip' (dict): parameters forleaf-area seasonal dynamics
                'lai_min': minimum LAI, fraction of annual maximum [-]
                'lai_ini': initial LAI fraction, if None lai_ini = Lai_min * LAImax
                'DDsum0': degreedays at initial time [days]
                'Tbase': base temperature [degC]
                'ddo': degreedays at bud burst [days]
                'ddur': duration of recovery period [days]
                'sdl':  daylength for senescence start [h]
                'sdur': duration of decreasing period [days]
        Returns:
            self (object)
        """
        self.LAImin = p['lai_min']  # minimum LAI, fraction of annual maximum
        self.ddo = p['ddo']
        self.ddmat = p['ddmat']
        self.ddur = p['ddur']

        # senescence starts at first doy when daylength < sdl
        doy = np.arange(1, 366)
        dl = daylength(lat=loc['lat'], lon=loc['lon'], doy=doy)

        ix = np.max(np.where(dl > p['sdl']))
        self.sso = doy[ix]  # this is onset date for senescence

        self.sdur = p['sdur']
        if p['lai_ini']==None:
            self.f = p['lai_min']  # current relative LAI [...1]
        else:
            self.f = p['lai_ini']

        # degree-day model
        self.Tbase = p['Tbase']  # [degC]
        self.DDsum = p['DDsum0']  # [degC]

    def run(self, doy, T, out=False):
        r"""Computes relative LAI based on seasonal cycle.

        Args:
            T (float): mean daily air temperature [degC]
            out (bool): if true returns LAI relative to annual maximum

        Note: Call once per day
        """
        # update DDsum
        if doy == 1:  # reset in the beginning of the year
            self.DDsum = 0.
        else:
            self.DDsum += np.maximum(0.0, T - self.Tbase)
        """
        # growth phase
        if self.DDsum <= self.ddo:
            f = self.LAImin
            self._growth_stage = 0.
            self._senesc_stage = 0.
        elif self.DDsum > self.ddo:
            self._growth_stage += 1.0 / self.ddur
            f = np.minimum(1.0, self.LAImin + (1.0 - self.LAImin) * self._growth_stage)

        # senescence phase
        if doy > self.sso:
            self._growth_stage = 0.
            self._senesc_stage += 1.0 / self.sdur
            f = 1.0 - (1.0 - self.LAImin) * np.minimum(1.0, self._senesc_stage)
        """
        # growth phase
        if self.DDsum <= self.ddo:
            f = self.LAImin
        elif self.DDsum > self.ddo:
            f = np.minimum(1.0, self.LAImin + (1.0 - self.LAImin) *
                 (self.DDsum - self.ddo) / (self.ddmat - self.ddo))

        # senescence phase
        if doy > self.sso:
            f = 1.0 - (1.0 - self.LAImin) * np.minimum(1.0,
                    (doy - self.sso) / self.sdur)

        # update LAI
        self.f = f
        if out:
            return f

def daylength(lat, lon, doy):
    """
    Computes daylength from location and day of year.

    Args:
        lat, lon - in deg, float or arrays of floats
        doy - day of year, float or arrays of floats

    Returns:
        dl - daylength (hours), float or arrays of floats
    """

    lat = lat * DEG_TO_RAD
    lon = lon * DEG_TO_RAD

    # ---> compute declination angle
    xx = 278.97 + 0.9856 * doy + 1.9165 * np.sin((356.6 + 0.9856 * doy) * DEG_TO_RAD)
    decl = np.arcsin(0.39785 * np.sin(xx * DEG_TO_RAD))

    # --- compute day length, the period when sun is above horizon
    # i.e. neglects civil twilight conditions
    cosZEN = 0.0
    dl = 2.0 * np.arccos(cosZEN - np.sin(lat)*np.sin(decl) / 
                         (np.cos(lat)*np.cos(decl))) / DEG_TO_RAD / 15.0  # hours

    return dl
