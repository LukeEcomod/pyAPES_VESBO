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
        self.f = np.maximum(self.fmin, np.minimum(0.065*S, 1.0))

        if out:
            return self.f

class LAI_cycle(object):
    r"""Dercribes seasonal cycle of leaf-area index (LAI)
    """
    def __init__(self, p):
        r""" Initializes LAI cycle model.

        Args:
            'laip' (dict): parameters forleaf-area seasonal dynamics
                'lai_min': minimum LAI, fraction of annual maximum [-]
                'lai_ini': initial LAI fraction, if None lai_ini = Lai_min * LAImax
                'DDsum0': degreedays at initial time [days]
                'Tbase': base temperature [degC]
                'ddo': degreedays at bud burst [days]
                'ddur': duration of recovery period [days]
                'sso': start doy of decrease, based on daylength [days]
                'sdur': duration of decreasing period [days]
        Returns:
            self (object)
        """
        self.LAImin = p['lai_min']  # minimum LAI, fraction of annual maximum
        self.ddo = p['ddo']
        self.ddur = p['ddur']
        self.sso = p['sso']
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

        # update LAI
        self.f = f
        if out:
            return f