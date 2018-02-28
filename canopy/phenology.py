# -*- coding: utf-8 -*-
"""
Functions and classes for plant phenology, seasonal cycle of photosynthetic 
capacity and leaf-area development.
"""

import numpy as np

class Photo_cycle():
    """
    Ceasonal cycle of photosynthetic machinery.
    Adapted from:
        Mäkelä et al. 2004. Tree Physiology 24, 369–376; Formulation as in
        Kolari et al. 2007 Tellus.
    """
    def __init__(self, p):
        """
        Args:
            p - parameter dict
        Returns:
            Phenology -instance
        """
        self.tau = p['tau']  # time constant (days)
        self.Tbase = p['Tbase']  # base temperature (degC)
        self.Smax = p['smax']  # threshold for full acclimation (degC)
        self.fmin = p['fmin']  # minimum photocapacity (-)

        # state variables
        self.X = p['Xo']  # initial delayed temperature (degC)
        self.f = 1.0  # relative photocapacity

    def _run(self, T, out=False):
        """
        computes new stage of temperature acclimation and phenology modifier.
        Kolari et al. 2007 Tellus
        Args:
            T - daily mean air temperature
        Returns:
            self.f - phenology modifier [0...1], updates object state
        Note: Call once per day
        """
        self.X = self.X + 1.0 / self.tau * (T - self.X)  # degC

        S = np.maximum(self.X - self.Tbase, 0.0)
        # self.f = np.maximum(self.fmin, np.minimum(S / self.Smax, 1.0))  # peltoniemi et al 2015
        self.f = np.maximum(self.fmin, np.minimum(0.065*S, 1.0))  # Kolari et al. 2007 Tellus

        if out:
            return self.f

class LAI_cycle():
    """
    Seasonal cycle of leaf-area index (LAI)
    """

    def __init__(self, p):
        """
        Args:
            p - parameter dict:
                 'lai_min': minimum LAI, fraction of annual maximum [-]
                 'lai_ini': initial LAI fraction, if None lai_ini = Lai_min * LAImax
                 'DDsum0': degreedays at initial time [days]
                 'Tbase': base temperature [degC]
                 'ddo': degreedays at bud burst [days]
                 'ddur': duration of recovery period [days]
                 'sso': start doy of decrease, based on daylength
                 'sdur': duration of decreasing period [days]
        Returns:
            Seasonal_LAI -instance
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

    def _run(self, doy, T, out=False):
        """
        computes relative LAI
        Args:
            self - object
            doy - day of year
            T - daily mean temperature
        Returns:
            updates state variables: self._growth_stage,
            self._senec_stage, self.DDsum, self.f
            returns self.f, LAI relative to annual maximum (if out=True)
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