# -*- coding: utf-8 -*-
"""
Functions and classes for micrometeorological  transport processes 
and physical processes in atmospheric surface layer

@author: L1656
"""
import numpy as np

# Constants used in the model calculations.
#: [-], von Karman constant
VON_KARMAN = 0.41

class Aerodynamics():
    """
    Computes wind speed at ground and canopy + boundary layer conductances
    Computes wind speed at ground height assuming logarithmic profile above and
    exponential within canopy.
    Refs:
       Cammalleri et al. 2010 Hydrol. Earth Syst. Sci
       Massman 1987, BLM 40, 179 - 197.
       Magnani et al. 1998 Plant Cell Env.
    """
    def __init__(self, p):
        """
        Args:
            p - parameter dict
                w - leaf length scale [m]
                zm - wind speed measurement height above canopy [m]
                zg - height above ground where Ug is computed [m]
                zos - forest floor roughness length, ~ 0.1*roughness element height [m]
        Returns:
            Aerodynomics -object
        """

        # parameters
        self.w = p['w']  # leaf length scale [m]
        self.zmeas = p['zmeas']  # wind speed measurement height above canopy [m]
        self.zg = p['zg']  # height above ground where Ug is computed [m]
        self.zos = p['zos']  # forest floor roughness length [m]

    def _run(self, LAI, hc, Uo):
        """
        Args:
            LAI - one-sided leaf-area /plant area index (m2m-2)
            hc - canopy height (m)
            Uo - mean wind speed at height zm (ms-1)
        Returns:
            ra - canopy aerodynamic resistance (s m-1)
            rb - canopy boundary layer resistance (s m-1)
            ras - forest floor aerod. resistance (s m-1)
            ustar - friction velocity (m s-1)
            Uh - wind speed at hc (m s-1)
            Ug - wind speed at zg (m s-1)
        """
        zm = hc + self.zmeas  # m
        beta = 285.0  # s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3
        alpha = LAI / 2.0  # wind attenuation coeff (Yi, 2008 eq. 23)
        d = 0.66 * hc  # m
        zom = 0.123 * hc  # m
        zov = 0.1 * zom
        zosv = 0.1 * self.zos

        # solve ustar and U(hc) from log-profile above canopy
        ustar = Uo * VON_KARMAN / np.log((zm - d) / zom) 
        Uh = ustar / VON_KARMAN * np.log((hc - d) / zom)
        
        # U(zg) from exponential wind profile
        zn = np.minimum(self.zg / hc, 1.0)  # zground can't be above canopy top
        Ug = Uh * np.exp(alpha * (zn - 1.0))

        # canopy aerodynamic & boundary-layer resistances (sm-1). Magnani et al. 1998 PCE eq. B1 & B5
        #ra = 1. / (kv*ustar) * np.log((zm - d) / zom)
        ra = 1./(VON_KARMAN**2.0 * Uo) * np.log((zm - d) / zom) * np.log((zm - d) / zov)
        rb = 1. / LAI * beta * ((self.w / Uh)*(alpha / (1.0 - np.exp(-alpha / 2.0))))**0.5
        ra = ra + rb

        # soil aerodynamic resistance (sm-1)
        ras = 1. / (VON_KARMAN**2.0*Ug) * np.log(self.zg / self.zos)*np.log(self.zg / zosv)

        return ra, rb, ras, ustar, Uh, Ug
