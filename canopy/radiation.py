# -*- coding: utf-8 -*-
"""
Created on Tue Mar 06 11:06:37 2018

@author: L1656
"""
import numpy as np

#  conversion deg -->rad
DEG_TO_RAD = np.pi / 180.0
#  conversion rad -->deg
RAD_TO_DEG = 180.0 / np.pi

class Radiation():
    """
    
    """
    def __init__(self, p):
        """
        Args:
            p - parameter dict
                'clump': clumping index [-]
                'kd': extinction coeff for diffuse radiation [-]
                ?? leaf-angle distribution [-], shoot Par-albedo [-], soil (moss) Par-albedo [-]
        Returns:
            Radiation -object
        """
        # parameters
        self.clump = p['clump']
        self.kd = p['kd']

    def layerwise_Rnet(self, LAI, Rnet):
        """
        Computes net radiation (Rnet) available at canopy and 
        ground based on above canopy Rnet
        Args:
            LAI (float): canopy leaf area index [m2/m2]
            Rnet (float): net radiation [W/m2]
        Returns:
            Rnet_c (float): net radiation available at canopy [W/m2]
            Rnet_gr (float): net radiation available at ground [W/m2]
        """
        # Rnet at canopy and forest floor
        Rnet_c = (1.0 - np.exp(-self.kd * self.clump * LAI)) * Rnet
        Rnet_gr = np.exp(-self.kd * LAI) * Rnet

        return Rnet_c, Rnet_gr

def solar_angles(lat, lon, jday, timezone=+2.0):
    """
    computes zenith, azimuth and declination angles for given location and time
    Args:
        lat, lon (deg)
        jday - decimal day of year (float or array)
        timezone - > 0 when east from Greenwich
    Returns:
        zen, azim, decl - rad
        sunrise, sunset, daylength (minutes)
    Algorithm based on NOAA solar calculator: https://www.esrl.noaa.gov/gmd/grad/solcalc/
    Equations: https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
    """
    lat0 = lat * DEG_TO_RAD
    jday = np.array(jday, ndmin=1)

    # fract. year (rad)
    if np.max(jday) > 366:
        y = 2*np.pi / 366.0 * (jday - 1.0)        
    else:
        y = 2*np.pi / 365.0 * (jday - 1.0)
    
    # declination angle (rad)
    decl = (6.918e-3 - 0.399912*np.cos(y) + 7.0257e-2*np.sin(y) - 6.758e-3*np.cos(2.*y)
        + 9.07e-4*np.sin(2.*y) - 2.697e-3*np.cos(3.*y) + 1.48e-3*np.sin(3.*y))
    
    # equation of time (min)
    et = 229.18*(7.5e-5 + 1.868e-3*np.cos(y) - 3.2077e-2*np.sin(y) 
        - 1.4615e-2*np.cos(2.*y) - 4.0849e-2*np.sin(2.*y))
    # print et / 60.
    # hour angle
    offset = et + 4.*lon - 60.*timezone
    fday = np.modf(jday)[0]  # decimal day
    ha = DEG_TO_RAD * ((1440.0*fday + offset) / 4. - 180.)  # rad
    
    # zenith angle (rad)
    aa = np.sin(lat0)*np.sin(decl) + np.cos(lat0)*np.cos(decl)*np.cos(ha)
    zen = np.arccos(aa)
    del aa

    # azimuth angle, clockwise from north in rad
    aa = -(np.sin(decl) - np.sin(lat0)*np.cos(zen)) / (np.cos(lat0)*np.sin(zen))        
    azim = np.arccos(aa)
    
    # sunrise, sunset, daylength
    zen0 = 90.833 * DEG_TO_RAD  # zenith angle at sunries/sunset after refraction correction

    aa = np.cos(zen0) / (np.cos(lat0)*np.cos(decl)) - np.tan(lat0)*np.tan(decl)
    ha0 = np.arccos(aa) * RAD_TO_DEG

    sunrise = 720.0 - 4.*(lon + ha0) - et  # minutes
    sunset = 720.0 - 4.*(lon - ha0) - et  # minutes 

    daylength = (sunset - sunrise)  # minutes

    sunrise = sunrise + timezone
    sunrise[sunrise < 0] = sunrise[sunrise <0] + 1440.0

    sunset = sunset + timezone
    sunset[sunset > 1440] = sunset[sunset > 1440] - 1440.0        

    return zen, azim, decl, sunrise, sunset, daylength