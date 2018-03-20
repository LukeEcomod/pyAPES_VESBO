# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 13:44:31 2018

@author: L1656
"""

import numpy as np
eps = np.finfo(float).eps  # machine epsilon

# Constants used in the model calculations.
#: [J kg-1 K-1], Specific heat of air
CP = 1004.67
#: [kg m-3], Density of water
RHO_WATER = 1000.0
#: [kg m-3], Density of air
RHO_AIR = 1.25
#: [K], zero degrees celsius in Kelvin
DEG_TO_KELVIN = 273.15
#: [kg mol\ :sup:`-1`\ ], molar mass of H\ :sub:`2`\ O
MOLAR_MASS_H2O = 18.015e-3
#: degrees to radians
DEG_TO_RAD =  np.pi / 180.0

class Canopy_Transpiration():
    """
    Computes dry canopy transpiration 
    """
    def __init__(self, p):
        """
        Args:
            p - parameter dict
                'kp': attenuation coefficient for PAR [-]
                'q50': half-sat. of leaf light response [W m-2]
                'rw': transpiration moisture response parameter
                'rwmin' transpiration moisture response parameter
        Returns:
            Canopy_Transpiration instance
        """
        self.kp = p['kp']
        self.q50 = p['q50']
        self.rw = p['rw']
        self.rwmin = p['rwmin']

    def _run(self, dt, LAI, gsref, T, AE, Qp, H2O, P, Ra, CO2, fPheno, Rew):
        """
        Computes ET from 2-layer canopy in absense of intercepted preciptiation,
        i.e. in dry-canopy conditions

        Args:
           dt: timestep [s]
           LAI: canopy leaf area index [m2/m2]
           gsref: canopy ...
           T: air temperature [degC]
           H2O: mixing ratio [mol/mol]
           P: ambient pressure [Pa]
           AE: available energy at canopy heigth [W m-2]
           Qp: PAR in [Wm-2]
           Ra: canopy aerodynamic resistance [s m-1]
           CO2: atm. CO2 mixing ratio [ppm]
           fPheno:  phenology modifier [-]
           Rew: relative extractable water [-]
        Returns:
           Tr: canopy transpiration [m]
        Refs:
        Launiainen et al. (2016). Do the energy fluxes and surface conductance
        of boreal coniferous forests in Europe scale with leaf area?
        Global Change Biol.
        Modified from: Leuning et al. 2008. A Simple surface conductance model
        to estimate regional evaporation using MODIS leaf area index and the
        Penman-Montheith equation. Water. Resources. Res., 44, W10419
        Original idea Kelliher et al. (1995). Maximum conductances for
        evaporation from global vegetation types. Agric. For. Met 85, 135-147

        Samuli Launiainen, Luke
        Last edit: 12 / 2017
        """

    # vapor pressure deficit [Pa]
        esat, _, _ = e_sat(T)  # [Pa]
        VPD = max(0.0, esat - H2O * P)

        """--- canopy conductance Gc (integrated stomatal conductance)----- """
        # fQ: Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
        fQ = 1.0 / self.kp * np.log((self.kp * Qp + self.q50) / (self.kp \
              * Qp * np.exp(-self.kp * LAI) + self.q50 + eps))
        # VPD -response
        fD = 1.0 / (np.sqrt(1e-3 * VPD) + eps)
        # soil moisture response: Lagergren & Lindroth, xxxx"""
        fRew = np.minimum(1.0, np.maximum(Rew / self.rw, self.rwmin))
        # CO2 -response of canopy conductance, derived from APES-simulations
        # (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
        fCO2 = 1.0 - 0.387 * np.log(CO2 / 380.0)
        # canopy conductance
        Gc = 1000 * gsref * fQ * fD * fRew * fCO2 * fPheno
        if np.isnan(Gc):
            Gc = eps

        """ --- transpiration rate [m] --- """
        Tr = np.maximum(0.0, penman_monteith(AE, VPD, T, Gc, 1./Ra, units='m')) * dt

        return Tr

def penman_monteith(AE, D, T, Gs, Ga, P=101300.0, units='W', type='evaporation'):
    """
    Computes latent heat flux LE (W m-2) i.e evapotranspiration rate ET (m s-1)
    from Penman-Monteith equation
    INPUT:
       AE - available energy [W m-2]
       VPD - vapor pressure deficit [Pa]
       T - ambient air temperature [degC]
       Gs - surface conductance [ms-1]
       Ga - aerodynamic conductance [ms-1]
       P - ambient pressure [Pa]
       units - W (Wm-2), m (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
       type - 'evaporation' or 'sublimation' (only when 'm' or 'mol')
    OUTPUT:
       x - evaporation rate in 'units'
    """

    _, s, g = e_sat(T, P)  # slope of sat. vapor pressure, psycrom const
    L = latent_heat(T, type)  # latent heat of vaporization or sublimation [J/kg]

    x = (s * AE + RHO_AIR * CP * Ga * D) / (s + g * (1.0 + Ga /(Gs + eps)))  # [Wm-2]

    if units is 'm':
        x = x / L / RHO_WATER  # m
    if units is 'mol':
        x = x / L / MOLAR_MASS_H2O  # mol m-2 s-1

    x = np.maximum(x, 0.0)
    return x

def e_sat(T, P=101300):
    """
    Computes saturation vapor pressure (Pa), slope of vapor pressure curve
    [Pa K-1]  and psychrometric constant [Pa K-1]
    IN:
        T - air temperature (degC)
        P - ambient pressure (Pa)
    OUT:
        esa - saturation vapor pressure in Pa
        s - slope of saturation vapor pressure curve (Pa K-1)
        g - psychrometric constant (Pa K-1)
    """

    Lambda = 1e3 * (3147.5 - 2.37 * (T + DEG_TO_KELVIN))  # lat heat of vapor [J/kg]
    esa = 1e3 * (0.6112 * np.exp((17.67 * T) / (T + 273.16 - 29.66)))  # Pa

    s = 17.502 * 240.97 * esa / ((240.97 + T) ** 2)
    g = P * CP / (0.622 * Lambda)
    return esa, s, g

def latent_heat(T, type='evaporation'):
    """
    Computes latent heat of vaporization or sublimation [J/kg]
    Args:
        T: ambient air temperature [degC]
        type: 'evaporation' or 'sublimation' (only when 'm' or 'mol')
    Returns:
        L: latent heat of vaporization or sublimation depending on 'type'[J/kg]
    """
    # latent heat of vaporizati [J/kg]
    Lv = 1e3 * (3147.5 - 2.37 * (T + DEG_TO_KELVIN))
    if type == 'evaporation':
        return Lv
    if type == 'sublimation':
        # latent heat sublimation [J/kg]
        Ls = Lv + 3.3e5
        return Ls
