# -*- coding: utf-8 -*-
"""
UTILITY FUNCTION FOR CANOPY MODEL
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

def dry_canopy_et(LAI, gsref, physpara, soilrp, D, Qp, AE, Ta, \
                  Ra=25.0, Ras=250.0, CO2=380.0, f=1.0, Rew=1.0, beta=1.0, fPheno=1.0):
    """
    Computes ET from 2-layer canopy in absense of intercepted preciptiation, i.e. in dry-canopy conditions
    IN:
       self - object
       D - vpd in [kPa]
       Qp - PAR in [Wm-2]
       AE - available energy in [W m-2]
       Ta - air temperature [degC]
       Ra - aerodynamic resistance [s m-1]
       Ras - soil aerodynamic resistance [s m-1]
       CO2 - atm. CO2 mixing ratio [ppm]
       f - franction of local Rnet available for evaporation at ground [-]
       Rew - relative extractable water [-]
       beta - term for soil evaporation resistance (Wliq/FC) [-]
       fPheno - phenology modifier [-]
    Args:
       Tr - transpiration rate [m s-1]
       Efloor - forest floor evaporation rate [m s-1]
       Gc - canopy conductance (integrated stomatal conductance)  [m s-1]
    SOURCES:
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

    kp = physpara['kp']  # (-) attenuation coefficient for PAR
    q50 = physpara['q50']  # Wm-2, half-sat. of leaf light response
    rw = physpara['rw']  # rew parameter
    rwmin = physpara['rwmin']  # rew parameter

    tau = np.exp(-kp * LAI)  # fraction of Qp at ground relative to canopy top

    """--- canopy conductance Gc (integrated stomatal conductance)----- """

    # fQ: Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
    fQ = 1.0 / kp * np.log((kp * Qp + q50) / (kp*Qp*np.exp(-kp * LAI) + q50 + eps) )

    # the next formulation is from Leuning et al., 2008 WRR for daily Gc; they refer to 
    # Kelliher et al. 1995 AFM but the resulting equation is not exact integral of K95.        
    # fQ = 1./ kp * np.log((Qp + q50) / (Qp*np.exp(-kp*self.LAI) + q50))

    # VPD -response
    fD = 1.0 / (np.sqrt(D) + eps)

    # soil moisture response: Lagergren & Lindroth, xxxx"""
    fRew = np.minimum(1.0, np.maximum(Rew / rw, rwmin))
    # fRew = 1.0

    # CO2 -response of canopy conductance, derived from APES-simulations
    # (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
    fCO2 = 1.0 - 0.387 * np.log(CO2 / 380.0)

    Gc = 1000 * gsref * fQ * fD * fRew * fCO2 * fPheno
    if np.isnan(Gc):
        Gc = eps

    """ --- transpiration rate --- """
    Tr = np.maximum(0.0, penman_monteith((1.-tau)*AE, 1e3*D, Ta, Gc, 1./Ra, units='m'))

    """--- forest floor evaporation rate--- """
    # soil conductance is function of relative water availability
    gcs = 1.0 / soilrp * beta**2.0

    # beta = Wliq / FC; Best et al., 2011 Geosci. Model. Dev. JULES

    Efloor = penman_monteith(tau * f * AE, 1e3*D, Ta, gcs, 1./Ras, units='m')
    if f==0:
        Efloor = 0.0  # snow on ground restricts Efloor

    return Tr, Efloor, Gc  # , fQ, fD, fRew

def penman_monteith(AE, D, T, Gs, Ga, P=101300.0, units='W'):
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
       units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    OUTPUT:
       x - evaporation rate in 'units'
    """

    _, s, g = e_sat(T, P)  # slope of sat. vapor pressure, psycrom const
    L = 1e3 *(3147.5 - 2.37 * (T + 273.15))

    x = (s * AE + RHO_AIR * CP * Ga * D) / (s + g * (1.0 + Ga /(Gs + eps)))  # Wm-2

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