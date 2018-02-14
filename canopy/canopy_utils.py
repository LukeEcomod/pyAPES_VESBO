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

def canopy_water_snow(W, SWE, LAI, cf, snowpara, dt, T, Prec, AE, VPD, Ra=25.0, U=2.0):
    """
    Calculates canopy water interception and SWE during timestep dt
    Args: 
        W : canopy water or snow storage [m]
        SWE (dict): snow state
            'SWE': snow water equivalent [m]
            'SWEi': water equivalent of ice in snow [m]
            'SWEl': water equivalent of liquid water in snow [m]
        LAI: leaf area index [m\ :sup:`2`\ m\ :sup:`-2`\]
        dt: timestep [s]
        T: air temperature [degC]
        Prec: precipitation rate during [m s\ :sup:`-1`\]
        AE: available energy (~net radiation) [W m\ :sup:`-2`\]
        VPD: vapor pressure deficit [kPa]
        Ra: canopy aerodynamic resistance [s m\ :sup:`-1`\]
    Returns:
        W : canopy water or snow storage [m]
        SWE (dict): snow state
            'SWE': snow water equivalent [m]
            'SWEi': water equivalent of ice in snow [m]
            'SWEl': water equivalent of liquid water in snow [m]
        Infil: potential infiltration to soil profile [m]
        Evap: evaporation / sublimation from canopy store [m]
        MBE: mass balance error [m]
    Samuli Launiainen & Ari Laurén 2014 - 2017
    Last edit 12 / 2017
    """

    # quality of precipitation [degC]
    Tmin = snowpara['Tmin']
    Tmax = snowpara['Tmax']
    Tmelt = snowpara['Tmelt']

    # interception storage capacities [m]
    Wmax = snowpara['wmax'] * LAI
    Wmaxsnow = snowpara['wmaxsnow'] * LAI

    # snow state [m]
    swei = SWE['SWEi']
    swel = SWE['SWEl']

    # melting and freezing coefficients [m/s]
    Kmelt = snowpara['kmelt'] - 1.64 / (24 * 3600 * 1000) * cf  # Kuusisto E, 'Lumi Suomessa'
    Kfreeze = snowpara['kfreeze']

    # fraction of Rn at ground
    kp = snowpara['kp']
    tau = np.exp(-kp * LAI)

    # inputs to arrays, needed for indexing later in the code
    gridshape = np.shape(LAI)  # rows, cols

    if np.shape(T) != gridshape:
        T = np.ones(gridshape) * T
        Prec = np.ones(gridshape) * Prec
        AE = np.ones(gridshape) * AE
        VPD = np.ones(gridshape) * VPD
        Ra = np.ones(gridshape) * Ra

    Prec = Prec * dt  # [m/s] -> [m]

    # latent heat of vaporization (Lv) and sublimation (Ls) [J/kg]
    Lv = 1e3 * (3147.5 - 2.37 * (T + DEG_TO_KELVIN))
    Ls = Lv + 3.3e5

    # compute 'potential' evaporation and sublimation rates for each grid cell  #### Onko GS ja GA näissä oikeinpäin?
    erate = np.zeros(gridshape)
    ixs = np.where((Prec == 0) & (T <= Tmin))
    ixr = np.where((Prec == 0) & (T > Tmin))
    Ga = 1. / Ra  # aerodynamic conductance [m/s]

    # resistance for snow sublimation adopted from:
    # Pomeroy et al. 1998 Hydrol proc; Essery et al. 2003 J. Climate;
    # Best et al. 2011 Geosci. Mod. Dev.
    # ri = (2/3*rhoi*r**2/Dw) / (Ce*Sh*W) == 7.68 / (Ce*Sh*W

    Ce = 0.01 * ((W + eps) / Wmaxsnow)**(-0.4)  # exposure coeff [-]
    Sh = (1.79 + 3.0 * U**0.5)  # Sherwood numbner [-]
    gi = Sh * W * 1000 * Ce / 7.68 + eps # [m/s]

    erate[ixs] = dt / Ls[ixs] / RHO_WATER * penman_monteith((1.0 - tau[ixs])*AE[ixs], 1e3*VPD[ixs], T[ixs], Ga[ixs], gi[ixs], units='W')
#        print('gi', gi, 'Ce', Ce, 'Sh', Sh)

    # evaporation of intercepted water, mm
    gs = 1e6
    erate[ixr] = dt / Lv[ixr] / RHO_WATER * penman_monteith((1.0 - tau[ixr])*AE[ixr], 1e3*VPD[ixr], T[ixr], Ga[ixr], gs, units='W')
    
    # print('erate', erate)

    # ---state of precipitation [as water (fW) or as snow(fS)]
    fW = np.zeros(gridshape)
    fS = np.zeros(gridshape)

    fW[T >= Tmax] = 1.0
    fS[T <= Tmin] = 1.0

    ix = np.where((T > Tmin) & (T < Tmax))
    fW[ix] = (T[ix] - Tmin) / (Tmax - Tmin)
    fS[ix] = 1.0 - fW[ix]
    del ix

    # --- Local fluxes (mm)
    Unload = np.zeros(gridshape)  # snow unloading
    Interc = np.zeros(gridshape)  # interception
    Melt = np.zeros(gridshape)   # melting
    Freeze = np.zeros(gridshape)  # freezing
    Evap = np.zeros(gridshape)

    """ --- initial conditions for calculating mass balance error --"""
    Wo = W  # canopy storage
    SWEo = SWE['SWE']  # Snow water equivalent m

    """ --------- Canopy water storage change -----"""
    # snow unloading from canopy, ensures also that seasonal LAI development does
    # not mess up computations
    ix = (T >= Tmax)
    Unload[ix] = np.maximum(W[ix] - Wmax[ix], 0.0)
    W = W - Unload
    del ix
    dW = W - Wo

    # Interception of rain or snow: asymptotic approach of saturation.
    # Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
    # Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
    ix = (T < Tmin)
    Interc[ix] = (Wmaxsnow[ix] - W[ix]) \
                * (1.0 - np.exp(-(cf[ix] / Wmaxsnow[ix]) * Prec[ix]))
    del ix
    
    # above Tmin, interception capacity equals that of liquid precip
    ix = (T >= Tmin)
    Interc[ix] = np.maximum(0.0, (Wmax[ix] - W[ix]))\
                * (1.0 - np.exp(-(cf[ix] / Wmax[ix]) * Prec[ix]))
    del ix
    W = W + Interc  # new canopy storage, mm

    Trfall = Prec + Unload - Interc  # Throughfall to field layer or snowpack

    # evaporate from canopy and update storage
    Evap = np.minimum(erate, W)  # mm
    W = W - Evap

    """ Snowpack (in case no snow, all Trfall routed to floor) """
    ix = np.where(T >= Tmelt)
    Melt[ix] = np.minimum(swei[ix], Kmelt[ix] * dt * (T[ix] - Tmelt))  # mm
    del ix
    ix = np.where(T < Tmelt)
    Freeze[ix] = np.minimum(swel[ix], Kfreeze * dt * (Tmelt - T[ix]))  # mm
    del ix

    # amount of water as ice and liquid in snowpack
    Sice = np.maximum(0.0, swei + fS * Trfall + Freeze - Melt)
    Sliq = np.maximum(0.0, swel + fW * Trfall - Freeze + Melt)

    PotInf = np.maximum(0.0, Sliq - Sice * snowpara['retention'])  # mm
    Sliq = np.maximum(0.0, Sliq - PotInf)  # mm, liquid water in snow

    # update Snowpack state variables
    SWE['SWEl'] = Sliq
    SWE['SWEi'] = Sice
    SWE['SWE'] = swel + swei
    
    # mass-balance error mm
    MBE = (W + SWE['SWE']) - (Wo + SWEo) - (Prec - Evap - PotInf)

    return W, SWE, PotInf, Trfall, Evap, Interc, MBE, erate, Unload, fS + fW

def dry_canopy_et(LAI, LAIconif, LAIdecid, physpara, soilrp, D, Qp, AE, Ta, \
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

    # --- gsref as LAI -weighted average of conifers and decid.
    gsref = 1.0 /LAI * (LAIconif * physpara['gsref_conif']
            + LAIdecid * physpara['gsref_decid'])

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
    Gc[np.isnan(Gc)] = eps

    """ --- transpiration rate --- """
    Tr = penman_monteith((1.-tau)*AE, 1e3*D, Ta, Gc, 1./Ra, units='m')
    Tr[Tr < 0] = 0.0

    """--- forest floor evaporation rate--- """
    # soil conductance is function of relative water availability
    gcs = 1.0 / soilrp * beta**2.0

    # beta = Wliq / FC; Best et al., 2011 Geosci. Model. Dev. JULES

    Efloor = penman_monteith(tau * f * AE, 1e3*D, Ta, gcs, 1./Ras, units='m')
    Efloor[f==0] = 0.0  # snow on ground restricts Efloor
    
    return Tr, Efloor, Gc  # , fQ, fD, fRew

def daylength(LAT, LON, DOY):
    """
    Computes daylength from location and day of year.

    Args:
        LAT, LON - in deg, float or arrays of floats
        doy - day of year, float or arrays of floats

    Returns:
        dl - daylength (hours), float or arrays of floats
    """
    CF = np.pi / 180.0  # conversion deg -->rad

    LAT = LAT*CF
    LON = LON*CF

    # ---> compute declination angle
    xx = 278.97 + 0.9856*DOY + 1.9165*np.sin((356.6 + 0.9856*DOY)*CF)
    DECL = np.arcsin(0.39785*np.sin(xx*CF))
    del xx

    # --- compute day length, the period when sun is above horizon
    # i.e. neglects civil twilight conditions
    cosZEN = 0.0
    dl = 2.0*np.arccos(cosZEN - np.sin(LAT)*np.sin(DECL) / (np.cos(LAT)*np.cos(DECL))) / CF / 15.0  # hours

    return dl

def aerodynamics(LAI, hc, Uo, w=0.01, zm=2.0, zg=0.5, zos=0.01):
    """
    computes wind speed at ground and canopy + boundary layer conductances
    Computes wind speed at ground height assuming logarithmic profile above and
    exponential within canopy
    Args:
        LAI - one-sided leaf-area /plant area index (m2m-2)
        hc - canopy height (m)
        Uo - mean wind speed at height zm (ms-1)
        w - leaf length scale (m)
        zm - wind speed measurement height above canopy (m)
        zg - height above ground where Ug is computed (m)
        zos - forest floor roughness length, ~ 0.1*roughness element height (m)
    Returns:
        ra - canopy aerodynamic resistance (s m-1)
        rb - canopy boundary layer resistance (s m-1)
        ras - forest floor aerod. resistance (s m-1)
        ustar - friction velocity (m s-1)
        Uh - wind speed at hc (m s-1)
        Ug - wind speed at zg (m s-1)
    SOURCE:
       Cammalleri et al. 2010 Hydrol. Earth Syst. Sci
       Massman 1987, BLM 40, 179 - 197.
       Magnani et al. 1998 Plant Cell Env.
    """
    zm = hc + zm  # m
    kv = 0.4  # von Karman constant (-)
    beta = 285.0  # s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3
    alpha = LAI / 2.0  # wind attenuation coeff (Yi, 2008 eq. 23)
    d = 0.66*hc  # m
    zom = 0.123*hc  # m
    zov = 0.1*zom
    zosv = 0.1*zos

    # solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * kv / np.log((zm - d) / zom) 
    Uh = ustar / kv * np.log((hc - d) / zom)
    
    # U(zg) from exponential wind profile
    zn = np.minimum(zg / hc, 1.0)  # zground can't be above canopy top
    Ug = Uh * np.exp(alpha*(zn - 1.0))

    # canopy aerodynamic & boundary-layer resistances (sm-1). Magnani et al. 1998 PCE eq. B1 & B5
    #ra = 1. / (kv*ustar) * np.log((zm - d) / zom)
    ra = 1./(kv**2.0 * Uo) * np.log((zm - d) / zom) * np.log((zm - d) / zov)    
    rb = 1. / LAI * beta * ((w / Uh)*(alpha / (1.0 - np.exp(-alpha / 2.0))))**0.5

    # soil aerodynamic resistance (sm-1)
    ras = 1. / (kv**2.0*Ug) * np.log(zg / zos)*np.log(zg / zosv)
    
    #print('ra', ra, 'rb', rb)
    ra = ra + rb
    return ra, rb, ras, ustar, Uh, Ug

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