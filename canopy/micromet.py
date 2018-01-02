# -*- coding: utf-8 -*-

"""
Module micromet 
Contains general methods related to Micrometeorological transport processes and physical processes in atmospheric surface layer.
Call all methods as e.g. "cMet.flog_Wind_profile(args)", where "cMet" is instance of cMicroMet -class. 
See documentation of interface of each method for inputs, outputs and datatypes.
Uses numpy for array manipulation & calculation.
    METHODS:

    AUTHOR:
        Samuli Launiainen, METLA 1/2011 - 4/2014 
    VERSION:
        15.04.2014: Codes converted from Matlab to Python. Not fully tested!! \n     
"""

import numpy as np
import matplotlib.pyplot as plt
eps = np.finfo(float).eps  # machine epsilon

""" define constants """
        
NT = 273.15  # 0 degC in Kelvin
NP = 101300.0  # Pa, sea level normal pressure
R = 8.314462175  # J mol-1 K-1, universal gas constant. One gets specific gas constant by R/M where M is molar mass
CP_AIR_MOLAR = 29.3  # J mol-1 K-1 molar heat capacity of air at constant pressure
CP_AIR_MASS = 1004.67  # J kg-1 K-1 heat capasity of the air at constant pressure
MAIR_DRY = 28.964e-3  # kg mol-1, molar mass of dry air
MH2O = 18.02e-3  # kg mol-1, molar mass of H2O
MCO2 = 44.01e-3  # kg mol-1, molar mass of CO2

SIGMA = 5.6697e-8  # Stefan-Boltzman constant W m-2 K-4
VON_KARMAN = 0.41  # von Kárman constant (-)
GRAV_CONST = 9.81  # ms-2 acceleration due gravity
DEG_TO_RAD = np.pi/180.0  # conversion deg -->rad
RAD_TO_DEG = 180.0/np.pi  # conversion rad -->deg


def test_canopyflow():
    """ testing closure_1_model_U """
    z = np.linspace(0, 20, 200)
    hc = 15.0
    LAI = 3.0
    b, c = 0.43, 1.83
    lad = generate_lad_weibul(z, LAI, hc, b, c)

    plt.figure()
    plt.plot(lad, z)

    Utop = 2.5
    Ubot = 0.0
    Cd = 0.15
    dPdx = 0.0
    tau, Un, Km, l_mix, d, zo = closure_1_model_U(z, Cd, lad, hc, Utop, Ubot, dPdx=dPdx, lbc_flux=None)


def generate_lad_weibul(z, LAI, h, b, c):
    """
    Generates leaf-area density profile from Weibull-distribution
    INPUT:
        z: height array (m), monotonic and constant steps
        LAI: leaf-area index (m2m-2)
        h: canopy height (m), scalar
        b: Weibull shape parameter 1, scalar
        c: Weibull shape parameter 2, scalar
    OUTPUT:
        LAD: leaf-area density (m2m-3), array \n
    SOURCE:
        Teske, M.E., and H.W. Thistle, 2004, A library of forest canopy structure for use in interception modeling.
        Forest Ecology and Management, 198, 341-350. 
        Note: their formula is missing brackets for the scale param. 
    AUTHOR:
        Gabriel Katul, 2009. Coverted to Python 16.4.2014 / Samuli Launiainen
    """
    z = np.array(z)
    dz = abs(z[1]-z[0])
    N = np.size(z)
    LAD = np.zeros(N)

    # dummy variables
    a = np.zeros(N)
    x = z[z <= h] / h  # normalized heigth
    n = np.size(x)

    # weibul-distribution
    a[0:n] = -(c/b)*(((1.0-x)/b)**(c-1.0))*(np.exp(-((1.0-x) / b)**c)) / (1.0-np.exp(-(1.0 / b)**c))
    # a[0] = 0.0  # no leaves at ground
    a = a / sum(a*dz)

    LAD = LAI*a
    # plt.figure(1)
    # plt.plot(LAD,z,'r-')      
    return LAD
    
def diffusivities_in_air(T):
    """
    Computes scalar diffusivities in still air as function of temperature.
    INPUT:
        T - air temperature (degC), scalar or array
    OUTPUT:
        Dt - heat diffusivity (m2s-1), scalars or arrays
        Dv - water vapor diffusivity(m2s-1)
        Dc - CO2 diffusivity (m2s-1)
        Do3 - O3 diffusivity (m2s-1)
    SOURCE:
        Based on tabulated diffusivities in Campbell and Norman, 1998.
        Introduction to Environmental Biophysics, Springer.
        Linear least squares fitting to data; valid at least in typical range ambient T (SL 25.4.2014).
    """    
    # diffusivities in m2 s-1
    Dt = 1e-6*(18.8 + 0.128*T)
    Dv = 1e-6*(21.2 + 0.1436*T)
    Dc = 1e-6*(13.8 + 0.096*T)
    Do3 = 1e-6*(17.6 + 0.12*T)
    
    return Dt, Dv, Dc, Do3
    
def pressure_from_altitude(ASL):
    """
    Approximates station pressure from site altitude
    INPUT:
        ASL - elevation above sea level (m)
    OUTPUT:
        Pamb - station pressure (Pa), assuming sea level at NP=101300 Pa
    SOURCE:
        Campbell & Norman, 1998. Introduction to Environmental biophysics, Springer.
    """
    ASL = np.array(ASL)
    Pamb = NP*np.exp(-(ASL/8200.0))
    return Pamb

"""Functions related to saturation vapor pressure and phase changes of H2O """

   
def latent_heat_vaporization(T, units="molar"):
    """
    Temperature dependency of latent heat of vaporization
    INPUT:
        T - temperature (degC)
        units - output units, "mass" = J kg-1 , "molar"= J mol-1
    OUTPUT:
        Lv - latent heat of vaporization in desired units
    """

    if np.any(T > 200):
        T = T - NT  #T must be in degC

    Lv = 1.0e6*(2.501 - 2.361e-3*T)*MH2O  # J mol-1
    
    if units=="mass":
        Lv = Lv / MH2O  # J kg-1
    return Lv
        
            
def saturation_vapor_pressure(T):
    """
    Computes saturation vapor pressure with respect to free and flat water surface for given temperature T
    INPUT:
        T - temperature (degC), scalar or array
    OUTPUT: 
        esat - saturation vapor pressure (Pa)
        delta - slope of saturation vapor pressure curve (Pa degC-1)
    SOURCE:
        Campbell & Norman, 1998. Introduction to Environmental Biophysics.
    """
    # constants
    a = 611.0  # Pa
    b = 17.502  # (-)
    c = 240.97  # degC

    esat = a*np.exp(b*T / (T+c))  # Pa
    delta = b*c*esat / ((c + T)**2)  # Pa degC-1
    return esat, delta


def psycrometric_constant(T, Pamb=101300):
    """
    Computes Psycrometric constant at temperature T
    INPUT:
        T - temperature (degC)
        Pamb - ambient pressure (Pa)
    OUTPUT:
        g - psychrometric constant
    USES:
        latent_heat_vaporization
    """
    Lv_mass = latent_heat_vaporization(T, units="mass")  # J kg-1
    g = Pamb*CP_AIR_MASS / (0.622*Lv_mass)  # Pa K-1
    return g


def vpd_from_rh(T, RH):
    """
    Computes vapor pressure deficit from temperature and relative humidity
    INPUT:
        T - temperature (degC), array or scalar
        RH - relative humidity (%), array or scalar
    OUTPUT:
        VPD - vapor pressure deficit (Pa), array
    USES:
        saturation_vapor_pressure
    """
    RH = np.array(RH)
    T = np.array(T)
    SVP, _ = saturation_vapor_pressure(T)
    VPD = (1.0 - RH / 100.0) * SVP  # Pa
    return VPD


def air_density(T, P=101300.0, h2o=0.0, units="mass"):
    """
    Computes air density at temperature T, pressure P and vapor pressure H2O
    INPUT:
        T - air temperature (degC), scalar or array
        P - ambient pressure (Pa), scalar or array, optional
        H2O - water vapor partial pressure (Pa), scalar or array, optional (default = dry air)
        units - units to return density: "mass" (default), "molar")
    OUTPUT:
        rhoa - density of dry (default) or moist air (kg m-3 or mol m-3), scalar or array
    Samuli Launiainen 28.4.2014
    """

    Tk = T + NT  # K

    # partial pressures of ideal gas are additive
    Pdry = P - h2o  # pressure of dry air

    if units == "mass":
        rhoa = (Pdry*MAIR_DRY + h2o*MH2O) / (R*Tk)  # kg m-3

    elif units == "molar":
        rho_d = Pdry / (R*Tk)  # dry air, mol m-3
        rho_v = h2o / (R*Tk)  # water vapor, mol m-3
        rhoa = rho_d + rho_v

    else:
        print "-----micromet.air_density: Error - check output units requested ---------"
        rhoa = np.nan
    return rhoa

# ---- evapotranspiration
def eq_evap(AE, T, P=101300.0, units='W'):
    """
    Calculates the equilibrium evaporation according to McNaughton & Spriggs, 1986. \n
    INPUT: 
        AE - Available energy (Wm-2)  
        T - air temperature (C)
        P - pressure (Pa)
        units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    OUTPUT: 
        equilibrium evaporation rate (Wm-2)
        constants
        NT=273.15; %0 degC in K
    """
    NT = 273.15
    Mw = 18e-3  # kg mol-1
    L = 1e3*(3147.5 - 2.37*(T + NT))  # latent heat of vaporization of water [J/kg]
    _, s = saturation_vapor_pressure(T)  # des / dT, Pa
  
    g = P*CP_AIR_MASS / (0.622*L)  # psychrom. const (Pa)

    x = np.divide((AE*s), (s+g))  # Wm-2 = Js-1m-2
    if units == 'mm':
        x = x / L  # kg m-2 s-1 = mm s-1
    elif units == 'mol':
        x = x / L / Mw  # mol m-2 s-1

    return x    


def penman_monteith(AE, D, T, Gs, Ga, P=101300.0, units='W'):
    """
    Computes latent heat flux LE (Wm-2) i.e evapotranspiration rate ET (mm/s)
    from Penman-Monteith equation
    INPUT:
       AE - available energy [Wm-2]
       VPD - vapor pressure deficit [Pa]
       T - ambient air temperature [degC]
       Gs - surface conductance [ms-1]
       Ga - aerodynamic conductance [ms-1]
       P - ambient pressure [Pa]
       units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    OUTPUT:
       x - evaporation rate in 'units'
    """
    # --- constants
    cp = 1004.67  # J kg-1 K-1
    rho = 1.25  # kg m-3
    Mw = 18e-3  # kg mol-1
    _, s = saturation_vapor_pressure(T)  # slope of sat. vapor pressure curve
    g = psycrometric_constant(T)

    L = 1e3 * (3147.5 - 2.37 * (T + 273.15))

    x = (s * AE + rho * cp * Ga * D) / (s + g * (1.0 + Ga / Gs))  # Wm-2

    if units is 'mm':
        x = x / L  # kgm-2s-1 = mms-1
    if units is 'mol':
        x = x / L / Mw  # mol m-2 s-1

    x = np.maximum(x, 0.0)
    return x

#------------- Flow within and above roughness sub-layer ------------------------
    
def wind_profile_exponential(z, hc, LAI, Uh, beta=1.0):
    """
    Computes exponential wind profile U(z) with respect to height z inside canopy
    INPUT:
        z - height array (m) 
        hc - canopy height (m), scalar
        LAI - one-sided leaf-area /plant area index (m2m-2)
        Uh - mean wind speed at canopy top, (ms-1)
        beta - attenuation coefficient = Ustar/Uh -ratio (-) at canopy top, optional
    OUTPUT:
        Uc - wind profile with respective to height, array
    SOURCE:
        Inoue, 1963: On the turbulent structure of air flow within crop 
        canopies. J. Meteorol Soc Jpn 41, 317-326
    AUTHOR:
        Samuli Launiainen 28.4.2014
    """
    z = np.array(z)       
    Uc = Uh*np.exp(beta*(z / hc - 1))  # m s-1
    return Uc

def wind_profile_loglinear(z, d, zom, zeta=1.0, ustar=None):
    """
    Log-linear wind profile in atm. surface layer.
    Args:
        z - height (scalar or array) (m)
        d - zero-place discplacement height (m)
        zom - roughness height for momentum (m)
        zeta - dimensionless stability parameter (z-d)/L (-)
        ustar - if None, normalized U/ustar is given as output
    Returns:
        Un - (m/s) is 'ustar' given as argument, else U/ustar (-)
    """
    
    zeta = np.array(zeta, ndmin=1)
    psi_m = np.zeros(len(zeta))

    # stable & neutral
    psi_m[zeta >=0] = 4.7 * zeta[zeta >= 0]
    
    # unstable
    x = (1.0 - 15*zeta[zeta < 0])**-0.25
    psi_m[zeta < 0] = -2.0*np.log(0.5*(1.0 + x)) - np.log(0.5*(1.0 + x**2)) + 2*np.arctan(x) - 0.5*np.pi
    
    Un = 1. / VON_KARMAN * (np.log((z-d) / zom) + psi_m)  # m/s
    
    if ustar:
        Un = Un * ustar
    return Un
 

def wind_profile_canopy(z, LAI, Uo, zm, hc, d, zom, beta=2.0):
    """
    Mean wind profile U(z) with respect to height z above and inside canopy in near-neutral conditions.
    Above canopy top logarithmic and within canopy exponential profile is assumed.
    INPUT:
        z - height array (m) 
        LAI - one-sided leaf-area /plant area index (m2m-2)
        Uo - mean wind speed at height zm (ms-1)
        zm - wind speed measurement height (m)
        hc - canopy height (m)
        d - zero-place discplacement height (m)
        zom - roughness length for momentum (m)
        beta - attenuation coefficient
    OUTPUT:
        Uc - wind profile with respective to height, array
    SOURCE:
        Inoue, 1963: On the turbulent structure of air flow within crop 
        canopies. J. Meteorol Soc Jpn 41, 317-326
    AUTHOR:
        Samuli Launiainen 7.9.2017. Follows Campbell & Norman, 1998 chapter 5.4 & 5.5.
    """    

    U = np.ones(len(z))*np.NaN

    # solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * VON_KARMAN / np.log((zm - d) / zom)  # m/s

    # above canopy top wind profile is logarithmic
    U[z >= hc] =  ustar / VON_KARMAN * np.log((z[z >= hc] - d) / zom) 

    # at canopy top, match log and exponential profiles
    Uh = ustar / VON_KARMAN * np.log((hc - d) / zom)  # m/s
    U[z <= hc] = Uh*np.exp(beta*(z[z <=hc] / hc - 1.0))  # m s-1
    
    return U

def most(zeta):
    """ 
    Monin-Obukhov stability functions for momentum in Atmospheric Surface Layer 
    (Businger & Dyer, from Stull chapter 9.7)
    INPUT:
        zeta - dimensionless stability parameter (z-d)/L, where z (m) is
        height above ground, d (m) zero-plane displacement height, and L (m) is Obukhov length, scalar
    OUTPUT:
        phi - dimensionless stability function, scalar        
    """
    zeta = float(zeta)
    if zeta < 0:
        phi = (1.0 - 15.0*zeta)**-0.25  # unstable
    if zeta >= 0:
        phi = 1.0 + 4.7*zeta  # neutral/stable

    return phi


def bulk_aerodynamic_conductance(U, zm, d, zom, zos=None):
    """
    Canopy bulk aerodynamic conductance in near-neutral conditions (ms-1)
    Args:
        U - mean wind speed (ms-1)
        zm - wind speed measurement height (m)
        d - zero place displacement height (m) 
        zom - roughness length for momentum (m)
        zov - roughness length of scalar s (m)
    Returns:
        ga - aerodynamic conductance (ms-1)
    Reference:
        Leuning et al. 2008. Water Resources Res. 44, W10419, doi:10.1029/2007WR006562.
    """
    if not zos:
        zos = zom

    ga = VON_KARMAN**2.0 * U / (np.log((zm - d) / zom) * np.log((zm - d) / zos))
    return ga


def bulk_aerodynamic_conductance_from_ust(Ust, U, Stanton):
    """
    Canopy bulk aerodynamic conductance (ms-1) from frict. velocity
    IN:
       Ust - friction velocity (ms-1)
       U - mean wind speed at flux measurement heigth (ms-1)
       Stanton - Stanton number (kB-1) for quasi-laminar boundary layer
           resistance. Typically kB=1...12, use 2 for vegetation ecosystems
           (Verma, 1989, Garratt and Hicks, 1973)
    OUT:
       Ga - aerodynamic conductance [ms-1]
    """
    ra = U / (Ust ** 2 + eps) + Stanton / (VON_KARMAN * (Ust + eps))  # sm-1
    Ga = 1.0 / ra  # ms-1
    return Ga
    
""" ----- 1st order closure models for canopy flow & scalar transport -----"""
def closure_1_model_U(z, Cd, lad, hc, Utop, Ubot, dPdx=0.0, lbc_flux=None):
    """
    Computes normalized mean velocity profile, shear stress and eddy diffusivity
    within and above horizontally homogenous plant canopies using 1st order closure.
    Accounts for horizontal pressure gradient force dPdx, assumes neutral diabatic stability.
    Soleves displacement height as centroid of drag force.
    IN:
       z - height (m), constant increments
       Cd - drag coefficient (0.1-0.3) (-)
       lad - plant area density, 1-sided (m2m-3)
       hc - canopy height (m)
       Utop - U /u* upper boundary
       Uhi - U /u* at ground (0.0 for no-slip)
       dPdx - u* -normalized horizontal pressure gradient
    OUT:
       tau - u* -normalized momentum flux 
       U - u* normalized mean wind speed (-) 
       Km - eddy diffusivity for momentum (m2s-1) 
       l_mix - mixing length (m)
       d - zero-plane displacement height (m)
       zo - roughness lenght for momentum (m)
    CODE:
        Gaby Katul, Samuli Launiainen 2008-2015. Converted to Python 17.5.2017
    """
    lad = 0.5*lad  # frontal plant-area density is half of one-sided
    dz = z[1] - z[2]
    N = len(z)
    U = np.linspace(Ubot, Utop, N)

    nn1 = max(2, np.floor(N/20))  # window for moving average smoothing

    # --- Start iterative solution
    err = 999.9
    eps1 = 0.1
    dPdx_m = 0.0

    while err > 0.01:
        Fd = Cd*lad*U**2  # drag force
        d = sum(z*Fd) / sum(Fd)  # displacement height
        l_mix = mixing_length(z, hc, d)  # m

        # --- dU/dz (m-1)
        y = central_diff(U, dz)
        # y = smooth(y, nn1)

        # --- eddy diffusivity & shear stress
        Km = l_mix**2*abs(y)
        # Km = smooth (Km, nn1)
        tau = -Km * y

        # ------ Set the elements of the Tri-diagonal Matrix
        a1 = -Km
        a2 = central_diff(-Km, dz)
        a3 = Cd*lad*U

        upd = (a1 / (dz*dz) + a2 / (2*dz))  # upper diagonal
        dia = (-a1*2 / (dz*dz) + a3)  # diagonal
        lod = (a1 / (dz*dz) - a2 / (2*dz))  # subdiagonal
        rhs = np.ones(N) * dPdx_m

        # upper BC
        upd[-1] = 0.
        dia[-1] = 1.
        lod[-1] = 0.
        rhs[-1] = Utop

        if not lbc_flux:  # --- lower BC, fixed Ubot
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = Ubot
        else:  # --- lower BC, flux-based
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = 0.  # zero-flux bc
            # how to formulate prescribed flux bc?
            # rhs[0] = lbc_flux

        # --- call tridiagonal solver
        Un = tridiag(lod, dia, upd, rhs)

        err = max(abs(Un - U))

        # --- Use successive relaxations in iterations
        U = eps1*Un + (1 - eps1)*U
        dPdx_m = eps1*dPdx + (1 - eps1)*dPdx_m

    # ---- return values
    tau = tau / tau[-1]  # normalized shear stress
    zo = (z[-1] - d)*np.exp(-0.4*U[-1])  # roughness length

    y = forward_diff(Un, dz)
    Kmr = l_mix**2 * abs(y)  # eddy diffusivity
    Km = smooth(Kmr, nn1)

    # --- for testing ----
    plt.figure(101)
    plt.subplot(221); plt.plot(Un, z, 'r-'); plt.title('U')
    plt.subplot(222); plt.plot(y, z, 'b-'); plt.title('dUdz')
    plt.subplot(223); plt.plot(l_mix, z, 'r-'); plt.title('l mix')
    plt.subplot(224); plt.plot(Km, z, 'r-', Kmr, z, 'b-'); plt.title('Km')    

    return tau, Un, Km, l_mix, d, zo
    

def mixing_length(z, h, d, l_min=None):
    """
    computes mixing length: linear above the canopy, constant within and
    decreases linearly close the ground (below z< alpha*h/kv)
    IN:
        z - computation grid, m
        h - canopy height, m
        d - displacement height, m
        l_min - intercept at ground, m
    OUT:
        lmix - mixing length, m
    """
    dz = z[1] - z[0]
    kv = 0.4  # von Karman constant

    if not l_min:
        l_min = dz / 2.0
    
    alpha = (h - d)*kv / h
    I_F = np.sign(z - h) + 1.0
    l_mix = alpha*h*(1 - I_F / 2) + (I_F / 2) * (kv*(z - d))
    
    sc = (alpha*h) / kv
    ix = np.where(z < sc)
    l_mix[ix] = kv*(z[ix] + dz / 2)
    l_mix = l_mix + l_min

    return l_mix

def forward_diff(y, dx):
    """
    computes gradient dy/dx using forward difference
    assumes dx is constatn
    """
    N = len(y)
    dy = np.ones(N) * np.NaN
    dy[0:-1] = np.diff(y)
    dy[-1] = dy[-2]
    return dy / dx


def central_diff(y, dx):
    """
    computes gradient dy/dx with central difference method
    assumes dx is constant
    """
    N = len(y)
    dydx = np.ones(N) * np.NaN
    # -- use central difference for estimating derivatives
    dydx[1:-1] = (y[2:] - y[0:-2]) / (2*dx)
    # -- use forward difference at lower boundary
    dydx[0] = (y[1] - y[0]) / dx
    # -- use backward difference at upper boundary
    dydx[-1] = (y[-1] - y[-2]) / dx

    return dydx

def smooth(a, WSZ):
    """
    smooth a by taking WSZ point moving average.
    NOTE: even WSZ is converted to next odd number.
    """
    WSZ = int(np.ceil(WSZ) // 2 * 2 + 1)
    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid') / WSZ
    r = np.arange(1, WSZ-1, 2)
    start = np.cumsum(a[:WSZ-1])[::2] / r
    stop = (np.cumsum(a[:-WSZ:-1])[::2] / r)[::-1]
    x = np.concatenate((start, out0, stop))
    return x

def closure_1_model_scalar(dz, Ks, Source, ubc, lbc, scalar,
                           T=20.0, P=101300.0, lbc_dirchlet=False):
    """
    Solves stationary scalar profiles in 1-D grid
    INPUT:
        dz - grid size (m), float
        Ks - eddy diffusivity (m2s-1), array
        Source - sink/source term, CO2 (umolm-3s-1), H2O (molm-3s-1), T (Wm-3),
            array
        ubc - upper boundary condition, value of CO2 (ppm), H2O (mol/mol)
            or T (degC) at highest gridpoint
        lbc - lower boundary condition, flux or value:
            flux: CO2 (umol m-2 s-1), H2O (mol m-2 s-1), T (Wm-2).
            value: CO2 (ppm), H2O (mol/mol), T (degC)
        scalar - 'CO2', 'H2O', 'T', string
        T - air temperature (degC), optional
        P - pressure (Pa)
        lbc_dirchlet - True for Dirchlet lower boundary
    OUT:
        Ca - mixing ratio profile (ppm)
    NOTE:
        assumes constant dz and Dirchlet upper boundary condition
    Gaby Katul & Samuli Launiainen, 2009 - 2017
    """
    epsi = 1e-5
    dz = float(dz)
    N = len(Ks)
    rho_a = P / (287.05*(T + 273.15))  # air density, kg/m3

    # -----Set elements of tridiagonal matrix
    a1 = Ks
    a2 = np.zeros(N)
    a2[1:] = np.diff(Ks) / dz
    a2[0] = a2[1]
    a3 = np.zeros(N)

    # super-diagonal, diagonal, subdiagonal elements and rhs vector
    upd = a1 / (dz*dz) + a2 / (2*dz)
    dia = -a1*2 / (dz*dz) + a3
    lod = a1 / (dz*dz) - a2 / (2*dz)
    rhs = np.zeros(N)

    # --- uppermost node, Dirchlet boundary
    upd[-1] = 0.
    dia[-1] = 1.
    lod[-1] = 0.

    # CO2
    if scalar.upper() == 'CO2':
        CF = 1e6*rho_a / (29e-3)  # umol m-3, molar conc. of air
        rhs = -Source / CF  # CO2 source/sink profile umol m-3 s-1 --> s-1.'
        rhs[-1] = 1e-6*ubc
        if not lbc_dirchlet:  # --- lower BC, flux-based
            print 'lbc'
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = (lbc / CF)*dz / (Ks[0] + epsi)
        else:  # --- lower BC, fixed concentration
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = 1e-6*lbc
        # --- call tridiagonal solver
        x = 1e6*tridiag(lod, dia, upd, rhs)  # ppm

    # H2O
    if scalar.upper() == 'H2O':
        CF = rho_a / (29e-3)  # mol m-3, molar conc. of air
        rhs = - Source / CF  # source/sink profile mol m-3 s-1 --> s-1.'
        rhs[-1] = ubc
        if not lbc_dirchlet:  # --- lower BC, flux-based
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = (lbc / CF)*dz / (Ks[0] + epsi)
        else:  # --- lower BC, fixed concentration
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = lbc
        # --- call tridiagonal solver
        x = tridiag(lod, dia, upd, rhs)  # mol/mol

    # Temperature
    if scalar.upper() == 'T':
        CF = 1004.0*rho_a  # cp*rho_a,  J m-3 K-1
        rhs = - Source / CF  # heat source/sink profile W m-3 --> K s-1.'
        rhs[-1] = ubc
        if not lbc_dirchlet:  # --- lower BC, flux-based
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = (lbc / CF)*dz / (Ks[0] + epsi)  # m s-1
            print rhs
        else:  # --- lower BC, fixed value
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = lbc
        # --- call tridiagonal solver
        x = tridiag(lod, dia, upd, rhs)  # degC
    return x


def tridiag(a, b, C, D):
    """
    tridiagonal matrix algorithm
    a=subdiag, b=diag, C=superdiag, D=rhs
    """
    n = len(a)
    V = np.zeros(n)
    G = np.zeros(n)
    U = np.zeros(n)
    x = np.zeros(n)

    V[0] = b[0].copy()
    G[0] = C[0] / V[0]
    U[0] = D[0] / V[0]

    for i in range(1, n):  # nr of nodes
        V[i] = b[i] - a[i] * G[i - 1]
        U[i] = (D[i] - a[i] * U[i - 1]) / V[i]
        G[i] = C[i] / V[i]

    x[-1] = U[-1]
    inn = n - 2
    for i in range(inn, -1, -1):
        x[i] = U[i] - G[i] * x[i + 1]
    return x
    
#def closure_1_model_CO2(dz, Kc, Csource, Ca0, lbc, T=20.0, P=101300.0, lbc_dirchlet=False):
#    """
#    INPUT:
#        dz - vertical grid size (m)
#        Kc - eddy diffusivity (m2s-1), array
#        Csource - CO2 sink/source term (umol m-3 s-1)
#        Ca0- CO2 mixing ratio at upper boundary (ppm)
#        lbc - lower boundary condition (umol m-2 ground s-1), OR ppm
#        T - air temperature (degC), optional
#        P - pressure (Pa)
#        lbc_dirchlet - True sets concentration boundary
#    OUT:
#        Ca - mixing ratio profile (ppm)
#    Gaby Katul & Samuli Launiainen, 2009 - 2017
#    """
#
#    # dz = z[1] - z[0]  # vertical increment (m)
#    N = len(Kc)
#
#    rho_a = P / (287.05*(T + 273.15))  # air density, kg/m3
#    CF = 1e6*rho_a / (29e-3)  # umol m-3, molar conc. of air
#
#    # -----Set up coefficients for ODE
#    a1 = Kc
#    a2 = np.zeros(N)
#    a2[1:] = np.diff(Kc) / dz
#    a2[0] = a2[1]
#    a3 = np.zeros(N)
#
#    # ---elements of the Tri-diagonal Matrix
#    upd = a1 / (dz*dz) + a2 / (2*dz)
#    dia = -a1*2 / (dz*dz) + a3
#    lod = a1 / (dz*dz) - a2 / (2*dz)
#    rhs = -Csource / CF  # CO2 source/sink profile umol m-3 s-1 --> s-1.'
#
#    if not lbc_dirchlet:
#        # --- lower BC, flux-based
#        upd[0] = -1.
#        dia[0] = 1.
#        lod[0] = 0.
#        rhs[0] = (lbc / CF)*dz / (Kc[0] + 0.00001)
#    else:
#        # --- lower BC, fixed concentration
#        upd[0] = 0.
#        dia[0] = 1.
#        lod[0] = 0.
#        rhs[0] = 1e-6*lbc
#
#    # --- upper BC
#    upd[-1] = 0.
#    dia[-1] = 1.
#    lod[-1] = 0.
#    rhs[-1] = 1e-6*Ca0  # umol mol-1
#
#    # --- Use tridiagonal solver to solve the tridiagonal matrix
#    Ca = 1e6*tridiag(lod, dia, upd, rhs)
#
#    return Ca
