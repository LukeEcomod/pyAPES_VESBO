# -*- coding: utf-8 -*-

"""
Package radiation: 
    Methods for computing short (SW) and long-wave (LW) radiation within and above plant canopies.
    See documentation of each method for inputs, outputs and datatypes
    Uses numpy for array manipulation & calculation
METHODS:
    test_radiation_functions - methods for testing SW and LW radiation functions
    generate_lad_weibul - weibul distribution function for creating leaf-area density distributions
    solar_dec_azimuth - solar zenith-, azimuth & inclination angle, daylength,sunrise,sunset
    cloudcover - retrospective determination of cloud fraction, fraction of diffuse radiation and atmospheric emissivity from surface observations \n       
    kbeam - beam attenuation coefficient
    kdiffuse - diffuse attenuation coefficient
    canopy_sw_ZhaoQualls - incident, reflected and absorbed SW in plant canopies, incl. multiple scattering
    canopy_sw_Spitters - incident and absorbed SW in plant canopies using simpler approach
    canopy-lw-ZhaoQualls - longwave radiation in plant canopies, incl. multiple scattering    
    canopy_lw - isothermal longwave radiation in plant canopies
DEPENDENCIES: numpy, matplotlib.pyplot
AUTHOR:
    Samuli Launiainen, METLA/Luke 1/2011 - 5/2017
VERSION:
    15.05.2017: Tested all functions   
"""

import numpy as np
import matplotlib.pyplot as plt
eps = np.finfo(float).eps  # machine epsilon

"""
Constants
"""
NT = 273.15  # 0 degC in K
SIGMA = 5.6697e-8  # Stefan-Boltzman constant W m-2 K-4
PAR_TO_UMOL = 4.56  # conversion from Wm-2 to micromol m-2 s-1
PAR_TO_WM2 = 1.0/4.56  # conversion from micromol m-2 s-1 to Wm-2
PAR_FRACTION = 0.45  # PAR fraction of global radiation (-), Jones H.(1992): Plants and Microclimate: ... Cambridge, 1992 
NIR_FRACTION = 0.55  # NIR fraction of global radiation (-)
DEG_TO_RAD = np.pi / 180.0  #  conversion deg -->rad
RAD_TO_DEG = 180.0 / np.pi  #  conversion rad -->deg
SOLAR_CONST = 1367.0  # Wm-2, solar constant at top of atm.

def test_radiation_functions(LAI, Clump, ZEN, method="canopy_sw_ZhaoQualls"):
    """
    Runs test script for SW and LW radiation methods.
    INPUT: 
        LAI: leaf area index (m2m-2)
        CLUMP: clumping index (-)
        ZEN: solar zenith angle (rad)            
        method: Name of function to be tested "functionname"
    OUTPUT:
        none, prints figures
    AUTHOR: 
        Samuli Launiainen, METLA 4/2014
    """

    """ define setup for testing models """

    #ZEN=35.0/180.0*np.pi        
    IbSky = 100.0
    IdSky = 100.0
    LeafAlbedo = 0.2
    SoilAlbedo = 0.2

    LAIz = np.ones(102)*float(LAI) / 100.0
    LAIz[0] = 0.
    LAIz[-1] = 0.
    N = len(LAIz)
    x = 1.0  # spherical leaf-angle distr.

    # for LW calculations
    T = np.linspace(15, 17, N) # Tair is 15degC at ground and 17 at upper boundary
    Tatm = 17
    Tsurf = 15
    LWdn0 = 0.85*SIGMA*(Tatm + NT)**4
    LWup0 = 0.98*SIGMA*(Tsurf + NT)**4
    print LWdn0, LWup0
    if method == "canopy_sw_ZhaoQualls":
        print "------TestRun of radiation.canopy_sw_ZhaoQualls with given LAI and CLUMP -----------"
        SWb, SWd, SWu, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_sl, alb = canopy_sw_ZhaoQualls(LAIz, Clump, x, ZEN, IbSky, IdSky, LeafAlbedo, SoilAlbedo, PlotFigs="True")                                                                         
        # print SWb,SWd,SWu,Q_sl,Q_sh,q_sl,q_sh,q_soil,f_sl,alb                                                                             

    if method == "canopy_sw_Spitters":
        print "------TestRun of radiation.canopy_sw_Spitters with given LAI and predefined lad profile-----------"
        SWb, SWd, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_sl, alb = canopy_sw_Spitters(LAIz, Clump, x, ZEN, IbSky, IdSky, LeafAlbedo, SoilAlbedo, PlotFigs="True")
        # print SWb, SWd, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_sl, alb  
    
    if method=="canopy_lw": 
        print "------TestRun of radiation.canopy_lw------------"
        LWnet, LWdn, WLup = canopy_lw(LAIz, Clump, T, LWdn0, LWup0, leaf_emi=1.0)

    if method == "canopy_lw_ZhaoQualls":
        print "------TestRun of radiation.canopy_lw_ZhaoQualls with given LAI and CLUMP -----------"
        LWnet, LWdn, WLup = canopy_lw_ZhaoQualls(LAIz, Clump, x, T, LWdn0, LWup0, leaf_emi=0.98, soil_emi=0.98, PlotFigs=True)   
            
                
def generate_lad_weibul(z, LAI, h, b, c):
    """
    Generates leaf-area density profile from Weibull-distribution
    INPUT:
        z: height array (m), monotonic and constant steps\n
        LAI: leaf-area index (m2m-2)\n
        h: canopy height (m), scalar \n
        b: Weibull shape parameter 1, scalar \n
        c: Weibull shape parameter 1, scalar \n
    OUTPUT:
        LAD: leaf-area density (m2m-3), array \n
    SOURCE:
        Teske, M.E., and H.W. Thistle, 2004, A library of forest canopy structure for use in interception modeling.
        Forest Ecology and Management, 198, 341-350. n.b. their formula is missing brackets for the scale param. \n
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
    a = a / sum(a*dz)

    LAD = LAI*a
    # plt.figure(1)
    # plt.plot(LAD,z,'r-')      
    return LAD

   
def solar_dec_azimuth(LAT, LON, DOY, HOUR, MINU, TimeZone="UTC+2"):
    """ 
    Computes solar angles, daylength and sunrise & sunset times
    INPUT:
        LAT: latitude (deg), + is northern hemisphere, scalar
        LON: longitude (deg), + is east from zero meridian, scalar
        DOY: julian day of year (scalar). 1st Jan = 1 
        HOUR: hour (local standard time), scalar
        MINU: minute, scalar
        TimeZone: optional parameter UTC time zone. Default = UTC+2 (Finland)
    OUTPUT:
        ZEN: solar zenith angle (rad)
        AZIM: solar asimuth angle (rad)
        DECL: solar declination angle (rad)
        DayLength: day length, i.e. sun above horizon (hours)
        t_sunrise: time of sunrise, i.e. sun at horizon (hours)
        t_sunset: time of sunset, i.e. sun at horizon
    SOURCE:
        Campbell & Norman, 1998: Introduction to environmental biophysics, Springer.   
    AUTHOR:
        Samuli Launiainen, METLA 04/2014
    """
    
    CF = DEG_TO_RAD  # conversion deg -->rad

    # dictionary of standard meridians of UTC timezones: note UTC key format in function call
    STMeridians = {'UTC+0': 0.0, 'UTC+1': 15.0, 'UTC+2': 30.0, 'UTC+3': 45.0, 
    'UTC+4': 60.0, 'UTC+5': 75.0, 'UTC+6': 90.0, 'UTC+7': 105.0, 'UTC+8': 120.0,
    'UTC+9': 135.0, 'UTC+10': 150.0, 'UTC+11': 165.0, 'UTC+12': 180.0,
    'UTC-1': -15.0, 'UTC-2': -30.0, 'UTC-3': -45.0,'UTC-4': -60.0, 'UTC-5': -75.0,
    'UTC-6': -90.0, 'UTC-7': -105.0, 'UTC-8': -120.0, 'UTC-9': -135.0,
    'UTC-10': -150.0, 'UTC-11': -165.0, 'UTC-12': -180.0}

    # input parameters
    LAT = float(LAT)*CF
    LON = float(LON)*CF
    DOY = float(DOY)
    tm = float(HOUR) + MINU / 60.0  # decimal hour of day
 
    LONo = STMeridians[TimeZone]  # standard meridian 
    # print LONo
    LC = 1.0 / 15.0*(LON / CF - LONo)  # longitude correction, hours

    # ---> compute declination angle
    xx = 278.97 + 0.9856*DOY + 1.9165*np.sin((356.6 + 0.9856*DOY)*CF)
    DECL = np.arcsin(0.39785*np.sin(xx*CF))
    del xx

    # ---> compute Zenith angle
    f = (279.575 + 0.9856*DOY)*CF
    ET = (-104.7*np.sin(f) + 596.2*np.sin(2.0*f) + 4.3*np.sin(3.0*f) - 12.7*np.sin(4*f) \
        - 429.3*np.cos(f) - 2*np.cos(2*f) + 19.3*np.cos(3.0*f)) / 3600.0

    to = 12.0 - LC - ET  # local solar noon
    aa = np.sin(LAT)*np.sin(DECL) + (np.cos(LAT))*np.cos(DECL)*np.cos(15.0*(tm-to)*CF)
    ZEN = np.arccos(aa)
    del aa,f

    # ---> compute azimuth angle	
    AZIM = np.arccos(-(np.sin(DECL)-(np.sin(LAT))*np.cos(ZEN)) / ((np.cos(LAT))*np.sin(ZEN)))

    # --- compute day length and sunrise/sunset times.
    # Consider period when sun is above horizon, i.e. neglect civil twilight conditions   
    cosZEN = 0
    DayLength = 2.0*np.arccos(cosZEN - np.sin(LAT)*np.sin(DECL) / (np.cos(LAT)*np.cos(DECL))) / CF / 15.0  # hours

    t_sunrise = to - DayLength / 2.0  # decimal hours
    t_sunset = to + DayLength / 2.0

    return ZEN, DECL, AZIM, DayLength, t_sunrise, t_sunset


def daylength(LAT, LON, DOY):
    """
    Computes daylength from location and day of year.
    """
    CF = DEG_TO_RAD  # conversion deg -->rad

    # input parameters
    LAT = float(LAT)*CF
    LON = float(LON)*CF
    DOY = float(DOY)

    # ---> compute declination angle
    xx = 278.97 + 0.9856*DOY + 1.9165*np.sin((356.6 + 0.9856*DOY)*CF)
    DECL = np.arcsin(0.39785*np.sin(xx*CF))
    del xx

    # --- compute day length and sunrise/sunset times.
    # Consider period when sun is above horizon, i.e. neglect civil twilight conditions   
    cosZEN = 0.0
    dl = 2.0*np.arccos(cosZEN - np.sin(LAT)*np.sin(DECL) / (np.cos(LAT)*np.cos(DECL))) / CF / 15.0  # hours

    return dl

def cloudcover(DOY, ZEN, Rglob, H2O):
    """
    Computes retrospective estimate on cloud cover fraction from surface observations
    INPUT: (all scalars or all Nx1 arrays) 
        DOY: julian day
        ZEN: sun zenith angle (rad)
        Rglob: measured global radiation (Wm-2) above canopy
        H2O: atmospheric vapor pressure (hPa)
    OUTPUT: (all scalars or all Nx1 arrays)    
        f_cloud - cloud cover fraction (0-1) [-]
        f_diff - fraction of diffuse to total radiation [-]
        emi_sky - atmospheric emissivity [-]
    SOURCE:
        Cloudiness estimate is based on Song et al. 2009 JGR 114, 2009, Appendix A & C.
        Clear-sky emissivity as in Niemelä et al. 2001 Atm. Res 58: 1-18., eq. 18 and cloud correction as Maykut & Church 1973.
        Reasonable fit against measured vs. estimated emi_sky in Hyytiälä data (tested 20.6.13)
    AUTHOR: 
        Samuli Launiainen, METLA, 2011-2013
    """        
    
            
    if np.size(DOY) == 1: #scalar input
        #inputs        
        DOY = float(DOY)
        ZEN = float(ZEN)
        Rglob = float(max(0, Rglob))
        H2O = float(H2O)
        
        Rclear = SOLAR_CONST*(1.0 + 0.033*np.cos(2*np.pi*(DOY-10) / 365))*np.cos(ZEN)  # clear sky Global radiation at surface
        tau_atm = max(0.0, Rglob / Rclear)  # atm transmissivity
        f_cloud = max(0.0, 1.0 - tau_atm / 0.7)  # cloud fraction
            
        if tau_atm > 0.7:  # fraction of diffuse rad
            f_diff = 0.2
        elif tau_atm >= 0.3 and tau_atm < 0.7:
            f_diff = 1.0-2*(tau_atm - 0.3)
        else:
            f_diff = 1.0
        
        if H2O > 2.0: #  Clear-sky atmospheric emissivity
            emi_clear = 0.72 + 0.009*(H2O - 2.0)
        else:
            emi_clear = 0.72 - 0.076*(H2O - 2.0)

        # all-sky emissivity (cloud-corrections)
        emi_sky = (1.0 + 0.22*f_cloud**2.75)*emi_clear  # Maykut & Church (1973)

    if np.size(DOY) > 1:  # array inputs 
        # inputs        
        DOY = np.array(DOY)
        ZEN = np.array(ZEN)
        Rglob = np.array(Rglob)
        H2O = np.array(H2O)

        # initialize local variables
        Rclear = np.zeros(np.size(Rglob))
        tau_atm = np.zeros(np.size(Rglob))       
        f_cloud = np.zeros(np.size(Rglob))
        f_diff = np.ones(np.size(Rglob))
        emi_clear = np.zeros(np.size(Rglob))
        emi_sky = np.zeros(np.size(Rglob))
        
        Rclear=SOLAR_CONST*(1.0 + 0.033*np.cos(2*np.pi*(DOY-10) / 365))*np.cos(ZEN)  # clear sky Global radiation at surface
        tau_atm = Rglob / Rclear  # atm transmissivity
        tau_atm[tau_atm < 0] = 0.0
        
        f_cloud = 1.0 - tau_atm / 0.7  # cloud cover fraction
        f_cloud[f_cloud < 0] = 0

        # compute  fraction of diffuse to total Global radiation: Song et al. 2009 JGR eq. A17.
        f_diff = 1.0 - 2.0*(tau_atm - 0.3)
        f_diff[tau_atm > 0.7] = 0.2
        f_diff[tau_atm < 0.3] = 1.0  # all diffuse

        #  Clear-sky atmospheric emissivity
        ind = H2O >= 2.
        emi_clear[ind] = 0.72 + 0.009*(H2O[ind] - 2.0)
        ind = H2O < 2.
        emi_clear[ind] = 0.72 - 0.076*(H2O[ind] - 2.0)
        del ind
        # print tau_atm
        # print emi_clear

        # all-sky emissivity (cloud-corrections)
        emi_sky = (1.0 + 0.22*f_cloud**2.75)*emi_clear  # Maykut & Church (1973)
        emi_sky[emi_sky > 1] = 0.98          

        # alternative formulations for all-sky emissivity
        # emi_sky=(1.0 + 0.2*f_cloud)*emi_clear
        # emi_sky=(1.0 + (1.0/emi_clear - 1.0)*0.87*f_cloud**3.49)*emi_sky # Niemelä et al. 2001 eq. 15 assuming Ts = Ta and surface emissivity = 1        
        # emi_sky=(1 - 0.84*f_cloud)*emi_sky + 0.84*f_cloud # Unsworth & Monteith (1975)

    return f_cloud, f_diff, emi_sky
        


def kbeam(ZEN, x=1.0):
    """
    COMPUTES BEAM ATTENUATION COEFFICIENT Kb (-) for given solar zenith angle ZEN (rad) and leaf angle distribution x (-)
    IN: 
        ZEN (rad): solar zenith angle 
        x (-): leaf-angle distr. parameter (optional)
                x = 1 : spherical leaf-angle distr. (default)
                x = 0 : vertical leaf-angle distr.
                x = inf. : horizontal leaf-angle distr
    OUT:
        Kb (-): beam attenuation coefficient (array of size(ZEN))
    SOURCE: Campbell & Norman (1998), Introduction to environmental biophysics
    """

    ZEN = np.array(ZEN)
    x = np.array(x)

    XN1 = (np.sqrt(x*x + np.tan(ZEN)**2))
    XD1 = (x + 1.774*(x + 1.182)**(-0.733))
    Kb = XN1 / XD1  # direct beam

    Kb = np.minimum(15, Kb)
    # if Kb<0:
    #     Kb=15

    return Kb


def kdiffuse(LAI, x=1.0):
    """
    COMPUTES DIFFUSE ATTENUATION COEFFICIENT Kd (-) BY INTEGRATING Kd OVER HEMISPHERE
    IN:
        LAI - stand leaf (or plant) are index (m2 m-2)
        x (-): leaf-angle distr. parameter (optional)
                x = 1 : spherical leaf-angle distr. (default)
                x = 0 : vertical leaf-angle distr.
                x = inf. : horizontal leaf-angle distr.
    OUT:
        Kd (-): diffuse attenuation coefficient
    USES: 
        kbeam(ZEN, x) for computing beam attenuation coeff.
    SOURCE: 
        Campbell & Norman, Introduction to environmental biophysics (1998, eq. 15.5)
    """

    LAI = float(LAI)
    x = np.array(x)

    ang = np.linspace(0, np.pi / 2, 90)  # zenith angles at 1deg intervals
    dang = ang[1]-ang[0]

    # beam attenuation coefficient - call kbeam
    Kb = kbeam(ang, x)
    # print(Kb)

    # integrate over hemisphere to get Kd, Campbell & Norman (1998, eq. 15.5)
    YY = np.exp(-Kb*LAI)*np.sin(ang)*np.cos(ang)

    Taud = 2.0*np.trapz(YY*dang)
    Kd = -np.log(Taud) / LAI  # extinction coefficient for diffuse radiation

    return Kd


def canopy_sw_ZhaoQualls(LAIz, Clump, x, ZEN, IbSky, IdSky, LeafAlbedo, SoilAlbedo, PlotFigs=False):
    """
    Computes incident (Wm-2 ground) SW radiation and absorbed (Wm-2 (leaf) radiation within canopies.
    INPUT:
        LAIz: layewise one-sided leaf-area index (m2m-2)
        Clump: element clumping index (0...1)
        x: param. of leaf angle distribution 
            (1=spherical, 0=vertical, inf.=horizontal) (-)
        ZEN: solar zenith angle (rad),scalar
        IbSky: incident beam radiation above canopy (Wm-2),scalar
        IdSky: downwelling diffuse radiation above canopy (Wm-2),scalar
        LAI: leaf (or plant-area) index, 1-sided (m2m-2),scalar
        LeafAlbedo: leaf albedo of desired waveband (-)
        SoilAlbedo: soil albedo of desired waveband (-)
        PlotFigs="True" plots figures.
    OUTPUT:
        SWbo: direct SW at z (Wm-2(ground))
        SWdo: downwelling diffuse at z (Wm-2(ground))
        SWuo: upwelling diffuse at z (Wm-2(ground))
        Q_sl: incident SW normal to sunlit leaves (Wm-2)
        Q_sh: incident SW normal to shaded leaves (Wm-2)
        q_sl: absorbed SW by sunlit leaves (Wm-2(leaf))
        q_sh: absorbed SW by shaded leaves (Wm-2(leaf))
        q_soil: absorbed SW by soil surface (Wm-2(ground))
        f_sl: sunlit fraction of leaves (-): Note: to get sunlit fraction below
            all vegetation f_sl[0] / Clump
        alb: canopy albedo (-)
    USES:
        kbeam(ZEN,x) & kdiffuse(LAI,x=1) for computing beam and diffuse attenuation coeff
    SOURCE:
        Zhao W. & Qualls R.J. (2005). A multiple-layer canopy scattering model
        to simulate shortwave radiation distribution within a homogenous plant 
         canopy. Water Resources Res. 41, W08409, 1-16.
    NOTE:
        At least for conifers NIR LeafAlbedo has to be decreased from leaf-scale  values to correctly model canopy albedo of clumped canopies.
        Adjustment from ~0.7 to 0.55 seems to be sufficient. This corresponds roughlty to 
        a=a_needle*[4*STAR / (1- a_needle*(1-4*STAR))], where a_needle is needle albedo 
        and STAR silhouette to total area ratio of a conifer shoot. STAR ~0.09-0.21 (mean 0.14) 
        for Scots pine (Smolander, Stenberg et al. papers)
    CODE: 
    Samuli Launiainen, Luke. Converted from Matlab & tested 15.5.2017
    """
    # --- check inputs and create local variables
    IbSky = max(IbSky, 0.0001)
    IdSky = max(IdSky, 0.0001)

    # original and computational grid
    LAI = Clump*sum(LAIz)  # effective LAI, corrected for clumping (m2 m-2)
    Lo = Clump*LAIz  # effective layerwise LAI (or PAI) in original grid

    Lcumo = np.cumsum(np.flipud(Lo), 0)  # cumulative plant area index from canopy top
    Lcumo = np.flipud(Lcumo)  # node 0 is canopy bottom, N is top

    # --- create computation grid
    N = np.size(Lo)  # nr of original layers
    M = np.minimum(100, N)  # nr of comp. layers
    L = np.ones([M+2])*LAI / M  # effective leaf-area density (m2m-3)
    L[0] = 0.
    L[M + 1] = 0.
    Lcum = np.cumsum(np.flipud(L), 0)  # cumulative plant area from top

    # ---- optical parameters
    aL = np.ones([M+2])*(1 - LeafAlbedo)  # leaf absorptivity
    tL = np.ones([M+2])*0.4  # transmission as fraction of scattered radiation
    rL = np.ones([M+2])*0.6  # reflection as fraction of scattered radiation

    # soil surface, no transmission
    aL[0] = 1. - SoilAlbedo
    tL[0] = 0.
    rL[0] = 1.

    # upper boundary = atm. is transparent for SW
    aL[M+1] = 0.
    tL[M+1] = 1.
    rL[M+1] = 0.

    # Compute black leaf extinction coefficients for direct beam and diffuse radiation
    Kb = kbeam(ZEN, x)
    Kd = kdiffuse(LAI, x)

    # fraction of sunlit & shad ground area in a layer (-)
    f_sl = np.flipud(np.exp(-Kb*Lcum))

    # beam radiation at each layer
    Ib = f_sl*IbSky

    # propability for beam and diffuse penetration through layer without interception
    taub = np.zeros([M+2])
    taud = np.zeros([M+2])
    taub[0:M+2] = np.exp(-Kb*L)
    taud[0:M+2] = np.exp(-Kd*L)

    # soil surface is non-transparent
    taub[0] = 0.
    taud[0] = 0.

    # backward-scattering functions (eq. 22-23) for beam rb and diffuse rd
    rb = np.zeros([M+2])
    rd = np.zeros([M+2])
    rb = 0.5 + 0.3334*(rL - tL) / (rL + tL)*np.cos(ZEN)
    rd = 2.0 / 3.0*rL/(rL + tL) + 1.0 / 3.0*tL / (rL + tL)

    rb[0] = 1.
    rd[0] = 1.
    rb[M+1] = 0.
    rd[M+1] = 0.

    # --- set up tridiagonal matrix A and solve SW without multiple scattering
    # from A*SW = C (Zhao & Qualls, 2006. eq. 39 & 42)
    A = np.zeros([2*M+2, 2*M+2])
  
    # lowermost rows: 0 = soil surface, M+2 is upper boundary
    A[0, 0] = 1.

    A[1, 0] = - (taud[1] + (1 - taud[1])*(1 - aL[1])*(1 - rd[1]))
    A[1, 1] = - 1*(taud[1] + (1 - taud[1])*(1 - aL[1])*(1 - rd[1]))*(1 - aL[0])     
    A[1, 2] = (1 - 1*rd[1]*(1 - aL[0])*1*(1 - aL[1])*(1 - taud[1])) 

    # middle rows
    for k in range(1, M + 1):
        A[2*k-1, 2*k-2] = - (taud[k] + (1 - taud[k])*(1 - aL[k])*(1 - rd[k]))
        A[2*k-1, 2*k-1] = - rd[k - 1]*(taud[k] + (1 - taud[k])*(1 - aL[k])*(1 - rd[k]))*(1 - aL[k-1])*(1 - taud[k-1])
        A[2*k-1, 2*k] = (1 - rd[k-1]*rd[k]*(1 - aL[k-1])*(1 - taud[k-1])*(1 - aL[k])*(1 - taud[k]))

        A[2*k, 2*k-1] = (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        A[2*k, 2*k] = - rd[k+1]*(taud[k] + (1 - taud[k])*(1 - aL[k])*(1 - rd[k]))*(1 - aL[k+1])*(1 - taud[k+1])
        A[2*k, 2*k+1] = -(taud[k] + (1 - taud[k])*(1 - aL[k])*(1 - rd[k]))

    # uppermost node2*M+2
    A[2*M+1, 2*M+1] = 1.
    del k
    # print A

    # --- RHS vector C
    C = np.zeros([2*M+2, 1])

    # lowermost row
    C[0] = SoilAlbedo*Ib[0]
    n = 1  # dummy
    for k in range(1, M+1):  # k=2:M-1,
        C[n] = (1 - rd[k-1]*rd[k]*(1 - aL[k-1])*(1 - taud[k-1])*(1 - aL[k])*(1 - taud[k]) )*rb[k]*(1 - taub[k])*(1 - aL[k])*Ib[k]
        C[n+1] = (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))*(1 - taub[k])*(1 - aL[k])*(1 - rb[k])*Ib[k] 
        # Ib(k+1) Ib(k):n sijaan koska tarvitaan kerrokseen tuleva
        n = n + 2

    # uppermost row
    C[2*M+1] = IdSky

    # ---- solve A*SW = C
    SW = np.linalg.solve(A, C)
    
    # upward and downward hemispherical radiation (Wm-2 ground)
    SWu0 = SW[0:2*M+2:2]
    SWd0 = SW[1:2*M+2:2]
    del A, C, SW, k  
  
    # ---- Compute multiple scattering, Zhao & Qualls, 2005. eq. 24 & 25.
    # downwelling diffuse after multiple scattering, eq. 24
    SWd = np.zeros([M+1])
    for k in range(M-1, -1, -1):  # downwards from layer k+1 to layer k
        X = SWd0[k+1] / (1 - rd[k]*rd[k+1]*(1-aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        Y = SWu0[k]*rd[k+1]*(1 - aL[k+1])*(1 - taud[k+1]) / (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        SWd[k+1] = X + Y
    SWd[0] = SWd[1] # SWd0[0]
    # print SWd  

    # upwelling diffuse after multiple scattering, eq. 25
    SWu = np.zeros([M+1])
    for k in range(0, M, 1):  # upwards from layer k to layer k+1
        X = SWu0[k] / (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        Y = SWd0[k+1]*rd[k]*(1 - aL[k])*(1 - taud[k]) / (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        SWu[k] = X + Y
    SWu[M] = SWu[M-1]

    # match dimensions of all vectors
    Ib = Ib[1:M+2]
    f_sl = f_sl[1:M+2]
    Lcum = np.flipud(Lcum[0:M+1])
    aL = aL[0:M+1]

    # --- NOW return values back to the original grid
    f_slo = np.exp(-Kb*Lcumo)
    SWbo = f_slo*IbSky  # Beam radiation


    # interpolate diffuse fluxes
    X = np.flipud(Lcumo)
    xi = np.flipud(Lcum)
    SWdo = np.flipud(np.interp(X, xi, np.flipud(SWd)))
    SWuo = np.flipud(np.interp(X, xi, np.flipud(SWu)))
    del X, xi

    # incident radiation on sunlit and shaded leaves Wm-2
    Q_sh = Kd*(SWdo + SWuo)  # normal to shaded leaves is all diffuse
    Q_sl = Kb*IbSky + Q_sh  # normal to sunlit leaves is direct and diffuse

    # absorbed components
    aLo = np.ones(len(Lo))*(1 - LeafAlbedo)
    aDiffo = aLo*Kd*(SWdo + SWuo)
    aDiro = aLo*Kb*IbSky

    # stand albedo
    alb = SWuo[-1] / (IbSky + IdSky + eps)
    # print alb
    # soil absorption (Wm-2 (ground))
    q_soil = (1 - SoilAlbedo)*(SWdo[0] + SWbo[0])

    # correction to match absorption-based and flux-based albedo, relative error <3% may occur
    # in daytime conditions, at nighttime relative error can be larger but fluxes are near zero
    aa = (sum(aDiffo*Lo + aDiro*f_slo*Lo) + q_soil) / (IbSky + IdSky + eps)
    F = (1. - alb) / aa
    # print('F', F)
    if F <= 0 or np.isfinite(F) is False:
        F = 1.

    aDiro = F*aDiro
    aDiffo = F*aDiffo
    q_soil = F*q_soil

    # sunlit fraction in clumped foliage; clumping means elements shade each other
    f_slo = Clump*f_slo

    # Following adjustment is required for energy conservation, i.e. total absorbed radiation
    # in a layer must equal difference between total radiation(SWup, SWdn) entering and leaving the layer. 
    # Note that this requirement is not always fullfilled in canopy rad. models.
    # now sum(q_sl*f_slo* + q_sh*(1-f_slo)*Lo = (1-alb)*(IbSky + IdSky) 
    q_sh = aDiffo*Clump  # shaded leaves only diffuse
    q_sl = q_sh + aDiro  # sunlit leaves diffuse + direct

    if PlotFigs:
        plt.figure(999)
        plt.subplot(221)
        plt.title("Source: radiation.canopy_sw_ZhaoQualls")

        plt.plot(f_slo, -Lcumo, 'r-', (1 - f_slo), -Lcumo, 'b-')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("sunlit & shaded fractions (-)")
        plt.legend(('f_{sl}', 'f_{sh}'), loc='best')

        # add input parameter values to fig
        plt.text(0.05, 0.75, r'$LAI$ = %1.1f m2 m-2' % (LAI))
        plt.text(0.05, 0.65, r'$ZEN$ = %1.3f rad' % (ZEN))
        plt.text(0.05, 0.55, r'$\alpha_l$ = %0.2f' % (LeafAlbedo))
        plt.text(0.05, 0.45, r'$\alpha_s$ = %0.2f' % (SoilAlbedo))

        plt.subplot(223)
        plt.plot(SWd, -Lcum, 'bo', SWdo, -Lcumo, 'b-', Ib, -Lcum, 'ro',
                 SWbo, -Lcumo, 'r-', SWu, -Lcum, 'go', SWuo, -Lcumo, 'g-')
        plt.legend(('SWd', 'Swdo', 'SWb', 'SWbo', 'SWu', 'SWuo'), loc='best')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("Incident SW (Wm-2 )")

        plt.subplot(224)
        plt.plot(q_sl, -Lcumo, 'ro-', q_sh, -Lcumo, 'bo-')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("Absorbed radiation (Wm-2 (leaf))")
        plt.legend(('sunlit', 'shaded'), loc='best')

    return SWbo, SWdo, SWuo, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_slo, alb
    
    
def canopy_sw_Spitters(LAIz, Clump, x, ZEN, IbSky, IdSky, LeafAlbedo, SoilAlbedo, PlotFigs="False"):
#ZEN, IbSky, IdSky, LAI, z, lad, x, CLUMP, LeafAlbedo, SoilAlbedo, PlotFigs="False"):
    """
    Computes profiles of incident and absorbed SW within plant canopies using analytic model of Spitters (1986) 
    without explicit treatment of upward and downward-scattered radiation.
    INPUT:
       LAIz: layewise one-sided leaf-area index (m2m-2)
        Clump: element clumping index (0...1)
        x: param. of leaf angle distribution 
            (1=spherical, 0=vertical, inf.=horizontal) (-)
        ZEN: solar zenith angle (rad),scalar
        IbSky: incident beam radiation above canopy (Wm-2),scalar
        IdSky: downwelling diffuse radiation above canopy (Wm-2),scalar
        LAI: leaf (or plant-area) index, 1-sided (m2m-2),scalar
        LeafAlbedo: leaf albedo of desired waveband (-)
        SoilAlbedo: soil albedo of desired waveband (-)
        PlotFigs="True" plots figures.
 
    OUTPUT:
        SWb: direct SW at z (Wm-2(ground))
        SWd: downwelling diffuse at z (Wm-2(ground))
        Q_sl: incident SW normal to sunlit leaves (Wm-2)
        Q_sh: incident SW normal to shaded leaves (Wm-2)
        q_sl: absorbed SW by sunlit leaves (Wm-2(leaf))
        q_sh: absorbed SW by shaded leaves (Wm-2(leaf))
        q_soil: absorbed SW by soil surface (Wm-2(ground))
        f_sl: sunlit fraction at z (-)
        alb: canopy albedo (-)
    USES:
        kbeam(ZEN, x) & kdiffuse(LAI, x) for computing beam and diffuse attenuation coeff.
    SOURCE:
        Attenuation coefficients and canopy reflectance based on Campbell & Norman (1998): An introduction to environmental 
        biophysics, Springer.
        Other algorithms from Spitters C.T.J. (1986): Separating the diffuse and direct component of global radiation 
        and its implications for modeling canopy photosynthesis part II: Calculation of canopy photosynthesis. 
        Agric. For. Meteorol. 38, 231-242.
    NOTE:
        At least for conifers NIR LeafAlbedo has to be decreased from leaf-scale values to correctly model canopy albedo of clumped canopies.
        Adjustment from ~0.7 to 0.55 seems to be sufficient. This corresponds to a=a_needle*[4*STAR / (1- a_needle*(1-4*STAR))], 
        where a_needle is needle albedo and STAR silhouetteto total area ratio of a conifer shoot. STAR ~0.09-0.21 (mean 0.14) for 
        Scots pine. Still then overestimates NIR absorption in upper canopy layers, compared to canopy_sw_ZhaoQualls with multiple scattering.
        Assumes isotropic scattering and does not explicitly compute upward reflected SW.        
        To compute incident total downward SW at z: SWdown=f_sl*SWbo + (1-f_sl)*SWdo       
        Arg. PlotFigs="True" plots figures.
    AUTHOR: 
        Samuli Launiainen, METLA 1/2011-4/2014 (converted to Python 16.04.2014)
    """
    # --- check inputs and create local variables
    IbSky = max(IbSky, 0.0001)
    IdSky = max(IdSky, 0.0001)

    L = Clump*LAIz  # effective layerwise LAI (or PAI) in original grid
    Lcum = np.flipud(np.cumsum(np.flipud(L), 0.0))  # cumulative plant area from the sky, index 0 = ground
    LAI = max(Lcum)    
    #N = np.size(L, 0)  # Nr of layers
    # print L,Lcum

    # attenuation coefficients
    Kb = kbeam(ZEN, x)
    Kd = kdiffuse(x, LAI)
 
    # sunlit fraction as a function of L
    f_sl = np.exp(-Kb*Lcum)

    # canopy hemispherical reflection coefficient Campbell & Norman (1998)
    rcpy1 = (1.0 - (1.0 - LeafAlbedo)**0.5) / (1.0 + (1.0 - LeafAlbedo)**0.5)

    # in case canopy is deep, soil reflections have small impact and canopy reflection coefficients becomes
    rb1 = 2.0*Kb / (Kb + 1.0)*rcpy1  # beam
    rd1 = 2.0*Kd / (Kd + 1.0)*rcpy1  # diffuse
    
    # but in sparser canopies soil reflectance has to be taken into account and this yields
    AA = ((rb1 - SoilAlbedo) / (rb1*SoilAlbedo - 1.0))*np.exp(-2.0*(1.0 - LeafAlbedo)**0.5*Kb*LAI)
    rb1 = (rb1 + AA) / (1.0 + rb1*AA)  # beam
    del AA

    AA = ((rd1 - SoilAlbedo) / (rd1*SoilAlbedo - 1.0))*np.exp(-2.0*(1.0 - LeafAlbedo)**0.5*Kd*LAI)
    rd1 = (rd1 + AA) / (1.0 + rd1*AA)  # diffuse
    del AA

    # canopy albedo
    alb = (rb1*IbSky + rd1*IdSky) / (IbSky + IdSky)

    # Incident SW as a function of Lcum

    qd1 = (1.0 - rd1)*IdSky*np.exp(-(1.0 - LeafAlbedo)**0.5*Kd*Lcum)  # attenuated diffuse
    qb1 = IbSky*np.exp(-Kb*Lcum)  # beam
    qbt1 = (1.0 - rb1)*IbSky*np.exp(-(1.0 - LeafAlbedo)**0.5*Kb*Lcum)  # total beam
    qsc1 = qbt1 - (1.0 - rb1)*qb1  # scattered part of beam
    # print Lcum, f_sl, qd1, qb1, qsc1

    # incident fluxes at each layer per unit ground area
    SWd = qd1 + qsc1  # total diffuse
    SWb = qb1  # total direct beam

    # incident to leaf surfaces (per m2 (leaf))
    Q_sh = Kd*SWd
    Q_sl = Q_sh + Kb*IbSky

    # absorbed components: A = -dq/dL (Wm-2 (leaf))
    # diffuse, total beam, direct beam
    Ad1 = (1.0 - rd1)*IdSky*(1.0 - LeafAlbedo)**0.5*Kd*np.exp(-(1.0 - LeafAlbedo)**0.5*Kd*Lcum)
    Abt1 = (1.0 - rb1)*IbSky*(1.0 - LeafAlbedo)**0.5*Kb*np.exp(-(1.0 - LeafAlbedo)**0.5*Kb*Lcum)
    Ab1 = (1.0 - rb1)*(1.0 - LeafAlbedo)*IbSky*Kb*np.exp(-(1.0 - LeafAlbedo)**0.5*Kb*Lcum)

    # absorbed at sunlit & shaded leaves (Wm-2(leaf))
    q_sh = Ad1 + (Abt1 - Ab1)  # total absorbed diffuse
    q_sl = q_sh + (1.0 - LeafAlbedo)*Kb*IbSky  # sunlit leaves recieve additional direct radiation, Spitters eq. 14

    # absorbed by soil surface (Wm-2(ground))
    q_soil = (1.0 - SoilAlbedo)*(SWb[-1] + SWd[-1])

    if PlotFigs == "True":         
        # print figures
        plt.figure(888)
        plt.title("Source: radiation.canopy_sw_Spitters")

        plt.subplot(121)
        plt.plot(SWb, -Lcum, 'ro-', SWd, -Lcum, 'bo-', qsc1, -Lcum, 'go-')
        plt.legend(('SWtotal', 'SWbeam', 'SWdiff', 'SWscatt'), loc='best')
        plt.ylabel("-Lcum")
        plt.xlabel("Incident SW (Wm-2 )")
        plt.title(r'$\alpha_{can}$ = %1.2f' % alb)

        # print parameter values to fig
        plt.text(0.05, -0.5, r'$LAI$ = %1.1f m2 m-2' % LAI)
        plt.text(0.05, -0.65, r'$ZEN$ = %1.3f rad' % ZEN)
        plt.text(0.05, -0.8, r'$\alpha_l$ = %0.2f' % LeafAlbedo)
        plt.text(0.05, -0.95, r'$\alpha_s$ = %0.2f' % SoilAlbedo)

        plt.subplot(122)
        plt.plot(q_sl, -Lcum, 'ro-', q_sh, -Lcum, 'bo-')
        plt.ylabel("-Lcum")
        plt.xlabel("Absorbed radiation (Wm-2 (leaf))")
        plt.legend(('sunlit', 'shaded'), loc='best')

    return SWb, SWd, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_sl, alb

 
def canopy_lw(Lz, Clump, T, LWdn0, LWup0, leaf_emi=1.0):
    """
    Estimates long-wave radiation budget and net isothermal LW radiation within the canopy
    Assumes canopy elements as black bodies (es=1.0) at local air temperature
    T(z). Neglects scattering etc.
    
    INPUT: 
       Lz - 1-sided LAI in each layer (m2m2(ground))
       CLUMP - clumping factor (-), [0...1]
       T: air temperature (K)
       LWdn0 - downwelling LW above uppermost gridpoint (Wm-2(ground)). LWdn0=eatm*b*Tatm^4
       LWup0 - upwelling LW at surface (Wm-2(ground)). LWup0=esurf*b*Tsurf^4
    OUTPUT:
       LWleaf - net isothermal long-wave absorbed by leaf (Wm-2 (leaf))
       LWdn - downwelling LW (Wm-2)
       LWup - upwelling LW (Wm-2)
    SOURCE:
       Modified from Flerchinger et al. 2009. Simulation of within-canopy radiation exchange, NJAS 57, 5-15.
       Assumes Tcan(z) = T(z) and neglects scattering (canopy elements black bodies)
    
    Samuli Launiainen (Luke). Last edit: 12.5.2017
    """

    if np.min(T) < 200:
        T = T + 273.15

    N = len(Lz)  # Node 0 = ground
    LWdn = np.zeros(N)
    LWup = np.zeros(N)

    LayerLAI = Clump*Lz  # plant-area m2m-2 in a layer addjusted for clumping
    cantop = max(np.where(LayerLAI>0)[0])  # node at canopy top
    
    # layerwise attenuation coeffcient
    Kd = 0.7815  # diffuse radiation, Flearchinger et al. 2009 NJAS
    tau = np.exp(-Kd*LayerLAI)

    # LW down
    LWdn[cantop+1:N] = LWdn0 
    for k in range(cantop, -1, -1):
        LWdn[k]=tau[k]*LWdn[k+1] +(1 - tau[k])*(leaf_emi*SIGMA*T[k]**4)
    del k

    # LW up
    LWup[0] = LWup0
    for k in range(1, cantop+1):
        LWup[k] = tau[k]*LWup[k-1] + (1 - tau[k])*(leaf_emi*SIGMA*T[k]**4)
    del k
    LWup[cantop+1:N] = LWup[cantop]

    # absorbed isothermal net radiation by the leaf (Wm-2(leaf))
    # Kd is mean projection of leaves 
    LWleaf = Kd*(1 - tau)*leaf_emi*(LWdn + LWup - 2*SIGMA*T**4)/ (Lz + eps)

    return LWleaf, LWdn, LWup
    
def canopy_lw_ZhaoQualls(LAIz, Clump, x, Tleaf, LWdn0, LWup0, leaf_emi=0.98, soil_emi=0.98, PlotFigs=False):
    """
    Estimates incident and absorbed LWradiation within plant canopies.
    Includes multiple scattering among the canopy layers and soil surface.
    INPUT:
      LAIz: leaf-area index per layer (m2(leaf) m-2(ground)), Nx1-array. LAIz[-1] MUST be 0!
      Clump: element clumping index (-), ]0..1[, scalar
      leaf_angle_para: leaf-angle distribution parameter (-), scalar: 1=spherical, 0=vertical, '1000'=horizontal
      Tleaf: leaf (or air) temperature (degC)
      LWdn0: downwelling LW at canopy top (Wm-2)
      LWup0: upwelling LW at soil surface (Wm-2)
      leaf_emi: leaf emissivity (-)
      soil_emi: soil surface emissivity (-)
    OUTPUT:
      LWleaf: longwave radiation absorbed per unit un-clumped leaf-area in a canopy layer(Wm-2(leaf))
      LWdn: downward long-wave (Wm-2(ground))
      LWup: upward long-wave (Wm-2(ground))
    USES:
      kdiffuse
    REFERENCES:
      Zhao & Qualls, 2005 Water Resources Research. Multi-layer Multiple scattering model
      Estimates SW attenuation inc. multiple scattering within the canopy.
      Campbell & Norman (1998): An introduction to environmental biophysics.
      Spitters, (1986): Separating the diffuse and direct component of global
      radiation and its implications for modeling canopy photosynthesis part
      II: Calculation of canopy photosynthesis. Agric. For. Meteorol. 38, 231-242.
      Wang & Leuning (1998): A two-leaf model for canopy conductance, photosynthesis and partitioning of
      available energy I: Model description and comparison with a multi-layered model. Agric. For. Meteorol. 91, 89-111.
      Juang, Katul et al. 2008: Investigating a hierarchy of Eulerian Closure Models for Scalar Transfer Inside Forested Canopies. BLM 128:1-32.
      Flerchinger et al. 2009. Simulation of within-canopy radiation exchange, NJAS 57, 5-15

    Samuli Launiainen (METLA) 2013-2015
    
    CONVERTED TO PYTHON / LAST EDIT: 15.5.2017
    """
    Tleaf += 273.15  # K

    # original and computational grid
    LAI = Clump*sum(LAIz)  # effective LAI, corrected for clumping (m2 m-2)
    Lo = Clump*LAIz + eps  # effective layerwise LAI (or PAI) in original grid

    Lcumo = np.cumsum(np.flipud(Lo), 0)  # cumulative plant area index from canopy top
    Lcumo = np.flipud(Lcumo)  # node 0 is canopy bottom, N is top

    # --- create computation grid
    N = np.size(Lo)  # nr of original layers
    M = np.maximum(10, N)  # nr of comp. layers

    L = np.ones([M+2])*LAI / M  # effective leaf-area density (m2m-3)
    L[0] = 0.
    L[M + 1] = 0.
    Lcum = np.cumsum(np.flipud(L), 0)  # cumulative plant area from top
    # interpolate T to comp. grid
    T = np.flipud(np.interp(Lcum, Lcumo, Tleaf))

    # ---- optical parameters
    # back-scattering fraction, approximation, eq. (6-7)
    if x == 1:  # spherical leaf distrib.
        rd = 2./3.
    elif x == 0:  # vertical leafs
        rd = 0.5
    elif x > 100:  # horizontal leafs
        rd = 1.
    else:
        print "radiation.canopy_lw_ZhaoQualls: check leaf angle distr. "

    rd = np.ones([M+2])*rd
    rd[0] = 1.
    rd[M+1] = 0.

    aL = np.ones([M+2])*leaf_emi  # leaf emissivity
    aL[0] = soil_emi
    # aL[M+1] = 0.

    # extinction coefficients for diffuse radiation
    Kd = kdiffuse(LAI, x)

    # propability of contact with canopy elements in each layer
    taud = np.exp(-Kd*L)  # diffuse
    taud[0] = 0.
    taud[M+1] = 1.

    # --- set up tridiagonal matrix A and solve LW without multiple scattering.
    # Zhao & Qualls eq's. (16 - 25)
    A = np.zeros([2*M+2, 2*M+2])
    
    #layer 0 = soil surface, M+1 is upper boundary 
    A[0, 0] = 1.
    
    # middle rows
    for k in range(1, M+1):
        A[2*k-1, 2*k-2] = - (taud[k] + (1 - taud[k])*(1 - aL[k])*(1 - rd[k]))
        A[2*k-1, 2*k-1] = - rd[k-1]*(taud[k] + (1 - taud[k])*(1 - aL[k])*(1 - rd[k]))*(1 - aL[k-1])*(1 - taud[k-1])
        A[2*k-1, 2*k] = (1 - rd[k-1]*rd[k]*(1 - aL[k-1])*(1 - taud[k-1])*(1 - aL[k])*(1 - taud[k]))

        A[2*k, 2*k-1] = (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        A[2*k, 2*k] = - rd[k+1]*(taud[k] + (1 - taud[k])*(1 - aL[k])*(1 - rd[k]))*(1 - aL[k+1])*(1 - taud[k+1])
        A[2*k, 2*k+1] = -(taud[k] + (1 - taud[k])*(1 - aL[k])*(1 - rd[k]))

    # uppermost node2*M+2
    A[2*M+1, 2*M+1] = 1.
    del k

    # --- RHS vector C
    C = np.zeros([2*M+2, 1])

    LWsource = aL*SIGMA*T**4

    # lowermost row
    C[0] = LWup0
    n = 1  # dummy
    for k in range(1, M+1):  # k=2:M-1,
        C[n] = (1 - rd[k-1]*rd[k]*(1 - aL[k-1])*(1 - taud[k-1])*(1 - aL[k]))*(1 - taud[k]) *LWsource[k]
        C[n+1] = (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1]))*(1 - taud[k+1])*LWsource[k]
        # Ib(k+1) Ib(k):n sijaan koska tarvitaan kerrokseen tuleva
        n = n + 2

    # uppermost rows
    C[2*M] = (1 - taud[M])*LWsource[M]
    C[2*M+1] = LWdn0

    # ---- solve A*SW = C
    LW = np.linalg.solve(A, C)

    # upward and downward hemispherical radiation (Wm-2 ground)
    LWu0 = LW[0:2*M+2:2]
    LWd0 = LW[1:2*M+2:2]
    del A, C, LW

    # ---- Compute multiple scattering, Zhao & Qualls, 2005. eq. (8 & 9)
    # downwelling diffuse after multiple scattering
    LWd = np.zeros([M+1])
    for k in range(M-1, -1, -1):  # downwards from layer k+1 to layer k
        X = LWd0[k+1] / (1 - rd[k]*rd[k+1]*(1-aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        Y = LWu0[k]*rd[k+1]*(1 - aL[k+1])*(1 - taud[k+1]) / (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        LWd[k+1] = X + Y
    LWd[0] = LWd0[0]

    # upwelling diffuse after multiple scattering
    LWu = np.zeros([M+1])
    for k in range(0, M, 1):  # upwards from layer k to layer k+1 
        X = LWu0[k] / (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        Y = LWd0[k+1]*rd[k]*(1 - aL[k])*(1 - taud[k]) / (1 - rd[k]*rd[k+1]*(1 - aL[k])*(1 - taud[k])*(1 - aL[k+1])*(1 - taud[k+1]))
        LWu[k] = X + Y
    LWu[M] = LWu[M-1]

    # --- NOW return values back to the original grid
    Lcum = Lcum[0:M+1]
    LWd = np.flipud(LWd)
    LWu = np.flipud(LWu)
    # plt.figure(99); plt.plot(LWd, -Lcum, 'r', LWu, -Lcum, 'g')

    X = np.flipud(Lcumo)
    xi = Lcum
    LWdn = np.interp(X, xi, LWd)
    LWup = np.interp(X, xi, LWu)
#    del X, xi

    #---------------------------------------------------------------------
    # check that interpolation is ok
    #plt.figure(100)
    #plt.plot(LWd,-xi,'r.-',LWdn,-X,'ro', LWu,-xi,'b-',LWup,-X,'bs')
    #plt.title('r = dn, b = up' ); plt.ylabel('LAI eff (m2m-2)')
    #---------------------------------------------------------------------

    # absorbed net LW per unit un-clumped leaf area (Wm-2(leaf)),
    # Flerchinger et al. 2009. NJAS 57, 5-15

    taud = np.exp(-Kd*Lo)
    eL = np.ones(len(taud))*leaf_emi
    # Laiz[-1] must be zero -> Lo[-1] is zero
    LWleaf = Clump*(1 - taud)*eL*(LWdn + LWup - 2*SIGMA*Tleaf**4) / Lo
    LWleaf[Lo==0] = 0.0
    if PlotFigs:
        plt.figure()
        plt.subplot(121)
        plt.title("radiation.canopy_lw_ZhaoQualls", fontsize=8)
    
        plt.plot(LWdn, -X, 'bo', label='LWdn')
        plt.plot(LWup, -X, 'ro', label='LWup')    
        plt.ylabel("-Lcum eff.")
        plt.xlabel("LW (Wm-2 )")
        plt.legend()
    
        plt.subplot(122)
        plt.plot(LWdn-LWup, -X, 'go',)
        plt.ylabel("-Lcum eff.")
        plt.xlabel("LWdn - LWup (Wm-2 )")
        plt.title('LWup0=%.1f, LWdn0=%.1f' % (LWup0, LWdn0))
    return LWleaf,  LWdn, LWup
