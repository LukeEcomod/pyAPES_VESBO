# -*- coding: utf-8 -*-
"""
.. module: radiation
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes distribution of within canopy radiation.
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Oct 02 09:04:05 2018

Note:
    migrated to python3
    - absolute import
    - added parenthesis to print
    - list(range(...)) forward-compatible: from builtins import range, import can be removed later

References:
Launiainen, S., Katul, G.G., Lauren, A. and Kolari, P., 2015. Coupling boreal
forest CO2, H2O and energy flows by a vertically structured forest canopy –
Soil model with separate bryophyte layer. Ecological modelling, 312, pp.385-405.
"""
from builtins import range
import numpy as np
import pandas as pd
import logging
from tools.utilities import tridiag
from matplotlib import pyplot as plt
from .constants import DEG_TO_RAD, DEG_TO_KELVIN, STEFAN_BOLTZMANN, SPECIFIC_HEAT_AIR, EPS
logger = logging.getLogger(__name__)

class Radiation(object):
    r""" Describes distribution of within canopy radiation.
    """
    def __init__(self, p, Ebal):
        """
        Args:
            p (dict):
                'clump': clumping index [-]
                'leaf_angle': leaf-angle distribution [-]
                'Par_alb': shoot Par-albedo [-]
                'Nir_alb': shoot NIR-albedo [-]
                'leaf_emi': 0.98
            Ebal (bool): True solves LW radiation
        Returns:
            self (object)
        """

        # parameters
        self.clump = p['clump']
        self.leaf_angle = p['leaf_angle']  # leaf-angle distribution [-]
        self.alb = {'PAR': p['Par_alb'],  # shoot Par-albedo [-]
                    'NIR': p['Nir_alb']}  # shoot Nir-albedo [-]
        self.leaf_emi = p['leaf_emi']

        # model functions to use
        self.SWmodel = 'ZHAOQUALLS'
        self.LWmodel = 'ZHAOQUALLS'

        logger.info('Shortwave radiation model: %s', self.SWmodel)
        if Ebal:
            logger.info('Longwave radiation model: %s', self.LWmodel)

    def shortwave_profiles(self, forcing, parameters):
        r""" Computes distribution of within canopy shortwave radiation
        using specified model.

        Args:
            forcing (dict):
                'zenith_angle': solar zenith angle [rad]
                'dir_par'/'dir_nir': direct Par [Wm-2]
                'diffSW': diffuse Par [Wm-2]
            parameters (dict):
                'LAIz': layewise one-sided leaf-area index [m2m-2]
                'radiation_type': 'NIR' or 'PAR'

        Returns:
            Q_sl (array): incident SW normal to sunlit leaves [W m-2]
            Q_sh (array): incident SW normal to shaded leaves [W m-2]
            q_sl (array): absorbed SW by sunlit leaves [W m-2(leaf)]
            q_sh (array): absorbed SW by shaded leaves [W m-2(leaf)]
            q_soil (array): absorbed SW by soil surface [W m-2(ground)]
            f_sl (array): sunlit fraction of leaves [-]: Note: to get sunlit fraction below
                all vegetation f_sl[0] / clump
            SW_gr (float): incident SW at ground level [W m-2]

        References:
            Zhao W. & Qualls R.J. (2005). A multiple-layer canopy scattering model
            to simulate shortwave radiation distribution within a homogenous plant
            canopy. Water Resources Res. 41, W08409, 1-16.
        """

        radtype = parameters['radiation_type'].upper()
        if radtype == 'PAR' or radtype == 'NIR':

            if self.SWmodel == 'ZHAOQUALLS':
                SWb, SWd, SWu, Q_sl, Q_sh, q_sl, q_sh, q_gr, f_sl, alb = canopy_sw_ZhaoQualls(
                    parameters['LAIz'],
                    self.clump, self.leaf_angle,
                    forcing['zenith_angle'],
                    forcing[radtype]['direct'],
                    forcing[radtype]['diffuse'],
                    self.alb[radtype],
                    parameters['ff_albedo'][radtype])

            results = {'sunlit':{'incident': Q_sl, 'absorbed': q_sl, 'fraction': f_sl},
                       'shaded':{'incident': Q_sh, 'absorbed': q_sh},
                       'ground': SWb[0] + SWd[0],
                       'up': SWu,
                       'down': SWb + SWd}

            # ADD OTHER MODELS??

            return results

        else:
            raise ValueError("Radiation type is not 'PAR' or 'NIR'")

    def longwave_profiles(self, forcing, parameters):
        r""" Computes distribution of within canopy longwave radiation
        using specified model.

        Args:
            forcing (dict):
                'leaf_temperature' (array): leaf temperature [degC]
                'lw_in' (float): downwelling longwave raditiona above uppermost gridpoint [W m-2(ground)]
                'lw_up' (float): upwelling longwave radiation below canopy (forest floor) [W m-2(ground)]
            parameters (dict):
                'LAIz' (array): layewise one-sided leaf-area index [m2 m-2]
                'ff_emissivity' (float):
        Returns:
            profiles (dict):
                'net_leaf_lw' (array): leaf net longwave radiation [W m-2(leaf)]
                'lw_dn' (array): downwelling LW [W m-2]
                'lw_up' (array): upwelling LW [W m-2]
                'radiative_conductance' (array): radiative conductance [mol m-2 s-1]

        References:
            Flerchinger et al. 2009. Simulation of within-canopy radiation exchange,
            NJAS 57, 5-15.
            Zhao, W. and Qualls, R.J., 2006. Modeling of long‐wave and net radiation
            energy distribution within a homogeneous plant canopy via multiple scattering
            processes. Water resources research, 42(8).
        """


        if self.LWmodel == 'FLERCHINGER':
            lw_leaf, lw_dn, lw_up, gr = canopy_lw(
                LAIz=parameters['LAIz'],
                Clump=self.clump,
                x=self.leaf_angle,
                T=forcing['leaf_temperature'],
                LWdn0=forcing['lw_in'],
                LWup0=forcing['lw_up'],
                leaf_emi=self.leaf_emi
            )

        if self.LWmodel == 'ZHAOQUALLS':
            lw_leaf, lw_dn, lw_up, gr = canopy_lw_ZhaoQualls(
                LAIz=parameters['LAIz'],
                Clump=self.clump,
                x=self.leaf_angle,
                Tleaf=forcing['leaf_temperature'],
                LWdn0=forcing['lw_in'],
                LWup0=forcing['lw_up'],
                leaf_emi=self.leaf_emi,
                soil_emi=parameters['ff_emissivity']
            )

        # ADD OTHER MODELS??

        results = {
            'net_leaf': lw_leaf,
            'down': lw_dn,
            'up': lw_up,
            'radiative_conductance': gr
        }

        return results

"""
stand-alone functions start here: these can be called with arguments only
"""

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
    ha0 = np.arccos(aa) / DEG_TO_RAD

    sunrise = 720.0 - 4.*(lon + ha0) - et  # minutes
    sunset = 720.0 - 4.*(lon - ha0) - et  # minutes

    daylength = (sunset - sunrise)  # minutes

    sunrise = sunrise + timezone
    sunrise[sunrise < 0] = sunrise[sunrise < 0] + 1440.0

    sunset = sunset + timezone
    sunset[sunset > 1440] = sunset[sunset > 1440] - 1440.0

    return zen, azim, decl, sunrise, sunset, daylength

def kbeam(ZEN, x=1.0):
    """
    COMPUTES BEAM ATTENUATION COEFFICIENT Kb (-) for given solar zenith angle
    ZEN (rad) and leaf angle distribution x (-)
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
    Kd = -np.log(Taud) / (LAI + EPS)  # extinction coefficient for diffuse radiation

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
        At least for conifers NIR LeafAlbedo has to be decreased from leaf-scale  values to
        correctly model canopy albedo of clumped canopies.
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

    # black leaf extinction coefficients for direct beam and diffuse radiation
    Kb = kbeam(ZEN, x)
    Kd = kdiffuse(LAI, x)

    # fraction of sunlit & shad ground area in a layer (-)
    f_sl = np.flipud(np.exp(-Kb*(Lcum)))

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
    f_slo = np.exp(-Kb*(Lcumo))
    SWbo = f_slo*IbSky  # Beam radiation


    # interpolate diffuse fluxes
    X = np.flipud(Lcumo)
    xi = np.flipud(Lcum)
    SWdo = np.flipud(np.interp(X, xi, np.flipud(SWd)))
    SWuo = np.flipud(np.interp(X, xi, np.flipud(SWu)))
    del X, xi

    # incident radiation on sunlit and shaded leaves Wm-2
    Q_sh = Clump*Kd*(SWdo + SWuo)  # normal to shaded leaves is all diffuse
    Q_sl = Kb*IbSky + Q_sh  # normal to sunlit leaves is direct and diffuse

    # absorbed components
    aLo = np.ones(len(Lo))*(1 - LeafAlbedo)
    aDiffo = aLo*Kd*(SWdo + SWuo)
    aDiro = aLo*Kb*IbSky

    # stand albedo
    alb = SWuo[-1] / (IbSky + IdSky + EPS)
    # print alb
    # soil absorption (Wm-2 (ground))
    q_soil = (1 - SoilAlbedo)*(SWdo[0] + SWbo[0])

    # correction to match absorption-based and flux-based albedo, relative error <3% may occur
    # in daytime conditions, at nighttime relative error can be larger but fluxes are near zero
    aa = (sum(aDiffo*Lo + aDiro*f_slo*Lo) + q_soil) / (IbSky + IdSky + EPS)
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

        plt.plot(f_slo, -Lcumo/Clump, 'r-', (1 - f_slo), -Lcumo/Clump, 'b-')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("sunlit & shaded fractions (-)")
        plt.legend(('f_{sl}, total = %.2f' % np.sum(f_slo*LAIz), 'f_{sh}, total = %.2f' % np.sum((1 - f_slo)*LAIz)), loc='best')

        # add input parameter values to fig
        plt.text(0.05, 0.75, r'$LAI$ = %1.1f m2 m-2' % (LAI))
        plt.text(0.05, 0.65, r'$ZEN$ = %1.3f rad' % (ZEN))
        plt.text(0.05, 0.55, r'$\alpha_l$ = %0.2f' % (LeafAlbedo))
        plt.text(0.05, 0.45, r'$\alpha_s$ = %0.2f' % (SoilAlbedo))

        plt.subplot(222)
        plt.plot(Q_sl, -Lcumo/Clump, 'ro-', Q_sh, -Lcumo/Clump, 'bo-')
        plt.plot(q_sl/(1-LeafAlbedo), -Lcumo/Clump, 'k-', q_sh/(1-LeafAlbedo), -Lcumo/Clump, 'k-')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("Incident radiation (Wm-2 (leaf))")
        plt.legend(('sunlit', 'shaded'), loc='best')

        plt.subplot(223)
        plt.plot(SWd, -Lcum/Clump, 'bo', SWdo, -Lcumo/Clump, 'b-', Ib, -Lcum/Clump, 'ro',
                 SWbo, -Lcumo/Clump, 'r-', SWu, -Lcum/Clump, 'go', SWuo, -Lcumo/Clump, 'g-')
        plt.legend(('SWd', 'Swdo', 'SWb', 'SWbo', 'SWu', 'SWuo'), loc='best')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("Incident SW (Wm-2 )")

        plt.subplot(224)
        plt.plot(q_sl, -Lcumo/Clump, 'ro-', q_sh, -Lcumo/Clump, 'bo-')
#        plt.plot((-SWdo[:-1]+SWdo[1:]-SWuo[1:]+SWuo[:-1])/(LAIz[:-1]+EPS),-Lcumo[1:]/Clump,'-k')
        plt.plot((1-np.exp(-Kd*Lo))*(SWdo + SWuo)/(LAIz+EPS),-Lcumo/Clump,'-k')
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

    # sunlit fraction in clumped foliage; clumping means elements shade each other
    f_sl = Clump*f_sl

    if PlotFigs:
        plt.figure(999)
        plt.subplot(221)
        plt.title("Source: radiation.canopy_sw_Spitters")

        plt.plot(f_sl, -Lcum/Clump, 'r-', (1 - f_sl), -Lcum/Clump, 'b-')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("sunlit & shaded fractions (-)")
        plt.legend(('f_{sl}, total = %.2f' % np.sum(f_sl*LAIz), 'f_{sh}, total = %.2f' % np.sum((1 - f_sl)*LAIz)), loc='best')

        # add input parameter values to fig
        plt.text(0.05, 0.75, r'$LAI$ = %1.1f m2 m-2' % (LAI))
        plt.text(0.05, 0.65, r'$ZEN$ = %1.3f rad' % (ZEN))
        plt.text(0.05, 0.55, r'$\alpha_l$ = %0.2f' % (LeafAlbedo))
        plt.text(0.05, 0.45, r'$\alpha_s$ = %0.2f' % (SoilAlbedo))

        plt.subplot(222)
        plt.plot(Q_sl, -Lcum/Clump, 'ro-', Q_sh, -Lcum/Clump, 'bo-')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("Incident radiation (Wm-2 (leaf))")
        plt.legend(('sunlit', 'shaded'), loc='best')

        plt.subplot(223)
        plt.plot(SWd, -Lcum/Clump, 'bo', SWb, -Lcum/Clump, 'ro')
        plt.legend(('SWd', 'SWb'), loc='best')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("Incident SW (Wm-2 )")

        plt.subplot(224)
        plt.plot(q_sl, -Lcum/Clump, 'ro-', q_sh, -Lcum/Clump, 'bo-')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("Absorbed radiation (Wm-2 (leaf))")
        plt.legend(('sunlit', 'shaded'), loc='best')

    return SWb, SWd, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_sl, alb

def compute_clouds_rad(doy, zen, Rg, H2O, Tair):
    """
    Estimates atmospheric transmissivity (tau_atm [-]), cloud cover fraction
    (f_cloud (0-1), [-]) and fraction of diffuse to total SW radiation (f_diff, [-])
    Args:
        doy - julian day
        zen - sun zenith angle (rad)
        Rg - total measured Global radiation above canopy (Wm-2)
        H2O - water vapor pressure (Pa)
    Returns:
        f_cloud - cloud cover fraction (0-1) [-]
        f_diff - fraction of diffuse to total radiation (0-1), [-]
        emi_sky - atmospheric emissivity (0-1) [-]

    ! values for Rg < 100 W/m2 linearly interpolated

    Cloudiness estimate is based on Song et al. 2009 JGR 114, 2009, Appendix A & C
    Clear-sky emissivity as in Niemelä et al. 2001 Atm. Res 58: 1-18.
    eq. 18 and cloud correction as Maykut & Church 1973.
    Reasonable fit against Hyytiälä data (tested 20.6.13)

    Samuli Launiainen, METLA, 2011-2013
    """

    # solar constant at top of atm.
    So = 1367.0
    # clear sky Global radiation at surface
    Qclear = np.maximum(0.0,
                        (So * (1.0 + 0.033 * np.cos(2.0 * np.pi * (np.minimum(doy, 365) - 10) / 365)) * np.cos(zen)))
    tau_atm = Rg / (Qclear + EPS)

    # cloud cover fraction
    f_cloud = np.maximum(0, 1.0 - tau_atm / 0.7)

    # calculate fraction of diffuse to total Global radiation: Song et al. 2009 JGR eq. A17.
    f_diff = np.ones(f_cloud.shape)
    f_diff[tau_atm > 0.7] = 0.2

    ind = np.where((tau_atm >= 0.3) & (tau_atm <= 0.7))
    f_diff[ind] = 1.0 - 2.0 * (tau_atm[ind] - 0.3)

    # clear-sky atmospheric emissivity
    ea = H2O / 100  # near-surface vapor pressure (hPa)
#    emi0 = np.where(ea >= 2.0, 0.72 + 0.009 * (ea - 2.0), 0.72 -0.076 * (ea - 2.0))
    emi0 = 1.24 * (ea/(Tair + 273.15))**(1./7.) # Song et al 2009

    # all-sky emissivity (cloud-corrections)
#    emi_sky = (1.0 + 0.22 * f_cloud**2.75) * emi0  # Maykut & Church (1973)
    emi_sky = (1 - 0.84 * f_cloud) * emi0 + 0.84 * f_cloud  # Song et al 2009 / (Unsworth & Monteith, 1975)

#    other emissivity formulas tested
#    emi_sky=(1 + 0.2*f_cloud)*emi0;
#    emi_sky=(1 + 0.3*f_cloud.^2.75)*emi0; % Maykut & Church (1973)
#    emi_sky=(1 + (1./emi0 -1)*0.87.*f_cloud^3.49).*emi0; % Niemelä et al. 2001 eq. 15 assuming Ts = Ta and surface emissivity = 1
#    emi_sky=(1-0.84*f_cloud)*emi0 + 0.84*f_cloud; % atmospheric emissivity (Unsworth & Monteith, 1975)

#    f_cloud[Rg < 100] = np.nan
#    f_diff[Rg < 100] = np.nan
#    emi_sky[Rg < 100] = np.nan

    f_cloud[Qclear < 10] = np.nan
    f_diff[Qclear < 10] = np.nan
    emi_sky[Qclear < 10] = np.nan

    df = pd.DataFrame({'f_cloud': f_cloud, 'f_diff': f_diff, 'emi_sky': emi_sky})
    df = df.interpolate()
    df = df.fillna(method='bfill')
    df = df.fillna(method='ffill')

    return df['f_cloud'].values, df['f_diff'].values, df['emi_sky'].values

def canopy_lw(LAIz, Clump, x, T, LWdn0, LWup0, leaf_emi=1.0, PlotFigs=False):
    """
    Estimates long-wave radiation budget and net isothermal LW radiation within the canopy
    Assumes canopy elements as black bodies (es=1.0) at local air temperature
    T(z). Neglects scattering etc.

    INPUT:
       LAIz - 1-sided LAI in each layer (m2m2(ground))
       CLUMP - clumping factor (-), [0...1]
       T: leaf temperature (degC)
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
    Kersti: modifications to indexing
    """

    N = len(LAIz)  # Node 0 = ground
    LWdn = np.zeros(N)
    LWup = np.zeros(N)
    LWleaf = np.zeros(N)
    LWnet = np.zeros(N)

    LayerLAI = Clump*LAIz  # plant-area m2m-2 in a layer addjusted for clumping
    cantop = max(np.where(LayerLAI>0)[0])  # node at canopy top

    # layerwise attenuation coeffcient
    Kd = kdiffuse(sum(Clump*LAIz), x)
    tau = np.exp(-Kd*LayerLAI)

    # LW down
    LWdn[cantop+1:N] = LWdn0  # downwelling LW entering layer i=cantop
    for k in range(cantop, -1, -1):
        LWdn[k]=tau[k]*LWdn[k+1] +(1 - tau[k])*(leaf_emi*STEFAN_BOLTZMANN*(T[k] + DEG_TO_KELVIN)**4)
    del k

    # LW up
    LWup[0] = LWup0  # upwelling LW entering layer i=0
    for k in range(1, cantop+2):
        LWup[k] = tau[k-1]*LWup[k-1] + (1 - tau[k-1])*(leaf_emi*STEFAN_BOLTZMANN*(T[k-1] + DEG_TO_KELVIN)**4)
    del k
    LWup[cantop+2:N] = LWup[cantop+1]

    # absorbed isothermal net radiation by the leaf (Wm-2(leaf))
    # Kd is mean projection of leaves
    LWleaf[0:cantop+1] = (1 - tau[0:cantop+1])*(
                        LWdn[1:cantop+2] + LWup[0:cantop+1] - 2*STEFAN_BOLTZMANN*leaf_emi*(T[0:cantop+1] + DEG_TO_KELVIN)**4)/(
                        LAIz[0:cantop+1] + EPS)
    LWnet[0:cantop+1] = (LWdn[1:cantop+2] - LWdn[0:cantop+1] + LWup[0:cantop+1] - LWup[1:cantop+2])/(LAIz[0:cantop+1]+EPS)

    gr = 2 * 4 * leaf_emi * STEFAN_BOLTZMANN * ( 1 - tau) * (T + DEG_TO_KELVIN) ** 3 / (LAIz + EPS) / SPECIFIC_HEAT_AIR

    if PlotFigs:
        Lcum = np.cumsum(np.flipud(LAIz))  # cumulative plant area index from canopy top
        Lcum = np.flipud(Lcum)
        plt.figure(99)
        plt.subplot(221)
        plt.title("radiation.canopy_lw", fontsize=8)

        plt.plot(LWdn, -Lcum, 'bo', label='LWdn')
        plt.plot(LWup, -Lcum, 'ro', label='LWup')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("LW (Wm-2 )")
        plt.legend()

        plt.subplot(222)
        plt.plot(LWnet, -Lcum, 'go',label='LWnet')
        plt.plot(LWleaf, -Lcum, 'ro',label='LWleaf')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("LW (Wm-2 )")
        plt.title('LWup0=%.1f, LWdn0=%.1f' % (LWup0, LWdn0))
        plt.legend()

        plt.subplot(223)
        plt.plot(LWdn, list(range(len(Lcum))), 'bo', label='LWdn')
        plt.plot(LWup, list(range(len(Lcum))), 'ro', label='LWup')
        plt.ylabel("N")
        plt.xlabel("LW (Wm-2 )")
        plt.legend()

        plt.subplot(224)
        plt.plot(LWleaf,list(range(len(Lcum))), 'ro',label='LWleaf')
#        plt.plot(LWnet,range(len(Lcum)), 'go',label='LWnet')
        plt.ylabel("N")
        plt.xlabel("LW (Wm-2 )")
        plt.title('LWup0=%.1f, LWdn0=%.1f' % (LWup0, LWdn0))
        plt.legend()
    return LWleaf, LWdn, LWup, gr

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
    Kersti

    CONVERTED TO PYTHON / LAST EDIT: 15.5.2017
    """
    # original and computational grid
    LAI = Clump*sum(LAIz)  # effective LAI, corrected for clumping (m2 m-2)
    Lo = Clump*LAIz  # effective layerwise LAI (or PAI) in original grid

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
    # index 0 at soil surface
    T = np.flipud(np.interp(Lcum, np.flipud(Lcumo), np.flipud(Tleaf)))  # for some reason x needs to be increasing..?
#    plt.figure()
#    plt.plot(T,Lcum,'og')
#    plt.plot(Tleaf,Lcumo,'or')

    # ---- optical parameters
    # back-scattering fraction, approximation, eq. (6-7)
    if x == 1:  # spherical leaf distrib.
        rd = 2./3.
    elif x == 0:  # vertical leafs
        rd = 0.5
    elif x > 100:  # horizontal leafs
        rd = 1.
    else:
        print("radiation.canopy_lw_ZhaoQualls: check leaf angle distr.")

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
    # initialize arrays: A=subdiag, B=diag, C=superdiag, D=rhs
    A = np.zeros(2*M+2)
    B = np.zeros(2*M+2)
    C = np.zeros(2*M+2)
    D = np.zeros(2*M+2)

    # subdiagonal
    A[1:2*M+1:2] = - (taud[1:M+1] + (1 - taud[1:M+1])*(1 - aL[1:M+1])*(1 - rd[1:M+1]))
    A[2:2*M+1:2] = 1 - rd[1:M+1]*rd[2:M+2]*(1 - aL[1:M+1])*(1 - taud[1:M+1])*(1 - aL[2:M+2])*(1 - taud[2:M+2])
    # diagonal
    B[0] = 1.0
    B[1:2*M+1:2] = - rd[0:M]*(taud[1:M+1] + (1 - taud[1:M+1])*(1 - aL[1:M+1])*(1 - rd[1:M+1]))*(
                    1 - aL[0:M])*(1 - taud[0:M])
    B[2:2*M+1:2] = - rd[2:M+2]*(taud[1:M+1] + (1 - taud[1:M+1])*(1 - aL[1:M+1])*(1 - rd[1:M+1]))*(
                    1 - aL[2:M+2])*(1 - taud[2:M+2])
    B[2*M+1] = 1.0
    # superdiagonal
    C[1:2*M+1:2] = 1 - rd[0:M]*rd[1:M+1]*(1 - aL[0:M])*(1 - taud[0:M])*(1 - aL[1:M+1])*(1 - taud[1:M+1])
    C[2:2*M+1:2] = - (taud[1:M+1] + (1 - taud[1:M+1])*(1 - aL[1:M+1])*(1 - rd[1:M+1]))

    # rhs
    LWsource = aL*STEFAN_BOLTZMANN*(T + DEG_TO_KELVIN)**4
    # lowermost row 0
    D[0] = LWup0
    # rows 1,3,5,...,M-3, M-1
    D[1:2*M+1:2] = (1 - rd[0:M]*rd[1:M+1]*(1 - aL[0:M])*(1 - taud[0:M])*(1 - aL[1:M+1])*(1 - taud[1:M+1]))*(
                    1 - taud[1:M+1]) *LWsource[1:M+1]
    # rows 2,4,6,..., M-2, M
    D[2:2*M+1:2] = (1 - rd[1:M+1]*rd[2:M+2]*(1 - aL[1:M+1])*(1 - taud[1:M+1])*(1 - aL[2:M+2])*(1 - taud[2:M+2]))*(
                    1 - taud[1:M+1])*LWsource[1:M+1]
    # uppermost row M+1
    D[2*M+1] = LWdn0

    # ---- solve a*LW = D
    if soil_emi < 1.0 and leaf_emi < 1.0:
        LW = tridiag(A,B,C,D)
    else:
        matrix = np.zeros([2*M+2, 2*M+2])
        row, col = np.diag_indices(matrix.shape[0])
        matrix[row, col] = B
        matrix[row[1:], col[:-1]] = A[1:]
        matrix[row[:-1], col[1:]] = C[:-1]
        LW = np.linalg.solve(matrix, D)
        del matrix, row, col

    # upward and downward hemispherical radiation (Wm-2 ground)
    LWu0 = LW[0:2*M+2:2]
    LWd0 = LW[1:2*M+2:2]
    del A, B, C, D, LW

    # ---- Compute multiple scattering, Zhao & Qualls, 2005. eq. (8 & 9)
    # downwelling diffuse after multiple scattering
    LWd = np.zeros(M+1)
    X = LWd0 / (1 - rd[0:M+1]*rd[1:M+2]*(1-aL[0:M+1])*(1 - taud[0:M+1])*(1 - aL[1:M+2])*(1 - taud[1:M+2]))
    Y = LWu0*rd[1:M+2]*(1 - aL[1:M+2])*(1 - taud[1:M+2]) / (1 - rd[0:M+1]*rd[1:M+2]*(1-aL[0:M+1])*(1 - taud[0:M+1])*(1 - aL[1:M+2])*(1 - taud[1:M+2]))
    LWd = X + Y

    # upwelling diffuse after multiple scattering
    LWu = np.zeros(M+1)
    X = LWu0 / (1 - rd[0:M+1]*rd[1:M+2]*(1-aL[0:M+1])*(1 - taud[0:M+1])*(1 - aL[1:M+2])*(1 - taud[1:M+2]))
    Y = LWd0*rd[0:M+1]*(1 - aL[0:M+1])*(1 - taud[0:M+1]) / (1 - rd[0:M+1]*rd[1:M+2]*(1-aL[0:M+1])*(1 - taud[0:M+1])*(1 - aL[1:M+2])*(1 - taud[1:M+2]))
    LWu = X + Y

    # --- NOW return values back to the original grid
    Lcum = Lcum[0:M+1]
    LWd = np.flipud(LWd)
    LWu = np.flipud(LWu)
    # plt.figure(99); plt.plot(LWd, -Lcum, 'r', LWu, -Lcum, 'g')

    X = Lcumo  # node 0 is canopy bottom
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
    LWleaf = np.zeros(len(taud))
    ic = np.where(LAIz > 0)[0]
    if len(ic) == 0:
        cantop = 0
    else:
        cantop = max(ic)
    LWleaf[0:cantop+1] = (1 - taud[0:cantop+1])*leaf_emi*(
                          LWdn[1:cantop+2] + LWup[0:cantop+1] - 2*STEFAN_BOLTZMANN*(Tleaf[0:cantop+1] + DEG_TO_KELVIN)**4)/(
                          LAIz[0:cantop+1] + EPS)

    gr = 2 * 4 * leaf_emi * STEFAN_BOLTZMANN * ( 1 - taud) * (Tleaf + DEG_TO_KELVIN) ** 3 / (LAIz + EPS) / SPECIFIC_HEAT_AIR

    if PlotFigs:
        plt.figure(99)
        plt.subplot(221)
        plt.title("radiation.canopy_lw_ZhaoQualls", fontsize=8)

        plt.plot(LWdn, -X/Clump, 'bo', label='LWdn')
        plt.plot(LWup, -X/Clump, 'ro', label='LWup')
#        plt.plot(LWd, -xi/Clump, 'go', label='LWdn_cg')
#        plt.plot(LWu, -xi/Clump, 'co', label='LWup_cg')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("LW (Wm-2 )")
        plt.legend()

        plt.subplot(222)
        plt.plot((-LWd[1:] + LWd[:-1] - LWu[:-1] + LWu[1:])/(L[1:-1]/Clump + EPS), -xi[1:]/Clump, 'go',label='LWnet')
        plt.plot(LWleaf, -X/Clump, 'ro', label='LWleaf')
        plt.ylabel("-Lcum eff.")
        plt.xlabel("LW (Wm-2 )")
        plt.title('LWup0=%.1f, LWdn0=%.1f' % (LWup0, LWdn0))
        plt.legend()

        plt.subplot(223)

        plt.plot(LWdn, list(range(len(X))), 'bo', label='LWdn')
        plt.plot(LWup, list(range(len(X))), 'ro', label='LWup')
        plt.ylabel("N")
        plt.xlabel("LW (Wm-2 )")
        plt.legend()

        plt.subplot(224)
        plt.plot(LWleaf, list(range(len(X))), 'ro', label='LWleaf')
#        plt.plot((LWdn[1:] - LWdn[:-1] + LWup[:-1] - LWup[1:])/(LAIz[:-1] + EPS), range(len(X)-1), 'go',label='LWnet')
        plt.ylabel("N")
        plt.xlabel("LW (Wm-2 )")
        plt.title('LWup0=%.1f, LWdn0=%.1f' % (LWup0, LWdn0))
        plt.legend()
    return LWleaf, LWdn, LWup, gr

def test_radiation_functions(LAI, Clump, ZEN, x=1.0, method="canopy_sw_ZhaoQualls", LAIz=None,
                             leaf_emi=0.98, soil_emi=0.98, LeafAlbedo=0.12, SoilAlbedo=0.1):
    """
    Runs test script for SW and LW radiation methods.
    INPUT:
        LAI: leaf area index (m2m-2)
        CLUMP: clumping index (-)
        ZEN: solar zenith angle (rad)
        x: spherical leaf-angle distr.
        method: Name of function to be tested "functionname"
    OUTPUT:
        none, prints figures
    AUTHOR:
        Samuli Launiainen, METLA 4/2014
    """

    # define setup for testing models

    #ZEN=35.0/180.0*np.pi
    IbSky = 100.0
    IdSky = 100.0

#    if LAIz == None:
#        LAIz = np.ones(102)*float(LAI) / 100.0
#        LAIz[0] = 0.
#        LAIz[-1] = 0.
    N = len(LAIz)

    # for LW calculations
    T = np.linspace(15, 17, N) # Tair is 15degC at ground and 17 at upper boundary
    Tatm = 17
    Tsurf = 16
    T = T * LAIz / (LAIz + EPS)
    LWdn0 = 0.85*STEFAN_BOLTZMANN*(Tatm + DEG_TO_KELVIN)**4
    LWup0 = 0.98*STEFAN_BOLTZMANN*(Tsurf + DEG_TO_KELVIN)**4
#    print LWdn0, LWup0
    if method == "canopy_sw_ZhaoQualls":
        print("------TestRun of radiation.canopy_sw_ZhaoQualls with given LAI and CLUMP -----------")
        SWb, SWd, SWu, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_sl, alb = canopy_sw_ZhaoQualls(LAIz, Clump, x, ZEN, IbSky, IdSky, LeafAlbedo, SoilAlbedo, PlotFigs="True")
        print(SWu[-1]/(SWb[-1]+SWd[-1]),alb)
#        print SWb,SWd,SWu,Q_sl,Q_sh,q_sl,q_sh,q_soil,f_sl,alb

    if method == "canopy_sw_Spitters":
        print("------TestRun of radiation.canopy_sw_Spitters with given LAI and predefined lad profile-----------")
        SWb, SWd, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_sl, alb = canopy_sw_Spitters(LAIz, Clump, x, ZEN, IbSky, IdSky, LeafAlbedo, SoilAlbedo, PlotFigs="True")
        # print SWb, SWd, Q_sl, Q_sh, q_sl, q_sh, q_soil, f_sl, alb

    if method == "canopy_lw":
        print("------TestRun of radiation.canopy_lw------------")
        LWnet, LWdn, LWup, gr = canopy_lw(LAIz, Clump, x, T, LWdn0, LWup0, leaf_emi=leaf_emi,PlotFigs=True)
        print(sum(LWnet*LAIz), LWdn[-1]-LWup[-1] - (LWdn[0]- LWup[0]))

    if method == "canopy_lw_ZhaoQualls":
        print("------TestRun of radiation.canopy_lw_ZhaoQualls with given LAI and CLUMP -----------")
        LWnet, LWdn, LWup, gr = canopy_lw_ZhaoQualls(LAIz, Clump, x, T, LWdn0, LWup0, leaf_emi=leaf_emi, soil_emi=soil_emi, PlotFigs=True)
        print(LWdn[-1],LWdn[-1]-LWup[-1])
