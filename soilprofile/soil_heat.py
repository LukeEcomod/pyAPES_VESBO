# -*- coding: utf-8 -*-
"""
soil.heat_flow contains functions for soil heat budget calculations

Created on Mon Apr 11 16:28:13 2016

@author: slauniai
Last edit: 3.7.2017 / Samuli
"""

import numpy as np
import matplotlib.pyplot as plt
from soil_water import wrc
from tools.utilities import tridiag as thomas

eps = np.finfo(float).eps  # machine epsilon

#: Constants

#: [m s\ :sup: `-2`\ ], standard gravity
GRAVITY = 9.81
#: [kg m\ :sup:`-3`\ ], ice density
ICE_DENSITY = 917.0
#: [K], normal temperature
NORMAL_TEMPERATURE = 273.15
#: [], latent heat of freezing
LATENT_HEAT_FREEZING = 333700.0
#: [J mol\ :sup: `-1`\ ], molar latent heat of vaporization at 20 \ :sup:`\circ'\C
L_MOLAR = 44100.0
#: [], freezing point of water
FREEZING_POINT_H2O = 0.0
#: [J mol\ :sup:`-1`\ ], universal gas constant
GAS_CONSTANT = 8.314
#: [kg mol\ :sup:`-1`\ ], molar mass of H\ :sub:`2`\ O
MOLAR_MASS_H2O = 18.015e-3
    
#: thermal condutivities [J m\ :sup:`-3`\ K \ :sup:`-1`\]
K_WATER = 0.57  # thermal cond of water W/(mK) = J/(smK)
K_ICE = 2.2  # thermal cond of ice J/(smK)
K_AIR = 0.025  # thermal conductivity of air J/(smK)
K_QUARTZ = 8.8  # thermal conductivity of Quarz  8.8 W/m/K
K_MINERAL = 2.9  # thermal conductivity ofother minerals W/m/K
K_ORG = 0.25  # thermal conductivity of organic matter W/m/K
K_SAND = 7.7 # thermal conductivity of sand W/m/K. Tian et al. 2007
K_SILT = 2.74 # thermal conductivity of silt W/m/K. Tian et al. 2007
K_CLAY = 1.93 # thermal conductivity of clay W/m/K. Tian et al. 2007

#: volumetric heat capacieties  [J m\ :sup:`-3`\ K \ :sup:`-1`\]
Cv_AIR = 1297.0  # air at 101kPa [Jm-3K-1]
Cv_WATER = 4.18e6  # water [Jm-3K-1]
Cv_ICE = 1.93e6  # ice [Jm-3K-1]
Cv_ORG = 2.50e6  # dry organic matter
Cv_MINERAL = 2.31e6  # soil minerals
Cv_GRANITE = 2.16e6  # granite
Cv_QUARTZ = 2.13e6  # quartz minerals

#: specific heats [J kg\ :sup:`-1`\ K \ :sup:`-1`\]
Cp_AIR = 1297.0
Cp_WATER = 4180.0
Cp_ICE = 2100.0
Cp_ORG = 1920.0
Cp_MINERAL = 870.0
Cp_GRANITE = 820.0
Cp_QUARTZ = 800.0
Cp_MINERALSOIL = 873.0  # Baldocchi's biometeorology notes

#: densities [kg m-3]
RHO_STONE = 2650.0  # solid mineral soils
RHO_ORG = 1300.0  # organic matter
RHO_WATER = 1000.0  # water
RHO_ICE = 920.0  # ice
RHO_AIR = 1.25  # air

"""
Volumetric heat capacities

"""

def volumetric_heat_capacity_mineral(bd):
    """
    volumetric heat capacity of dry mineral soil computed from bulk density
    Args:
        bd: bulk density [kgm-3]
    Returns:
        cv - volumetric heat capacity
    """
    return Cp_MINERALSOIL / bd


def volumetric_heat_capacity(poros, wliq=0.0, wice=0.0, cs=None, vOrg=0.0):
    """
    Computes volumetric heat capacity of soil
    IN:
        poros - porosity [m3 m-3]
        wliq - vol. liquid water content [m3 m-3]
        wice - vol. ice content [m3 m-3]
        cs - dry soil heat capacity [J m-3 K-1]
        orgfract - fraction of organic matter of solid volume, used to estimate cs if not given
    OUT:
        cv - volumetric heat capacity of soil [J m-3 K-1]
    """
    poros = np.array(poros, ndmin=1)
    
    # heat capacity of solid constituents in total solil volume [J m-3 K-1]
    if cs is None:
        vOrg = (1. - poros)*vOrg
        vMin = (1. - poros - vOrg)
        cs = vMin*Cv_MINERAL + vOrg*Cv_ORG

    wair = poros - wliq - wice

    cv = cs*(1 - poros) + Cv_WATER*wliq + Cv_ICE*wice + Cv_AIR*wair

    return cv

"""
Thermal conductivity of soils

"""


def thermal_conductivity_deVries(poros, wliq, wice=0.0, h=None, pF=None, T=15.0, ks=None, vOrg=0.0, vQuartz=0.0):
    """
    Thermal conductivity of soil
    IN:
        poros - porosity [m3m-3]
        wliq - liquid water content [m3m-3]
        wice - ice content [m3m-3]
        h - soil water tension [m]
        pF - vanGenuchten pF parameter dict
        T - temperature [degC]
        ks - soil parent material thermal conductivity
        vOrg - organic matter fraction (of total vol. of soil)
        vQuartz - quartz fraction (of total vol. of soil)
    """
    # K_WATER # thermal cond of water W/(mK) = J/(smK)
    # K_ICE  # thermal cond of ice J/(smK)
    # K_AIR  # thermal conductivity of air J/(smK)
    Dv = 24e-6  # molecular diffusivity [m2 s-1] of water vapor in air at 20degC
    P = 101300.0  # Pa

    # --------
    cf = np.divide(P, 287.05*(T + 273.15)) / 29e-3  # molar concentration of air [m3 mol-1]
    wliq = np.array(wliq, ndmin=1)
    wice = np.array(wice, ndmin=1)
    poros = np.array(poros, ndmin=1)

    N = len(wliq)
    wair = poros - wliq

    # vol. weighted solid conductivty [W m-1 K-1]
    if ks is None:
        ks = np.zeros(N) + (poros - vQuartz - vOrg)*K_MINERAL + vQuartz*K_QUARTZ + vOrg*K_ORG

    ga = np.ones(N)*0.0125
    ix = np.where(vOrg > 0.7*(1.0-poros))  # organic soil: organic matter content >40% of solid volume
    ga[ix] = 0.33  # organic soil

    # corrections for latent heat flux & vapor phase conductivity
    if h is not None:
        rh = soilRH(T, h)
    else:
        rh = 1.0

    es, ss = e_sat(T)
    wo = 0.10
    fw = np.divide(1.0, 1.0 + (wliq / wo)**-4.0)  # correction for latent heat flux

    # vapor phase apparent thermal conductivity inc. Stefan correction (W m-1 K-1)
    kg = K_AIR + L_MOLAR*ss*rh*fw*cf*Dv / (P - es)
    kf = kg + fw*(K_WATER - kg)  # fluid thermal cond.
    del fw, es, ss, cf

    # --- weight factors [-]------
    r = 2.0 / 3.0
    gg = r*(1.0 + ga*(kg / kf - 1.0))**(-1) + (1.0 - r)*(1.0 + (1.0 - 2*ga)*(kg / kf - 1.0))**(-1)  # gas phase                
    gw = r*(1.0 + ga*(K_WATER / kf - 1.0))**(-1) + (1.0 - r)*(1.0 + (1.0-2*ga)*(K_WATER / kf - 1.0))**(-1)  # liquid phase
    gs = r*(1.0 + ga*(ks / kf - 1.0))**(-1) + (1.0 - r)*(1.0 + (1.0-2*ga)*(ks / kf - 1.0))**(-1)  # solid phase

    L = (wliq*K_WATER*gw + kg*wair*gg + (1.0 - poros)*ks*gs + wice*K_ICE)\
        / (wliq*gw + wair*gg + (1.0 - poros)*gs + wice)  # W m-1 K-1                  

    return L


def thermal_conductivity(poros, wliq, wice=0.0, vOrg=0.0, vSand=0.0, vSilt=0.0, vClay=0.0, bedrockL=3.0):
    """
    Thermal conductivity of soil using Tian et al. 2016 simplified deVries-model. For organic layers
    (vOrg>0.9) uses o'Donnell et al. 2009. deVries-type model gives reasonable approximation also in
    organic layers when wliq <0.4, at ~0.9 ~25% underestimate compated to o'Donnell.
    
    IN:
        poros - porosity [m3m-3]
        wliq - liquid water content [m3m-3]
        wice - ice content [m3m-3]
        vOrg - organic matter volume fraction (of total vol. of soil)
        vClay - clay fraction (of total solid volume)
        vSilt - silt -"-
        vSand - sand -"-
        bedrockL - thermal conductivity of non-porous bedrock W m-1 K-1
    Refs:
        Tian et al. 2016: A simplified de Vries-based model to estimate thermal
        conductivity of unfrozen and frozen soil. E.J.Soil Sci 67, 564-572.
        o'Donnell et al. 2009. Thermal conductivity of organic soil horizons. Soil Sci. 174, 646-651.
    """
    poros = np.array(poros, ndmin=1)
    wliq = np.array(wliq, ndmin=1)
    #wliq = np.array(wliq, ndmin=1)
    wice = np.array(wice, ndmin=1)
    wair = poros - wliq - wice

    # component thermal conductivities W m-1 K-1
    ks = K_SAND**vSand * K_SILT**vSilt * K_CLAY**vClay * K_ORG**vOrg
    kw = K_WATER
    ki = K_ICE
    kg = K_AIR

    # shape factors
    g_sand = 0.182
    g_silt = 0.0534
    g_clay = 0.00775
    g_org = 0.33  # campbell & norman, 1998
    g_min = g_sand * vSand + g_silt * vSilt + g_clay * vClay + g_org * vOrg

    g_air = 0.333 * (1.0 - wair / poros)
    g_ice = 0.333 * (1.0 - wice / poros)
    # print('soil_heat(L226) g_air: {}, wair: {}, and poros: {}'.format(g_air, wair, poros))
     # --- weight factors [-]------
    r = 2.0 / 3.0
    # solid minerals
    f_min = r*(1.0 + g_min*(ks / kw - 1.0))**(-1) + (1.0 - r)*(1.0 + (1.0-2*g_min)*(ks / kw - 1.0))**(-1.)  
    # gas phase
    f_gas = r*(1.0 + g_air*(kg / kw - 1.0))**(-1) + (1.0 - r)*(1.0 + (1.0 - 2*g_air)*(kg / kw - 1.0))**(-1.)
    # liquid phase                
    f_ice = r*(1.0 + g_ice*(ki / kw - 1.0))**(-1) + (1.0 - r)*(1.0 + (1.0-2*g_ice)*(ki / kw - 1.0))**(-1.)
    
    #print f_min, f_gas
    # thermal conductivity W m-1 K-1
    L = (wliq*kw + f_gas*wair*kg + f_min*poros*ks + f_ice*wice*ki) \
        / (wliq + f_gas*wair + f_min*poros + f_ice*wice)                  

    # in fully organic layer, use o'Donnell et al. 2009
    L[vOrg >= 0.9] = 0.032 + 5e-1 * wliq[vOrg >= 0.9]

    # in bedrock
    L[poros == 0] = bedrockL

    return L
    

def thermal_conductivity_Campbell(poros, wliq, wice=0.0, vQuartz=0.0, vClay=0.0, vMineral=1.0):
    """
    Soil thermal conductivity using Campbell (1995); extended by Hansson et al. (2004) to frozen conditions
                    % c1-c5 are from Hydrus 1-D code
    CHECK THAT THIS IS CORRECT!!!
    """

    c1 = 0.57 + 1.73*vQuartz + 0.93*vMineral / (1.0 - 0.74*vQuartz)  # W/(mK)
    c2 = 2.8*(1.0-poros)  # W/(mK)
    c3 = 1 + 2.6*np.sqrt(np.maximum(0.005, vClay))  # [-]
    c4 = 0.03 + 0.7*(1.0 - poros)  # W/(mK)
    c5 = 4.0  # [-]

    #  Hansson et al, 2004 Kanagawa sandy loam        
    #         c1=0.55; %W/(mK)
    #         c2=0.80; %W/(mK)
    #         c3=3.07; %[-]
    #         c4=0.13; %W/(mK)
    #         c5=4; %[-])      

    f1 = 13.0  # [-]
    f2 = 1.06  # [-]
    F = 1.0 + f1*wice**f2  # [-]

    L = c1 + c2*(wliq + F*wice) - (c1 - c4)*np.exp(-(c3*(wliq + F*wice))**c5)  # W/(mK)

    return L


def thermal_conductivity_solids(vMineral=0.0, vQuartz=0.0, vOrg=0.0):
    """
    Soil solids thermal conductivity from mineral composition.
    Note: vMineral + vQuartz + vOrg = 1. floats or np.arrays
    Args:
        vMineral: mineral fraction of volume [-]
        vQuartz: quartz fraction -"-
        vOrg: organic fraction -"-
    Returns:
        thermal conductivity of solid material [Wm-1K-1]
    """
    vtot = vMineral + vQuartz + vOrg
    ks = (K_MINERAL*vMineral + K_QUARTZ*vQuartz + K_ORG*vOrg) / vtot
    
    return ks


def thermal_conductivity_organic(wliq):
    """
    Organic matter thermal conductivity based on O'Donnel et al. (2009, Soil Sci. 174).
    Args:
        wliq: vol. water content [-]
    Returns:
        thermal conductivity [Wm-1K-1]
    """
    L = 3.0e-2 + 5.0e-1*wliq
    L = np.maximum(L, 4e-2)

    return L


def thermal_conductivity_simple(poros, wliq, wice, ks=None, soilComp=None):
    """
    Simple estimate for thermal conductivity in soil as volume-weighted average of constituents.
    Args:
        poros (float/np.array): porosity [-]
        wliq (float/np.array): vol. water content [-]
        wice (float/np.array: vol. ice content [-]
        ks (float/np.array): dry soil (solids only) thermal conductivity [Wm-1K-1]
        soil_comp (dict):
            * 'vMineral', fraction of total volume
            * 'vQuartz', -"-
            * 'vOrg', -"-
    Returns:
        thermal conductivity (float/np.array) [Wm-1K-1]
    """
    wair = poros - wliq - wice
       
    if ks is not None:  # dry soil thermal cond. given as input
        L = ks + K_WATER*wliq + K_ICE*wice + K_AIR*wair
    elif soilComp is not None:
        ks = K_MINERAL*soilComp['vMineral'] + K_QUARTZ*soilComp['vQuartz'] + K_ORG*soilComp['vOrg']
        L = ks + K_WATER*wliq + K_ICE*wice + K_AIR*wair

    return L

"""
Soil heat transfer in 1D

"""
def heatflow_1D_new(t_final, grid, poros, T0, Wtot, ubc, lbc, spara, S=0.0, steps=10):
    """
     Solves soil heat flow in 1-D using implicit, backward finite difference solution
     of heat equation (Hansson et al., 2004; Saito et al., 2006):

    .. math::
        \\frac{\\partial C_p T}{\\partial t} =
        \\frac{\\partial}{\\partial z} \\left[\\lambda(\\theta)\\frac{\\partial T}{\\partial z}\\right] +
        C_w\\frac{\\partial q T}{\\partial z}

    where :math:`C_p` is volumetric heat capacity of soil, :math:`C_w` is volumetric heat
    capacity of liquid water, :math:`\\lambda(z)` is heat conductivity in soil, and
    :math:`q(z)` is the liquid water flow.

    Backward finite difference solution adapted from the on of Van Dam and Feddes (2000) for
    solving Richars equation.

    Reference:
        Hansson et al. (2004). Vadose Zone Journal 3:693-704
        Saito et al. (2006). Vadose Zone Journal 5:784-800
        van Dam and Feddes (2000). J. Hydrology xxxx

    Args:
        t_final: solution timestep [s]
        z: grid,<0, monotonically decreasing [m]
        poros: porosity
        T0: initial temperature profile [degC]
        Wtot: volumetric water content [m3m-3]
        ubc: upper bc: {'type': (give str 'flux','temperature'), 'value': ubc value}.
            Downward flux <0 (i.e. surface heating / bottom cooling)
        lbc: lower bc, formulate as ubc
        spara: soil type parameter dict with necessary fields for
            * 'cs':
                dry soil vol. heat capacity [Jm\ :sup:`-3` K\ :sup:`-1`\ ]
            * 'fp':
                freezing curve shape parameter
                    * 2...4 for clay soils
                    * 0.5...1.5 for sandy soils
                    * < 0.5 for organic soils
            * 'vOrg': organic matter volume fraction (of solid volume)
            * 'vSand': sand -"-
            * 'vSilt': silt -"-
            * 'vClay': clay -"-
        S:
            local heat sink/source array [Wm-3 =  Js-1m3], <0 for sink
        steps:
            initial subtimesteps used to proceed to 't_final'

    Returns:
        T: new temperature profile [m]
        Wliq: new liquid water content [m3m-3]
        Wice: new ice ontent [m3m-3]
        Fheat: heat flux array [Wm-2]
        R: thermal conductivity [Wm-1K-1]
    CODE:
        Samuli Launiainen, Luke 19.4.2016. Converted from Matlab (APES SoilProfile.Soil_HeatFlow)
        Kersti, 26.7.2018 edits to avoid computation getting stuck during freezing conditions,
            still problems with convergence around zero temperatures.
            Check heat balance error to see problematic periods
    NOTE:
        (19.4.2016): Tested ok; solved 'coupled' water&heatflow solve heatflow and waterflow
        sequentially for small timesteps
    TODO:
        1. think how heat advection with infiltration should be accounted for
        2. gradient boundary should be implemented as flux-based
    """

    Conv_crit1 = 1.0e-3  # degC
    Conv_crit2 = 1.0e-5  # ice content m3/m3

    LH = 0.0    # heat sink/source to 1st node due infiltration or evaporation
                # LH=-cw*infil_rate *T_infil/ dz[0] [Wm-3]

    # -------------- constants ---------------

    rhoi = RHO_ICE  # density of ice
    Lf = LATENT_HEAT_FREEZING  # latent heat of freezing, J/kg; % Lv is latent heat of vaportization
    Tfr = FREEZING_POINT_H2O  # freezing point of water (degC)

    # ------------------- computation grid -----------------------

    # grid
    z = grid['z']
    dz = grid['dz']
    dzu = grid['dzu']
    dzl = grid['dzl']

    N = len(z)

    # -----get parameters from dict
    fp = spara['fp']       # freezing-curve parameter, 2...4 for clay and 0.5-1.5 for sandy soils
                            # (Femma-code/Harri Koivusalo)

    cs = spara['cs']        # dry soil vol. heat capacity [J m-3 K-1]
    vorg = spara['vOrg']    # organic matter vol. fraction [-]
    vsand = spara['vSand']  # sand vol. fraction [-]
    vsilt = spara['vSilt']  # silt vol. fraction [-]
    vclay = spara['vClay']  # clay vol. fraction [-]
    bedrockL = spara['bedrockL']

    S = np.ones(N)*S  # sink-source array [Wm-3 =  Js-1m3], < 0 is heat sink

    """ -------- find solution at t_final --------- """

    # ---initial conditions
    T = T0
    # ---initiate heat flux array [W/m3]
    Fheat = np.zeros(N+1)

    # Liquid and ice content, and dWliq/dTs
    Wliq, Wice, gamma = frozen_water(T, Wtot, fp=fp, To=Tfr)

    # running time [s]
    t = 0.0
    # initial and computational time step [s]
    dto = t_final / steps
    dt = dto  # adjusted during solution

    while t < t_final:  # loop until solution timestep
        # these will stay cduring iteration over time step "dt"
        T_old = T
        Wice_old = Wice

        # old bulk soil heat capacity [Jm-3K-1]
        CP_old = volumetric_heat_capacity(poros, Wliq, wice=Wice, cs=cs)

        # Fraction of water and ice assumed constant during dt, so these only computed here
        # Thermal conductivity - this remains constant during iteration
        R = thermal_conductivity(poros, Wliq, wice=Wice, vOrg=vorg, vSand=vsand, vSilt=vsilt,
                                 vClay=vclay, bedrockL=bedrockL)
        R = spatial_average2(R, method='arithmetic')

        # changes during iteration
        T_iter = T.copy()
        Wice_iter = Wice.copy()

        err1 = 999.0
        err2 = 999.0
        iterNo = 0

        # start iterative solution of heat equation

        while err1 > Conv_crit1 or err2 > Conv_crit2:
            # print 'err1=' +str(err1) +'   err2=' + str(err2)
            iterNo += 1

            # bulk soil heat capacity [Jm-3K-1]
            CP = volumetric_heat_capacity(poros, Wliq, wice=Wice, cs=cs)
            # heat capacity due to freezing/thawing [Jm-3K-1]
            A = rhoi*Lf*gamma

            # --- set up tridiagonal matrix
            a = np.zeros(N)
            b = np.zeros(N)
            g = np.zeros(N)
            f = np.zeros(N)

            # ---- intermediate nodes
            b[1:-1] = CP[1:-1] + A[1:-1] + dt / dz[1:-1] * (R[1:N-1] / dzu[1:-1] + R[2:N] / dzl[1:-1])
            a[1:-1] = - dt / (dz[1:-1]*dzu[1:-1]) * R[1:N-1]
            g[1:-1] = - dt / (dz[1:-1]*dzl[1:-1]) * R[2:N]
            f[1:-1] = CP_old[1:-1] * T_old[1:-1] + A[1:-1] * T_iter[1:-1] \
                    + Lf*rhoi*(Wice_iter[1:-1] - Wice_old[1:-1]) - S[1:-1] * dt

            # ---------- top node (n=0)
            # LH is heat input by infiltration. loss by evaporation not currently implemented

            if ubc['type'] == 'flux':  # or ubc['type'] is 'grad':
                F_sur = ubc['value']
                b[0] = CP[0] + A[0] + dt / (dz[0] * dzl[0]) * R[1]
                a[0] = 0.0
                g[0] = -dt / (dz[0] * dzl[0]) * R[1]
                f[0] = CP_old[0]*T_old[0] + A[0]*T_iter[0] + Lf*rhoi*(Wice_iter[0] - Wice_old[0])\
                        - dt / dz[0] * F_sur - dt*LH - dt*S[0]

            if ubc['type'] == 'temperature':   # fixed T at imaginary node at surface
                T_sur = ubc['value']
                b[0] = CP[0] + A[0] + dt / dz[0]* (R[0] / dzu[0] + R[1] / dzl[0])
                a[0] = 0.0
                g[0] = -dt / (dz[0]*dzl[0]) * R[1]
                f[0] = CP_old[0]*T_old[0] + A[0]*T_iter[0] + Lf*rhoi*(Wice_iter[0] - Wice_old[0])\
                        + dt / (dz[0]*dzu[0])*R[0]*T_sur - dt*LH - dt*S[0]

            # ------ bottom node (n=N)
            if lbc['type'] == 'flux':  # or lbc['type'] is 'grad':
                F_bot = lbc['value']
                b[-1] = CP[-1] + A[-1] + dt / (dz[-1]*dzu[-1]) * R[N-1]
                a[-1] = -dt / (dz[-1]*dzu[-1]) * R[N-1]
                g[-1] = 0.0
                f[-1] = CP_old[-1]*T_old[-1] + A[-1]*T_iter[-1] + Lf*rhoi*(Wice_iter[-1] - Wice_old[-1])\
                        - dt / dz[-1]*F_bot - dt*S[-1]
    
            if lbc['type'] == 'temperature':  # fixed temperature, Tbot "at node N+1"
                T_bot = lbc['value']
                b[-1] = CP[-1] + A[-1] + dt / dz[-1] * (R[N-1] / dzu[-1] + R[N] / dzl[-1])
                a[-1] = -dt / (dz[-1]*dzu[-1]) * R[N-1]
                g[-1] = 0.0
                f[-1] = CP_old[-1]*T_old[-1] + A[-1]*T_iter[-1] + Lf*rhoi*(Wice_iter[-1] - Wice_old[-1])\
                        + dt/(dz[-1]*dzl[-1]) * R[N]*T_bot - dt*S[-1]

            # --------- update iteration values
            T_iterold = T_iter.copy()
            Wice_iterold = Wice_iter.copy()

            # --- call tridiagonal solver
            T_iter = thomas(a, b, g, f)

            Wliq_iter, Wice_iter, gamma = frozen_water(T_iter, Wtot, fp=fp)

            # if problems reaching convergence devide time step and retry
            if iterNo == 20:
                dt = dt / 3.0
                if dt > 10:
                    dt = max(dt, 30)
                    iterNo = 0
                    T_iter = T_old.copy()
                    Wice_iter = Wice_old.copy()
                    #print 'soil_heat.heatflow_1D: More than 20 iterations, new dt = ' + str(dt) + 'Terr1 = ' + str(err1)+ 'Terr2 = ' + str(err2)
                    continue
                else:  # convergence not reached with dt=30s, break
                    #print 'soil_heat.heatflow_1D: Solution not converging' + 'Terr1 = ' + str(err1)+ 'Terr2 = ' + str(err2) + 'Tair = ' + str(T_sur)
                    break
            # check solution, if problem continues break
            elif any(np.isnan(T_iter)):
                dt = dt / 3.0
                if dt > 10:
                    dt = max(dt, 30)
                    iterNo = 0
                    T_iter = T_old.copy()
                    Wice_iter = Wice_old.copy()
                    #print 'soil_heat.heatflow_1D: Solution blowing up, new dt = ' + str(dt)
                    continue
                else:  # convergence not reached with dt=30s, break
                    #print 'soil_heat.heatflow_1D: Problem with solution, blow up'
                    break

            err1 = np.max(abs(T_iter - T_iterold))
            err2 = np.max(abs(Wice_iter - Wice_iterold))
#            print 'Terr1 = ' + str(err1) + ' Terr2 = ' + str(err2)

        # ------- ending iteration loop -------
        # print 'Pass: t = ' +str(t) + ' dt = ' +str(dt) + ' iterNo = ' + str(iterNo)

        # update state
        T = T_iter.copy()
        # Liquid and ice content, and dWice/dTs
        Wliq, Wice, gamma = frozen_water(T, Wtot, fp=fp, To=Tfr)
        # Heat flux [W m-2]
        Fheat[1:-1] += -R[1:-1]*(T[1:] - T[:-1])/dzl[:-1] * dt / t_final
        Fheat[0] += -R[0]*(T[0] - T_sur)/dzu[0] * dt / t_final
        Fheat[-1] += -R[-1]*(T_bot - T[-1])/dzl[-1] * dt / t_final

        # ----------- solution time and new timestep ------------

        t += dt
#        print 'soil_heat.heatflow_1D: t = ' + str(t) + ' dt = ' + str(dt) + ' iterNo = ' + str(iterNo) + 'Terr1 = ' + str(err1)

        dt_old = dt  # save temporarily

        # select new time step based on convergence
        if iterNo <= 3:
            dt = dt * 2
        elif iterNo >= 6:
            dt = dt / 2

        # limit to minimum of 30s
        dt = max(dt, 30)

        # save dto for output to be used in next run of heatflow1D()
        if dt_old == t_final or t_final > t:
            dto = min(dt, t_final)

        # limit by time left to solve
        dt = min(dt, t_final-t)

    # ----- now at t = t_final, end while loop

    def heat_content(x):
        wliq, wice, _ = frozen_water(x, Wtot, fp=fp, To=Tfr)
        Cv = volumetric_heat_capacity(poros, wliq=wliq, wice=Wice, cs=cs)
        return Cv * x - Lf * rhoi * wice

    def heat_balance(x):
        return heat_content(T0) + t_final * (Fheat[:-1] - Fheat[1:]) / dz - heat_content(x)

    heat_be = sum(dz * heat_balance(T))
    return T, Wliq, Wice, Fheat, R, dto, heat_be

def heatflow_1D_new2(t_final, grid, poros, T0, Wtot, ubc, lbc, spara, S=0.0, steps=10):
    """
     Solves soil heat flow in 1-D using implicit, backward finite difference solution
     of heat equation (Hansson et al., 2004; Saito et al., 2006):

    .. math::
        \\frac{\\partial C_p T}{\\partial t} =
        \\frac{\\partial}{\\partial z} \\left[\\lambda(\\theta)\\frac{\\partial T}{\\partial z}\\right] +
        C_w\\frac{\\partial q T}{\\partial z}

    where :math:`C_p` is volumetric heat capacity of soil, :math:`C_w` is volumetric heat
    capacity of liquid water, :math:`\\lambda(z)` is heat conductivity in soil, and
    :math:`q(z)` is the liquid water flow.

    Backward finite difference solution adapted from the on of Van Dam and Feddes (2000) for
    solving Richars equation.

    Reference:
        Hansson et al. (2004). Vadose Zone Journal 3:693-704
        Saito et al. (2006). Vadose Zone Journal 5:784-800
        van Dam and Feddes (2000). J. Hydrology xxxx

    Args:
        t_final: solution timestep [s]
        z: grid,<0, monotonically decreasing [m]
        poros: porosity
        T0: initial temperature profile [degC]
        Wtot: volumetric water content [m3m-3]
        ubc: upper bc: {'type': (give str 'flux','temperature'), 'value': ubc value}.
            Downward flux <0 (i.e. surface heating / bottom cooling)
        lbc: lower bc, formulate as ubc
        spara: soil type parameter dict with necessary fields for
            * 'cs':
                dry soil vol. heat capacity [Jm\ :sup:`-3` K\ :sup:`-1`\ ]
            * 'fp':
                freezing curve shape parameter
                    * 2...4 for clay soils
                    * 0.5...1.5 for sandy soils
                    * < 0.5 for organic soils
            * 'vOrg': organic matter volume fraction (of solid volume)
            * 'vSand': sand -"-
            * 'vSilt': silt -"-
            * 'vClay': clay -"-
        S:
            local heat sink/source array [Wm-3 =  Js-1m3], <0 for sink
        steps:
            initial subtimesteps used to proceed to 't_final'

    Returns:
        T: new temperature profile [m]
        Wliq: new liquid water content [m3m-3]
        Wice: new ice ontent [m3m-3]
        Fheat: heat flux array [Wm-2]
        R: thermal conductivity [Wm-1K-1]
    CODE:
        Samuli Launiainen, Luke 19.4.2016. Converted from Matlab (APES SoilProfile.Soil_HeatFlow)
        Kersti, 26.7.2018, computation following method in FEMMA but correction of T based 
            on heat balance is not working, bisection method may work?
            Now computes temperature profile with constant Wice (results in heat balance 
            error during freezing/thawing)
    NOTE:
        (19.4.2016): Tested ok; solved 'coupled' water&heatflow solve heatflow and waterflow
        sequentially for small timesteps
    TODO:
        1. think how heat advection with infiltration should be accounted for
        2. gradient boundary should be implemented as flux-based
    """

    Conv_crit1 = 1.0e-3  # degC

    LH = 0.0    # heat sink/source to 1st node due infiltration or evaporation
                # LH=-cw*infil_rate *T_infil/ dz[0] [Wm-3]

    # -------------- constants ---------------

    rhoi = RHO_ICE  # density of ice
    Lf = LATENT_HEAT_FREEZING  # latent heat of freezing, J/kg; % Lv is latent heat of vaportization
    Tfr = FREEZING_POINT_H2O  # freezing point of water (degC)

    # ------------------- computation grid -----------------------

    # grid
    z = grid['z']
    dz = grid['dz']
    dzu = grid['dzu']
    dzl = grid['dzl']

    N = len(z)

    # -----get parameters from dict
    fp = spara['fp']       # freezing-curve parameter, 2...4 for clay and 0.5-1.5 for sandy soils
                            # (Femma-code/Harri Koivusalo)

    cs = spara['cs']        # dry soil vol. heat capacity [J m-3 K-1]
    vorg = spara['vOrg']    # organic matter vol. fraction [-]
    vsand = spara['vSand']  # sand vol. fraction [-]
    vsilt = spara['vSilt']  # silt vol. fraction [-]
    vclay = spara['vClay']  # clay vol. fraction [-]
    bedrockL = spara['bedrockL']

    S = np.ones(N)*S  # sink-source array [Wm-3 =  Js-1m3], < 0 is heat sink

    """ -------- find solution at t_final --------- """

    # ---initial conditions
    T = T0
    # ---initiate heat flow array
    Fheat = np.zeros(N+1)

    # Liquid and ice content, and dWliq/dTs
    Wliq, Wice, gamma = frozen_water(T, Wtot, fp=fp, To=Tfr)

    # running time [s]
    t = 0.0
    # initial and computational time step [s]
    dto = t_final / steps
    dt = dto  # adjusted during solution

    while t < t_final:  # loop until solution timestep
        # these will stay cduring iteration over time step "dt"
        T_old = T

        # Fraction of water and ice assumed constant during dt, so these only computed here
        # Thermal conductivity
        R = thermal_conductivity(poros, Wliq, wice=Wice, vOrg=vorg, vSand=vsand, vSilt=vsilt,
                                 vClay=vclay, bedrockL=bedrockL)
        R = spatial_average2(R, method='arithmetic')
        # bulk soil heat capacity [Jm-3K-1]
        CP = volumetric_heat_capacity(poros, Wliq, wice=Wice, cs=cs)
        # heat capacity of freezing/thawing [Jm-3K-1]
        A = rhoi*Lf*gamma

        # changes during iteration
        T_iter = T.copy()

        err1 = 999.0
        iterNo = 0

        # start iterative solution of heat equation

        while err1 > Conv_crit1:
            # print 'err1=' +str(err1) +'   err2=' + str(err2)
            iterNo += 1

            # --- set up tridiagonal matrix
            a = np.zeros(N)
            b = np.zeros(N)
            g = np.zeros(N)
            f = np.zeros(N)

            # ---- intermediate nodes
            b[1:-1] = CP[1:-1] + A[1:-1] + dt / dz[1:-1] * (R[1:N-1] / dzu[1:-1] + R[2:N] / dzl[1:-1])
            a[1:-1] = - dt / (dz[1:-1]*dzu[1:-1]) * R[1:N-1]
            g[1:-1] = - dt / (dz[1:-1]*dzl[1:-1]) * R[2:N]
            f[1:-1] = CP[1:-1] * T_old[1:-1] + A[1:-1] * T_old[1:-1] - S[1:-1] * dt

            # ---------- top node (n=0)
            # LH is heat input by infiltration. loss by evaporation not currently implemented

            if ubc['type'] == 'flux':  # or ubc['type'] is 'grad':
                F_sur = ubc['value']
                b[0] = CP[0] + A[0] + dt / (dz[0] * dzl[0]) * R[1]
                a[0] = 0.0
                g[0] = -dt / (dz[0] * dzl[0]) * R[1]
                f[0] = CP[0]*T_old[0] + A[0]*T_old[0] - dt / dz[0] * F_sur - dt*LH - dt*S[0]

            if ubc['type'] == 'temperature':   # fixed T at imaginary node at surface
                T_sur = ubc['value']
                b[0] = CP[0] + A[0] + dt / dz[0]* (R[0] / dzu[0] + R[1] / dzl[0])
                a[0] = 0.0
                g[0] = -dt / (dz[0]*dzl[0]) * R[1]
                f[0] = CP[0]*T_old[0] + A[0]*T_old[0] + dt / (dz[0]*dzu[0])*R[0]*T_sur - dt*LH - dt*S[0]

            # ------ bottom node (n=N)
            if lbc['type'] == 'flux':  # or lbc['type'] is 'grad':
                F_bot = lbc['value']
                b[-1] = CP[-1] + A[-1] + dt / (dz[-1]*dzu[-1]) * R[N-1]
                a[-1] = -dt / (dz[-1]*dzu[-1]) * R[N-1]
                g[-1] = 0.0
                f[-1] = CP[-1]*T_old[-1] + A[-1]*T_old[-1] - dt / dz[-1]*F_bot - dt*S[-1]
    
            if lbc['type'] == 'temperature':  # fixed temperature, Tbot "at node N+1"
                T_bot = lbc['value']
                b[-1] = CP[-1] + A[-1] + dt / dz[-1] * (R[N-1] / dzu[-1] + R[N] / dzl[-1])
                a[-1] = -dt / (dz[-1]*dzu[-1]) * R[N-1]
                g[-1] = 0.0
                f[-1] = CP[-1]*T_old[-1] + A[-1]*T_old[-1] + dt/(dz[-1]*dzl[-1]) * R[N]*T_bot - dt*S[-1]

            # --------- update iteration values
            T_iterold = T_iter.copy()

            # --- call tridiagonal solver
            T_iter = thomas(a, b, g, f)

            # if problems reaching convergence devide time step and retry
            if iterNo == 20:
                dt = dt / 3.0
                if dt > 10:
                    dt = max(dt, 30)
                    iterNo = 0
                    print 'soil_heat.heatflow_1D: More than 20 iterations, new dt = ' + str(dt) + 'Terr1 = ' + str(err1)
                    continue
                else:  # convergence not reached with dt=30s, break
                    print 'soil_heat.heatflow_1D: Solution not converging' + 'Terr1 = ' + str(err1)
                    break
            # check solution, if problem continues break
            elif any(np.isnan(T_iter)):
                dt = dt / 3.0
                if dt > 10:
                    dt = max(dt, 30)
                    iterNo = 0
                    T_iter = T_old.copy()
                    print 'soil_heat.heatflow_1D: Solution blowing up, new dt = ' + str(dt)
                    continue
                else:  # convergence not reached with dt=30s, break
                    print 'soil_heat.heatflow_1D: Problem with solution, blow up'
                    break

            err1 = np.max(abs(T_iter - T_iterold))
#            print 'Terr1 = ' + str(err1) + ' Terr2 = ' + str(err2)

        # ------- ending iteration loop -------
        # print 'Pass: t = ' +str(t) + ' dt = ' +str(dt) + ' iterNo = ' + str(iterNo)

        # update state
        T = T_iter.copy()
        # Liquid and ice content, and dWice/dTs
        Wliq, Wice, gamma = frozen_water(T, Wtot, fp=fp, To=Tfr)
        # Heat flux [W m-2]
        Fheat[1:-1] += -R[1:-1]*(T[1:] - T[:-1])/dzl[:-1] * dt / t_final
        Fheat[0] += -R[0]*(T[0] - T_sur)/dzu[0] * dt / t_final
        Fheat[-1] += -R[-1]*(T_bot - T[-1])/dzl[-1] * dt / t_final

        # ----------- solution time and new timestep ------------

        t += dt
#        print 'soil_heat.heatflow_1D: t = ' + str(t) + ' dt = ' + str(dt) + ' iterNo = ' + str(iterNo) + 'Terr1 = ' + str(err1)

        dt_old = dt  # save temporarily

        # select new time step based on convergence
        if iterNo <= 3:
            dt = dt * 2
        elif iterNo >= 6:
            dt = dt / 2

        # limit to minimum of 30s
        dt = max(dt, 30)

        # save dto for output to be used in next run of heatflow1D()
        if dt_old == t_final or t_final > t:
            dto = min(dt, t_final)

        # limit by time left to solve
        dt = min(dt, t_final-t)

    # ----- now at t = t_final, end while loop

    def heat_content(x):
        wliq, wice, _ = frozen_water(x, Wtot, fp=fp, To=Tfr)
        Cv = volumetric_heat_capacity(poros, wliq=wliq, wice=Wice, cs=cs)
        return Cv * x - Lf * rhoi * wice

    def dheat_content(x):
        dx = 1e-3
        return -(heat_content(x+dx) - heat_content(x-dx))/(2*dx)
    
    def heat_balance(x):
        return heat_content(T0) + t_final * (Fheat[:-1] - Fheat[1:]) / dz - heat_content(x)

#    This part didn't work out..
#    T_new = newtons_method(heat_balance, dheat_content, T, Conv_crit2)
#
#    T = T_new.copy()
#    Wliq, Wice, _ = frozen_water(T, Wtot, fp=fp, To=Tfr)

    heat_be = sum(dz * heat_balance(T))
    return T, Wliq, Wice, Fheat, R, dto, heat_be

def newtons_method(f, df, x0, e):
    alfa=1.0
    i = 0
    while  abs(sum(f(x0))) > e:# and i < 20:
        x0 = x0 - alfa*f(x0)/df(x0)
        i += 1
        if i % 20 == 0:
            alfa *= 0.5
        if i>100:
            print 'Newtons method reached 100 iterations!'
            break
    print 'Newtons method iterations: ' + str(i)
    return x0

def bisection_method(f, x1, x2, e):
    x0 = (x1 + x2) / 2.0
    delta = np.max((x2 - x1) / 2.0)
    i = 0
    while delta > e and i < 50:
        ix = np.where(f(x1)*f(x0) < 0)
        ixx = np.where(f(x1)*f(x0) > 0)
        x2[ix] = x0
        x1[ixx] = x0
        x0 = (x1 + x2) / 2.0
        delta = np.max((x2 - x1) / 2.0)
        i += 1
    print 'Bisection method iterations: ' +  str(i)
    return x0

def heatflow_1D(t_final, grid, poros, T0, Wliq0, Wice0, ubc, lbc, spara, S=0.0, steps=10):
    """
     Solves soil heat flow in 1-D using implicit, backward finite difference solution
     of heat equation (Hansson et al., 2004; Saito et al., 2006):

    .. math::
        \\frac{\\partial C_p T}{\\partial t} =
        \\frac{\\partial}{\\partial z} \\left[\\lambda(\\theta)\\frac{\\partial T}{\\partial z}\\right] +
        C_w\\frac{\\partial q T}{\\partial z}

    where :math:`C_p` is volumetric heat capacity of soil, :math:`C_w` is volumetric heat
    capacity of liquid water, :math:`\\lambda(z)` is heat conductivity in soil, and
    :math:`q(z)` is the liquid water flow.

    Backward finite difference solution adapted from the on of Van Dam and Feddes (2000) for
    solving Richars equation.

    Reference:
        Hansson et al. (2004). Vadose Zone Journal 3:693-704
        Saito et al. (2006). Vadose Zone Journal 5:784-800
        van Dam and Feddes (2000). J. Hydrology xxxx

    Args:
        t_final: solution timestep [s]
        z: grid,<0, monotonically decreasing [m]
        poros: porosity
        T0: initial temperature profile [degC]
        Wliq0: liquid water content [m3m-3]
        Wice0: ice content [m3m-3]
        ubc: upper bc: {'type': (give str 'flux','temperature'), 'value': ubc value}.
            Downward flux <0 (i.e. surface heating / bottom cooling)
        lbc: lower bc, formulate as ubc
        spara: soil type parameter dict with necessary fields for
            * 'cs':
                dry soil vol. heat capacity [Jm\ :sup:`-3` K\ :sup:`-1`\ ]
            * 'fp':
                freezing curve shape parameter
                    * 2...4 for clay soils
                    * 0.5...1.5 for sandy soils
                    * < 0.5 for organic soils
            * 'vOrg': organic matter volume fraction (of solid volume)
            * 'vSand': sand -"-
            * 'vSilt': silt -"-
            * 'vClay': clay -"-
        S:
            local heat sink/source array [Wm-3 =  Js-1m3], <0 for sink
        steps:
            initial subtimesteps used to proceed to 't_final'

    Returns:
        T: new temperature profile [m]
        Wliq: new liquid water content [m3m-3]
        Wice: new ice ontent [m3m-3]
        Fheat: heat flux array [Wm-2]
        R: thermal conductivity [Wm-1K-1]
    CODE:
        Samuli Launiainen, Luke 19.4.2016. Converted from Matlab (APES SoilProfile.Soil_HeatFlow)
    NOTE:
        (19.4.2016): Tested ok; solved 'coupled' water&heatflow solve heatflow and waterflow
        sequentially for small timesteps
    TODO:
        1. think how heat advection with infiltration should be accounted for
        2. gradient boundary should be implemented as flux-based
    """

    Conv_crit1 = 1.0e-3  # degC
    Conv_crit2 = 1e-4  # wliq [m3m-3]

    LH = 0.0    # heat sink/source to 1st node due infiltration or evaporation
                # LH=-cw*infil_rate *T_infil/ dz[0] [Wm-3]
    
    # -------------- constants ---------------
    # grav = GRAVITY  # acceleration due gravity, kgm/s2
    rhoi = 917.0  # ice density, kg/m3
    # NT = NORMAL_TEMPERATURE 
    Lf = LATENT_HEAT_FREEZING  # latent heat of freezing, J/kg; % Lv is latent heat of vaportization
    Tfr = FREEZING_POINT_H2O  # freezing point of water (degC)

    # ------------------- computation grid -----------------------

    # grid
    z = grid['z']
    dz = grid['dz']
    dzu = grid['dzu']
    dzl = grid['dzl']

    N = len(z)

    # -----get parameters from dict
    fp = spara['fp']       # freezing-curve parameter, 2...4 for clay and 0.5-1.5 for sandy soils
                            # (Femma-code/Harri Koivusalo)
#    if type(poros) is float:
#        poros = np.empty(N) + poros

    cs = spara['cs']        # dry soil vol. heat capacity [J m-3 K-1]
    vorg = spara['vOrg']    # organic matter vol. fraction [-]
    vsand = spara['vSand']  # sand vol. fraction [-]
    vsilt = spara['vSilt']  # silt vol. fraction [-]
    vclay = spara['vClay']  # clay vol. fraction [-]
    bedrockL = spara['bedrockL']

    S = np.ones(N)*S  # sink-source array [Wm-3 =  Js-1m3], < 0 is heat sink

    """ -------- find solution at t_final --------- """

    # ---initial conditions
    T = T0
#    Wliq = Wliq0
#    Wice = Wice0
    Wliq, Wice, gamma = frozen_water(T, Wliq0 + Wice0, fp=fp)

    # h0 = wrc(pF, x=Wliq + Wice, var='Th')  # get hydraulic head [m]
    # C = diff_capa(pF, h0)  # dWtot/dh

    dt = t_final / float(steps)  # initial time step [s]
    t = 0.0  # running time

    while t < t_final:  # loop until solution timestep
        # these will stay cduring iteration over time step "dt"
        T_old = T
        # Wliq_old = Wliq
        Wice_old = Wice

        CP_old = volumetric_heat_capacity(poros, Wliq, wice=Wice, cs=cs)  # [Jm-3K-1]

#        R = Cond_deVries(poros, Wliq, wice=Wice, h=h0, pF=pF, T=T, ks=ks)
#        R = spatialAverage(R, method='arithmetic')

        # these change during iteration
        T_iter = T.copy()
        Wliq_iter = Wliq.copy()
        Wice_iter = Wice.copy()

        err1 = 999.0
        err2 = 999.0
        iterNo = 0

        # start iterative solution of heat equation

        while err1 > Conv_crit1 or err2 > Conv_crit2:  # and pass_flag is False: 
            # print 'err1=' +str(err1) +'   err2=' + str(err2)
            iterNo += 1

            CP = volumetric_heat_capacity(poros, Wliq_iter, wice=Wice_iter, cs=cs)
            R = thermal_conductivity(poros, Wliq, wice=Wice, vOrg=vorg, vSand=vsand, vSilt=vsilt, 
                                     vClay=vclay, bedrockL=bedrockL)
            R = spatial_average(R, method='arithmetic')

            # A = Lf**2.0*rhoi / (grav*T_iter+NT)*C  # term due to phase-changes [J m-3 K-1]
            #print gamma
            
            gamma[T_iter > Tfr] = 0.0

            A = -Lf*rhoi*gamma
            A[T_iter > Tfr] = 0.0
            
            # --- set up tridiagonal matrix
            a = np.zeros(N)
            b = np.zeros(N)
            g = np.zeros(N)
            f = np.zeros(N)

            # ---- intermediate nodes
            b[1:-1] = CP[1:-1] + A[1:-1] + dt / dz[1:-1] * (R[1:-1] / dzu[1:-1] + R[2:] / dzl[1:-1])
            a[1:-1] = - dt / (dz[1:-1]*dzu[1:-1]) * R[1:-1]
            g[1:-1] = - dt / (dz[1:-1]*dzl[1:-1]) * R[2:]
            f[1:-1] = CP_old[1:-1] * T_old[1:-1] + A[1:-1] * T_iter[1:-1]\
                      + Lf*rhoi* (Wice_iter[1:-1] - Wice_old[1:-1]) - dt*S[1:-1]

            # ---------- top node (n=0)
            # LH is heat input by infiltration. loss by evaporation not currently implemented

            if ubc['type'] == 'flux':  # or ubc['type'] is 'grad':
                F_sur = ubc['value']
                b[0] = CP[0] + A[0] + dt / (dz[0]*dzl[0]) * R[1]
                a[0] = 0.0
                g[0] = -dt / (dz[0]*dzl[0]) * R[1]
                f[0] = CP_old[0]*T_old[0] + A[0]*T_iter[0] + Lf*rhoi* (Wice_iter[0] - Wice_old[0])\
                       - dt / dz[0]*F_sur - dt*LH - dt*S[0]

            if ubc['type'] == 'temperature':   # fixed T at imaginary node at surface
                T_sur = ubc['value']
                b[0] = CP[0] + A[0] + dt / dz[0]* (R[0] / dzu[0] + R[1] / dzl[0])
                a[0] = 0.0
                g[0] = -dt / (dz[0]*dzl[0]) * R[1]
                f[0] = CP_old[0]*T_old[0] + A[0]*T_iter[0] + Lf*rhoi* (Wice_iter[0] - Wice_old[0])\
                       + dt / (dz[0]*dzu[0])*R[0]*T_sur - dt*LH - dt*S[0]

            # ------ bottom node (n=N)
            if lbc['type'] == 'flux':  # or lbc['type'] is 'grad':
                F_bot = lbc['value']
                b[-1] = CP[-1] + A[-1] + dt / (dz[-1]*dzu[-1]) * R[-2]
                a[-1] = -dt / (dz[-1]*dzu[-1]) * R[-2]
                g[-1] = 0.0
                f[-1] = CP_old[-1]*T_old[-1] + A[-1]*T_iter[-1]\
                        + Lf*rhoi* (Wice_iter[-1] - Wice_old[-1]) - dt / dz[-1]*F_bot - dt*S[-1]
    
            if lbc['type'] == 'temperature':  # fixed temperature, Tbot "at node N+1"
                T_bot = lbc['value']
                b[-1] = CP[-1] + A[-1] + dt / dz[-1] * (R[-2] / dzu[-1] + R[-1] / dzl[-1])
                a[-1] = -dt / (dz[-1]*dzu[-1]) * R[-2]
                g[-1] = 0.0
                f[-1] = CP_old[-1]*T_old[-1] + A[-1]*T_iter[-1]\
                        + Lf*rhoi*(Wice_iter[-1] - Wice_old[-1])\
                        + dt/(dz[-1]*dzl[-1]) * R[-1]*T_bot - dt*S[-1]

            # --------- update iteration values
            T_iterold = T_iter.copy()
            Wliq_iterold = Wliq_iter.copy()
            # Wice_iterold = Wice_iter.copy()

            # --- call tridiagonal solcver
            T_iter = thomas(a, b, g, f)

            Wliq_iter, Wice_iter, gamma = frozen_water(T_iter, Wliq_iter + Wice_iter, fp=fp)

            if iterNo == 7:
                dt = dt / 3.0
                iterNo = 0
                continue  # re-try with smaller timestep
            elif any(np.isnan(T_iter)):
                # print 'soil_heat.heatflow_1D: nan found, breaking'
                break  # break while loop
#            elif np.isnan(Wliq_iter).any():
#                print 'soil_heat.heatflow_1D: nan found in Wliq_iter, breaking'
#                break

            err1 = np.max(abs(T_iter - T_iterold))
            err2 = np.max(abs(Wliq_iter - Wliq_iterold))

        # ------- ending iteration loop -------
        # print 'Pass: t = ' +str(t) + ' dt = ' +str(dt) + ' iterNo = ' + str(iterNo)

        # update state
        T = T_iter.copy()
        Wliq = Wliq_iter.copy()
        Wice = Wice_iter.copy()

        # solution time & new initial timestep
        t += dt

        if iterNo < 2:
            dt = dt*1.25
        elif iterNo > 4:
            dt = dt / 1.25

        dt = min(dt, t_final - t)
        # print 'new dt= ' +str(dt)

    # ----- now at t = t_final, end while loop, compute heat flux profile
    Fheat = nodal_fluxes(z, T, R)  # Wm-2

    return T, Wliq, Wice, Fheat, R


def heatflow_homogenous(t_final, z, To, k, S, ubcType, ubc, lbcType, lbc):
    """
        Solves 1D heat equation in homogenous soil using explicit Eulerian method.
        Forward difference in time, centered in space

        Args:
            t_final: solution time [s]            
            z: grid [m], expects constant increments
            To: initial T profile
            k: thermal diffusivity (m2/s); thermal diffusivity is L/(rho_s*c_s)
            S: heat sink/source 
            ubcType: type of upper bc 'dirchlet' for fixed value, 'neumann' for flux
            ubc: value of upper bc
            lbcType: type of lower bc
            lbc : value of lower bc

        Returns:
            U: new temperature profile

            Samuli Launiainen 14.4.2016
    """
    N = len(z)
    dz = abs(z[1] - z[0])

    dt = 0.4*dz**2.0 / k
    R = k*dt / (dz**2)  # mesh fourier number
    # print dt, R
    U = To.copy()

    t = 0.0
    while t < t_final:
        Uxo = U.copy()
        # upper bc
        if ubcType is 'dirchlet':
            U[0] = ubc
        else:
            U[0] = Uxo[0] + R*(2*U[1] - 2*U[0] - 2*dz*ubc)
        # U[0]=ubc; U[-1] = lbc
        for m in range(1, N-1):
            U[m] = Uxo[m] + R*(U[m-1] - 2*U[m] + U[m+1]) + dt*S[m]
        # lower bc
        if lbcType is 'dirchlet':
            U[-1] = lbc
        else:
            U[-1] = Uxo[-1] + R*(2*U[-2] - 2*U[-1] + 2*dz*lbc)

        t += dt
        # print 't= ', t
    return U


def soil_temperature_Rankinen(t_final, z, To, T_sur, Ds, para, steps=10):
    """
    Soil temperature profile at depths z by simple approximation given in Rankinen et al. (2004).
    Assumes: 
        * depth-constant soil heat capacity and thermal conductivity
        * no-flow conditions below depth in consideration
        * explicit solution
    Not suited for computing realistic temperature profile evolution; use heatflow_1D instead.
    Args:
        t_final: solution time [s]
        z: grid [m]
        To: initial temperature profile [degC]
        T_sur: surface temperature (or air temperature) [degC]
        Ds: snow depth [m]
        para: dictionary (ranges from Rankinen et al. 2004. Table 3:)
            * 'Cs': soil heat capacity Jm-3K-1. Range 1.0e6...1.3e6 for mineral soils, xxx for organic soils
            * 'Cice': latent heat storage/release correction term,
                        accounted for T<0.01degC. Range 4.1e6...9e6; larger for peat soils?
            * 'Ktherm': soil thermal conductivity Wm-1K-1, range 0.5...0.8 Wm-1K-1 (mineral),
                        0.3...0.45 Wm-1K-1 for peat soils
            * 'fs': damping parameter for snow depth impacts. Range 2.1 ...7.1
        steps: number of subtimesteps
    Returns:
        Tnew: new soil temperature profile

    Reference:
        Rankinen et al. 2004. Hydrol. Earth Sys. Sci. 8(4) 706 - 716.
    """
    # -----input params
    R = para['Ktherm']  # Wm-1K-1
    fs = para['fs']  # m-1, snow damping parameter
    Cice = para['Cice']  # apparent ice vol. heat capacity J m-3 K-1
    Cs = para['Cs']  # soil heat capacity Jm-3K-1, range 1e6...1.3e6. LOWER FOR organic soil!

    dt = t_final / float(steps)
    N = len(z)

    Told = To
    Tnew = np.empty(N)
    t = 0.0
    while t < t_final:
        for n in range(0, N):
            if Told[n] <= 0.01:
                Ca = Cs + Cice
            else:
                Ca = Cs

            Tnew[n] = Told[n] + dt*R / (Ca*(2*z[n])**2)*(T_sur - Told[n])
            Tnew[n] = np.exp(-fs*Ds)*Tnew[n]  # snow correction

        Told = Tnew
        t += dt

        # plt.plot(Told,z,'r-')
    return Tnew


""" ----Utility functions ---- """

def find_index(a, func):
    """
    finds indexes or array elements that fill the condition
    call as indices(a, lambda x: criteria)
    """
    return [i for (i, val) in enumerate(a) if func(val)]


def spatial_average(y, x=None, method='arithmetic'):
    """
    calculates spatial average of quantity y
    Args:
        y: variable array
        x: grid array
        method (str): 'arithmetic', 'geometric', 'dist_weighted'
    Returns:
        f: averaged y
    """

    N = len(y)
    f = np.empty(N)
    if method is 'arithmetic':
        f[1:-1] = 0.5*(y[0:-2] + y[1:-1])
        f[0] = y[0]
        f[-1] = y[-1]

    elif method is 'geometric':
        f[1:-1] = np.sqrt(y[1:-2]*y[2:])
        f[0] = y[0]
        f[-1] = y[-1]

    elif method is 'dist_weighted':
        a = (x[0:-2] - x[2:]) * y[:-2]*y[1:-1]
        b = y[1:-1] * (x[:-2] - x[1:-1]) + y[:-2]*(x[1:-1] - x[2:])

        f[1:-1] = a / b
        f[0] = y[0]
        f[-1] = y[-1]

    return f

def nodal_fluxes(x, y, K):
    """
    Calculates fluxes between nodal points in 1-D grid
    Args:
        x: grid array, monotonic!
        y: variable array
        K: conductivity, float or array
    Returns:
        f: flux array
    """
    f = np.empty(np.shape(y))
    yprim = central_difference(x, y)

    f = -K*yprim  # flux
    return f


def central_difference(x, y):
    """
    derivative by central difference
    IN:
        x - grid
        y - value vector
    OUT:
        yprim - derivative
    """
    yprim = np.empty(np.shape(y))

    yprim[1:-1] = (y[2:] - y[0:-2]) / (x[2:] - x[0:-2])  # central difference
    yprim[0] = (y[1] - y[0]) / (x[1] - x[0])  # forward difference at left boundary
    yprim[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])  # backward difference at right boundary

    return yprim

def soilRH(T, Psi):
    """ 
    Soil relative humidity (-) according Philip and de Vries, 1957.
    Args:
        float or np.array:
            T: temperature (degC)
            Psi: soil water tension (m H2O)

    Returns:
        float or np.array:
            relative humidity [-]
    Note: accepts scalar of list inputs
    """   
    # relative humidity in soil pores [-]
    rh = np.exp(-(MOLAR_MASS_H2O*GRAVITY*abs(Psi)) / (GAS_CONSTANT*(T + 273.15)))
    return rh


def e_sat(T):
    """"
    Saturation vapor pressure [Pa] over water surface and its derivative
    d(es)/dT [Pa K-1]
    Args:
        float/np.array: T [degC]
    Returns:
        float/np.array:
        * es [Pa]
        * d(es)/dT [Pa K-1]
    """
    if isinstance(T, float) == False:
        T = np.array(T)
    esat = 611.0*np.exp((17.502*T)/(T + 240.97))  # (Stull 1988)
    ss = 17.502*240.97*esat / ((240.97 + T)**2)  # [Pa K-1], from Baldocchi's Biometeorology -notes
    return esat, ss
    

def frozen_water(T, wtot, fp=2.0, To=0.0):
    """
    Approximates ice content from soil temperature and total water content
    Args:
        T: soil temperature [degC]
        wtot: total vol. water content [m3 m-3]
        fp: parameter of freezing curve [-];  2...4 for clay and 0.5-1.5 for sandy soils, 
            <0.5 for peat soils (Nagare et al. 2012 HESS)
        To: freezing temperature of soil water [degC]
    Returns:
        wliq: vol. water content [m3 m-3]
        wice: vol. ice content [m3 m-3]
    For peat soils, see experiment of Nagare et al. 2012:
    http://scholars.wlu.ca/cgi/viewcontent.cgi?article=1018&context=geog_faculty
    """
    #if np.isnan(wtot).any():
        #print('NaN encountered in soil_heat wtot in frozen water')

    wtot = np.array(wtot, ndmin=1)
    T = np.array(T, ndmin=1)

    wice = wtot*(1.0 - np.exp(-(To - T) / fp))
    # derivative dwliq/dT
    gamma = (wtot - wice) / fp

    wice[T > To] = 0.0
    gamma[T > To] = 0.0

    wliq = wtot - wice

    return wliq, wice, gamma

def spatial_average2(y, x=None, method='arithmetic'):
    """
    Calculates spatial average of quantity y, from node points to soil compartment edges
    Args: 
        y (array): quantity to average
        x (array): grid,<0, monotonically decreasing [m]
        method (str): flag for method 'arithmetic', 'geometric','dist_weighted'
    Returns: 
        f (array): averaged y, note len(f) = len(y) + 1
    """

    N = len(y)
    f = np.empty(N+1)  # Between all nodes and at surface and bottom
    if method is 'arithmetic':
        f[1:-1] = 0.5*(y[:-1] + y[1:])
        f[0] = y[0]
        f[-1] = y[-1]

    elif method is 'geometric':
        f[1:-1] = np.sqrt(y[:-1] * y[1:])
        f[0] = y[0]
        f[-1] = y[-1]

    elif method is 'dist_weighted':                                             # En ymmärrä, ei taida olla käyttössä
        a = (x[0:-2] - x[2:])*y[:-2]*y[1:-1]
        b = y[1:-1]*(x[:-2] - x[1:-1]) + y[:-2]*(x[1:-1] - x[2:])

        f[1:-1] = a / b
        f[0] = y[0]
        f[-1] = y[-1]

    return f