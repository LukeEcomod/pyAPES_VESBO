# -*- coding: utf-8 -*-
"""
.. module: soilprofile.heat
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Represents soil heat balance.

Created on Thu Oct 04 09:04:05 2018
"""

import numpy as np
import matplotlib.pyplot as plt
from water import wrc
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
K_MINERAL = 2.9  # thermal conductivity of other minerals W/m/K
K_ORG = 0.25  # thermal conductivity of organic matter W/m/K
K_SAND = 7.7 # thermal conductivity of sand W/m/K. Tian et al. 2007
K_SILT = 2.74 # thermal conductivity of silt W/m/K. Tian et al. 2007
K_CLAY = 1.93 # thermal conductivity of clay W/m/K. Tian et al. 2007

#: volumetric heat capacieties  [J m\ :sup:`-3`\ K \ :sup:`-1`\]
CV_AIR = 1297.0  # air at 101kPa [Jm-3K-1]
CV_WATER = 4.18e6  # water [Jm-3K-1]
CV_ICE = 1.93e6  # ice [Jm-3K-1]
CV_ORGANIC = 2.50e6  # dry organic matter
CV_MINERAL = 2.31e6  # soil minerals
CV_GRANITE = 2.16e6  # granite
CV_QUARTZ = 2.13e6  # quartz minerals

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

class HeatModel(object):
    r""" Represents soil heat balance.
    """
    def __init__(self, grid, profile_propeties, model_specs, volumetric_water_content):
        r""" Initializes soil heat balance model.

        Args:
            grid (dict):
                'z': node elevations, soil surface = 0.0 [m]
                'dz': thickness of computational layers [m]
                'dzu': distance to upper node [m]
                'dzl': distance to lower node [m]
            profile_propeties (dict): arrays of len(z)
                'porosity' (array): soil porosity (=ThetaS) [m\ :sup:`3` m\ :sup:`-3`\ ]
                'solid_heat_capacity' (array): [J m-3 (solid volume) K-1]
                    ! if nan, estimated from organic/mineral composition
                'fractions' (dict): fractions of solid volume [-]
                    'organic' (array)
                    'sand' (array)
                    'silt' (array)
                    'clay' (array)
                'freezing_curve' (array): freezing curve parameter [-]
                'bedrock_thermal_conductivity' (array): nan above bedrock depth[W m-1 K-1]
            model_specs (dict):
                'solve': True/False
                'initial_condition':
                    'temperature' (array/float): initial temperature [degC]
                'lower_boundary':
                    'type' (str): 'flux' or 'temperature'
                    'value' (float): value of flux [W m-2] or temperature [degC]
        Returns:
            self (object)
        """
        # lower boundary condition
        self.lbc = model_specs['lower_boundary']

        # grid
        self.grid = grid
        # dummy
        self.zeros = np.zeros(len(self.grid['z']))

        # profile parameters
        self.porosity = profile_propeties['porosity']
        self.bedrock_thermal_conductivity = profile_propeties['bedrock_thermal_conductivity']
        self.solid_composition = profile_propeties['fractions']
        self.freezing_curve = profile_propeties['freezing_curve']
        self.solid_heat_capacity = profile_propeties['solid_heat_capacity']
        ix = np.where(np.isnan(self.solid_heat_capacity))[0]
        self.solid_heat_capacity[ix] = solid_volumetric_heat_capacity(self.solid_composition['organic'])

        # initialize state
        self.update_state(model_specs['initial_condition'], volumetric_water_content)

        # model specs
        if model_specs['solve']:
            # Keep track of dt used in solution
            self.subdt = None

    def run(self, dt, forcing, heat_sink=None, lower_boundary=None):
        r""" Runs soil heat balance.
        Args:
            dt: time step [s]
            forcing (dict):
                'ground_heat_flux' (float): heat flux from soil surface [W m-2]
                    OR 'temperature' (float): soil surface temperature [degC]
            volumetric_water_content (array):  [m3 m-3]
            heat_sink (array or None): layerwise sink term [W m-3]
                ! if None set to zero
            lower_boundary (dict or None): allows boundary to be changed in time
                'type' (str): 'flux' or 'temperature'
                'value' (float): value of flux [W m-2] or temperature [degC]
        Returns:
            fluxes (dict):[W m-2]
                'energy_closure'
        """
        if self.subdt is None:
            self.subdt = dt

        if lower_boundary is not None:  # allows lower boundary to change in time
            self.lbc = lower_boundary

        if 'ground_heat_flux' in forcing:
            upper_boundary ={'type': 'flux', 'value': forcing['ground_heat_flux']}
        else:
            upper_boundary ={'type': 'temperature', 'value': forcing['temperature']}

        fluxes, state = heatflow1D(t_final=dt,
                                   grid=self.grid,
                                   T_ini=self.T,
                                   Wtot=volumetric_water_content,
                                   poros=self.porosity,
                                   ubc=upper_boundary,
                                   lbc=self.lbc,
                                   spara=self.solids_prop,  # dict of properties of solid part
                                   S=heat_sink,
                                   steps= dt / self.subdt)

        self.update_state(state)

        return fluxes

    def update_state(self, state={}, Wtot=0.0):
        r""" Updates state to HeatModel object.
        Args:
            state (dict):
                'temperature' (float):  [degC]
        Returns:
            updates .T, .Wice, .Wliq, .thermal_conductivity
        """
        if 'temperature' in state:
            self.T = state['temperature']
        if 'volumetric_ice_content' in state:
            self.Wice = state['volumetric_ice_content']
        else:
            _, self.Wice, _ = frozen_water(self.T, Wtot, self.freezing_curve, To=FREEZING_POINT_H2O)
        if 'volumetric_liquid_water_content' in state:
            self.Wliq = state['volumetric_liquid_water_content']
        else:
            self.Wliq = Wtot - self.Wice
        self.thermal_conductivity = thermal_conductivity(self.porosity, self.Wliq, self.Wice,
                                                         solid_composition=self.solid_composition,
                                                         bedrockL=self.bedrock_thermal_conductivity)

def solid_volumetric_heat_capacity(f_org):
    r""" Computes volumetric heat capacity of solids in soil
    based on organic/mineral compsition.
    Args:
        f_org: fraction of organic matter of solid volume [-]
    Returns:
        cv_solids: volumetric heat capacity of solids [J m-3 (solids) K-1]
    """
    cv_solids = f_org * CV_ORGANIC + (1. - f_org) * CV_MINERAL

    return cv_solids

def volumetric_heat_capacity(poros, cv_solids, wliq=0.0, wice=0.0):
    r""" Computes volumetric heat capacity of soil.
    Args:
        poros - porosity [m3 m-3]
        cv_solids - heat capacity of solids [J m-3(solids) K-1]
        wliq - vol. liquid water content [m3 m-3]
        wice - vol. ice content [m3 m-3]
    Returns:
        cv - volumetric heat capacity of soil [J m-3 K-1]
    """
    wair = poros - wliq - wice

    cv = cv_solids*(1. - poros) + CV_WATER * wliq + CV_ICE * wice + CV_AIR * wair

    return cv

def thermal_conductivity(poros, wliq, wice, solid_composition, bedrockL):
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
    f_sand = solid_composition['sand']
    f_silt = solid_composition['silt']
    f_clay = solid_composition['clay']
    f_org = solid_composition['organic']
    
    poros = np.array(poros, ndmin=1)
    wliq = np.array(wliq, ndmin=1)
    wice = np.array(wice, ndmin=1)
    wair = poros - wliq - wice

    # component thermal conductivities W m-1 K-1
    ks = K_SAND**f_sand * K_SILT**f_silt * K_CLAY**f_clay * K_ORG**f_org
    kw = K_WATER
    ki = K_ICE
    kg = K_AIR

    # shape factors
    g_sand = 0.182
    g_silt = 0.0534
    g_clay = 0.00775
    g_org = 0.33  # campbell & norman, 1998
    g_min = g_sand * f_sand + g_silt * f_silt + g_clay * f_clay + g_org * f_org

    g_air = 0.333 * (1.0 - wair / poros)
    g_ice = 0.333 * (1.0 - wice / poros)

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
    L[f_org >= 0.9] = 0.032 + 5e-1 * wliq[f_org >= 0.9]

    # in bedrock
    L[~np.isnan(bedrockL)] = bedrockL[~np.isnan(bedrockL)]

    return L

""" Soil heat transfer in 1D """

def heatflow1D(t_final, grid, poros, T0, Wtot, ubc, lbc, spara, S=0.0, steps=10):
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

    Code:
        Samuli Launiainen, Luke 19.4.2016. Converted from Matlab (APES SoilProfile.Soil_HeatFlow)
        Kersti Haahti, Luke 26.7.2018
            edits to avoid computation getting stuck during freezing conditions,
            still problems with convergence around zero temperatures.
            Check heat balance error to see problematic periods
    Notes:
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
        CP_old = volumetric_heat_capacity(poros, Wliq, wice=Wice, cs=cs, vOrg=vorg)

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
            CP = volumetric_heat_capacity(poros, Wliq, wice=Wice, cs=cs, vOrg=vorg)
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
        if ubc['type'] == 'temperature':
            Fheat[0] += -R[0]*(T[0] - T_sur)/dzu[0] * dt / t_final
        if ubc['type'] == 'flux':
            Fheat[0] += -F_sur * dt / t_final
        if lbc['type'] == 'temperature':
            Fheat[-1] += -R[-1]*(T_bot - T[-1])/dzl[-1] * dt / t_final
        if lbc['type'] == 'flux':
            Fheat[-1] += -F_bot * dt / t_final

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
        Cv = volumetric_heat_capacity(poros, wliq=wliq, wice=Wice, cs=cs, vOrg=vorg)
        return Cv * x - Lf * rhoi * wice

    def heat_balance(x):
        return heat_content(T0) + t_final * (Fheat[:-1] - Fheat[1:]) / dz - heat_content(x)

    heat_be = sum(dz * heat_balance(T))

    fluxes = {'energy_closure': heat_be / t_final}

    state = {'temperature': T,
             'volumetric_ice_content': Wice}

    return fluxes, state, dto


""" utility functions """

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

    wtot = np.array(wtot)
    T = np.array(T)
    fp = np.array(fp)

    wice = wtot*(1.0 - np.exp(-(To - T) / fp))
    # derivative dwliq/dT
    gamma = (wtot - wice) / fp

    ix = np.where(T > To)[0]
    wice[ix] = 0.0
    gamma[ix] = 0.0

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