# -*- coding: utf-8 -*-
"""
.. module: interception
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes interteption in multilayer canopy.
Based on MatLab implementation by Samuli Launiainen.

Created on Thu Mar 01 13:21:29 2018

Note:
    migrated to python3
    - absolute imports
    - in for-loop range() does not need to be wrapped in list()

References:
Tanaka, K., 2002. Multi-layer model of CO2 exchange in a plant community
coupled with the water budget of leaf surfaces. Ecological Modelling, 147(1), pp.85-104.

Last edit:
SL 12.11.2019: removed correction factors for rain and snow. State of Prec
now computed based on T at highest gridpoint.

SL 13.1.2020: changed units from [m] to [kg m-2] and [m s-1] to [kg m-2 s-1]

"""

import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from .micromet import e_sat, leaf_boundary_layer_conductance
from .constants import WATER_DENSITY, MOLAR_MASS_H2O, SPECIFIC_HEAT_AIR, DEG_TO_KELVIN, EPS

import logging
logger = logging.getLogger(__name__)

class Interception(object):
    r"""Describes interception in multilayer canopy.
    """
    def __init__(self, p, LAIz):
        r""" Initializes interception object.
        Args:
            p (dict):
                'wmax': maximum interception storage capacity for rain [kg water m-2 per unit of LAI]
                'wmaxsnow': maximum interception storage capacity for snow [kg water m-2 per unit of LAI]
                'Tmin': temperature below which all is snow [degC]
                'Tmax': temperature above which all is water [degC]
                'w_ini': initial canopy storage [kg water m-2(ground)]
            LAIz (array): leaf area index per canopy layer [m\ :sup:`2`\ m\ :sup:`-2`\]
        Returns:
            self (object)
        """
        # parameters:
        # maximum storage capacities [m H2O per unit of LAI]
        self.wmax = p['wmax']  # for rainfall
        self.wmaxsnow = p['wmaxsnow']  # for snowfall
        # quality of precipitation [degC]
        self.Tmin = p['Tmin']
        self.Tmax = p['Tmax']

        #self.c_rain = p['c_rain']
        #self.c_snow = p['c_snow']

        # Leaf orientation factor with respect to incident Prec (horizontal leaves -> 1)
        self.leaf_orientation = p['leaf_orientation']

        # initial state [kg m-2 s-1]
        self.W = np.minimum(p['w_ini'], p['wmax'] * LAIz)

        self.update()

    def run(self, dt, forcing, parameters, controls):
        r""" Computes interception and unloading of rain or snow,
        evaporation and condensation are computed based on wet leaf water balance.

        Rate of change of water stored at each canopy layer (W) is as in Tanaka, 2002:
            (1) dW(z)/dt = F(1-W(z)/Wmax(z))I(z) - (W(z)/Wmax(z))E(z), I(z)=dPrec/dz=Fa(z)dz(1-W(z)/Wmax(z))Prec(z+dz) when E>0 (evaporation)
            (2) dW(z)/dt = (1-W(z)/Wmax(z))(FI(z)-E(z)), I(z)=dPrec/dz=Fa(z)dz(1-W(z)/Wmax(z))Prec(z+dz) + a(z)dz(W(z)/Wmax(z))E(z) when E<0(condensation)
        a(z) is one-sided plant area density (m2m-3), F fraction of leaf area (0-1)
        perpendicular to Prec and Wmax(z) maximum water storage (mm) at each layer.

        Args:
            dt: timestep [s]
            forcing (dict):
                'net_lw_leaf' (array): net radiation balance at each layer [W m\ :sup:`-2`\]
                'sw_absorbed' (array): absorbed shortwave radiation at each layer [W m\ :sup:`-2`\]
                'precipitation' (float): precipitation rate above canopy [kg  s\ :sup:`-1`\]
                'air_pressure' (float): ambient pressure [Pa]
                'leaf_temperature' (array): average leaf temperature used in LW computation [degC]
                'radiative_conductance' (array): radiative conductance [mol m-2 s-1]
                'h2o' (array): [mol mol-1]
                'wind_speed' (array): [m s-1]
                'air_temperature' (array): [degC]
            parameters (dict):
                LAIz (array): leaf area index per canopy layer [m\ :sup:`2`\ m\ :sup:`-2`\]
                leaf_length (array): leaf length scale for aerodynamic resistance [m]
            controls (dict):
                'energy_balance': boolean

        Returns:
            'throughfall' (float): total throughfall [kg m-2 s-1], multiply by WATER_DENSITY to get [kg m-2 s-1]
            'throughfall_rain' (float): rain throughfall [kg m-2 s-1]
            'throughfall_snow' (float): snow throughfall [kg m-2 s-1]
            'interception' (float): [kg m-2 s-1]
            'evaporation'(float): [kg m-2 s-1]
            'condensation' (float): [kg m-2 s-1]
            'condensation_drip' (float): [kg m-2 s-1]
            'water_closure' (float): [kg m-2 s-1]
            'sources': {'h2o': dqsource,
                        'sensible_heat': Heat / dt,
                        'fr': Fr / dt,
                        'latent_heat': dqsource * L},
            'evaporation_ml' (array): evaporation/condensation rate in layers [kg m-2 s-1]
            'throughfall_ml' (array): throughfall rate in layers [kg m-2 s-1]
            'condensation_drip_ml' (array): condensation drip rate in layers [kg m-2 s-1]

        """
        lt = np.maximum(EPS, parameters['leaf_length'])

        Prec = forcing['precipitation']
        P = forcing['air_pressure']
        Tl_ave = forcing['leaf_temperature']
        H2O = forcing['h2o']
        U = forcing['wind_speed']
        T = forcing['air_temperature']
        LAIz = parameters['LAIz']

        Ebal = controls['energy_balance']

        if controls['energy_balance']:
            gr = forcing['lw_radiative_conductance']
            Rabs = (forcing['sw_absorbed'] +
                    forcing['net_lw_leaf'])
        else:
            gr = 0.0

        # number of canopy layers
        N = len(LAIz)
        ic = np.where(LAIz > 0)

        # initial guess for wet leaf temperature
        Tl_wet = self.Tl_wet.copy()
        Told = Tl_wet.copy()

        # latent heat of vaporization/sublimation at temperature T [J/mol]
        L = latent_heat(T) * MOLAR_MASS_H2O

        # Leaf orientation factor with respect to incident Prec (horizontal leaves -> 1)
        F = self.leaf_orientation

        # --- state of precipitation (based on T at highest gridpoint)---
        # fraction as water [-]
        if T[-1] < self.Tmin:
            fW = 0.0
        else:
            fW = np.minimum(1.0, (T[-1] - self.Tmin) / (self.Tmax - self.Tmin))

#        fW = np.ones(N)
#        ix = np.where(T < self.Tmin)
#        fW[ix] = 0.0
#        ix = np.where((T >= self.Tmin) & (T <= self.Tmax))
#        fW[ix] = (T[ix] - self.Tmin) / (self.Tmax - self.Tmin)

        # correction of precipitation
        #Prec = Prec * fW[-1] * self.c_rain + Prec * (1 - fW[-1]) * self.c_snow

        # maximum interception storage capacities layerwise [m]
        Wmax = (fW * self.wmax + (1 - fW) * self.wmaxsnow) * LAIz + eps

#        # boundary layer conductances for H2O and heat [mol m-2 s-1]
#        gb_h, _, gb_v = leaf_boundary_layer_conductance(U, lt, T, 0.0, P)

        # vapor pressure deficit between leaf and air, and slope of vapor pressure curve at T
        es, s = e_sat(Tl_wet)
        Dleaf = es / P - H2O  #np.maximum(0.0, es / P - H2O)  # [mol/mol]
        s = s / P  # [mol/mol/degC]

        """ --- wet Leaf temperature --- """
        itermax = 20
        err = 999.0
        iterNo = 0
        while err > 0.01 and iterNo < itermax:
            iterNo += 1

            # boundary layer conductances for H2O and heat [mol m-2 s-1]
            gb_h, _, gb_v = leaf_boundary_layer_conductance(U, lt, T, 0.5*(Tl_wet + Told) - T, P)  # OK to assume dt = 0.0?? convergence problems otherwise

            Told = Tl_wet.copy()

            if Ebal:
                # solve leaf temperature [degC]
                Tl_wet[ic] = (Rabs[ic] + SPECIFIC_HEAT_AIR*gr[ic]*Tl_ave[ic] + SPECIFIC_HEAT_AIR*gb_h[ic]*T[ic] - L[ic]*gb_v[ic]*Dleaf[ic]
                  + L[ic]*s[ic]*gb_v[ic]*Told[ic]) / (SPECIFIC_HEAT_AIR*(gr[ic] + gb_h[ic]) + L[ic]*s[ic]*gb_v[ic])
                err = np.nanmax(abs(Tl_wet - Told))

                if (err < 0.01 or iterNo == itermax) and abs(np.mean(T) - np.mean(Tl_wet)) > 20.0:
                    logger.debug(controls['logger_info'] + ',%s Unrealistic wet leaf temperature %.2f set to air temperature %.2f, %.2f, %.2f, %.2f',
                         iterNo,
                         np.mean(Tl_wet), np.mean(T),
                         np.mean(forcing['net_lw_leaf']), np.mean(Tl_ave), np.mean(self.Tl_wet))
                    Tl_wet = T.copy()

                elif iterNo == itermax:
                    logger.debug(controls['logger_info'] + ' Maximum number of iterations reached: Tl_wet = %.2f, err = %.2f',
                             np.mean(Tl_wet), err)

                es, s = e_sat(Tl_wet)
                Dleaf = es / P - H2O  #np.maximum(0.0, es / P - H2O)  # [mol/mol]
                s = s / P  # [mol/mol/degC]

            else:
                err = 0.0

        # --- energy and water fluxes for wet leaf ---
        # or sublimation/deposition ????????? GITHUB SPATHY!!

        # sensible heat flux [W m-2(wet leaf)]
        Hw = SPECIFIC_HEAT_AIR * gb_h * (Tl_wet - T)
        # non-isothermal radiative flux [W m-2 (wet leaf)]
        Frw = SPECIFIC_HEAT_AIR * gr *(Tl_wet - Tl_ave)
        # evaporation rate from wet leaf [kg m-2 s-1] (negative for condensation)
        Ep = gb_v * Dleaf * MOLAR_MASS_H2O

        # Assume no evaporation during rain when energy balance not solved
        if Ebal == False and Prec > 0.0:
            Ep = np.zeros(N)

        # --- canopy water storage change ---
        W = self.oldW.copy()  # layerwise canopy storage [kg m-2 s-1]

        # Unloading in canopy, ensures also that seasonal
        # LAI development does not mess up computations
        for n in reversed(range(N)):  # start from highest grid point
            Unload = max(W[n] - Wmax[n], 0.0)  # unloading from layer n
            W[n] -= Unload  # update storage of layer n
            if n != 0:
                W[n-1] += Unload  # unloading added to layer below (if not lower layer)
        # Unload = unloading below canopy [kg m-2]

        # timestep subdivision to calculate change in canopy water store, no impact??
        Nsteps = 1  # number of subtimesteps
        subdt = dt / Nsteps  # [s]

        # initiate cumulative variables
        Interc = np.zeros(N)  # interception [kg m-2]
        Evap = np.zeros(N)  # evaporation [kg m-2]
        Cond = np.zeros(N)  # condesation [kg m-2]
        Heat = np.zeros(N)  # sensible heat flux [W m-2(ground)]
        Fr = np.zeros(N)  # sensible heat flux [W m-2(ground)]
        wf = np.zeros(N)  # wetness ratio
        Tr = np.zeros(N)  # throughfall within canopy [kg m-2]
        Trfall = 0.0  # throughfall below canopy [kg m-2]

        if Prec > 0 or np.any(np.less(Ep, 0)) or np.any(np.greater(W, 0)):
            for t in range(Nsteps):
                Ir = np.zeros(N)  # interception rate
                dW = np.zeros(N)  # change in storage [m]
                P = np.zeros(N+1)  # precipitation rate to layer
                P[-1] = Prec  # above canopy equals precipitation rate 
                for n in reversed(range(N)):  # start from highest grid point
                    if Ep[n] >= 0:  # evaporation case
                        # change in storage
                        dW[n] = (F * P[n+1] / (F * P[n+1] + Ep[n] + eps) * Wmax[n] - W[n]) \
                                * (1.0 - np.exp(-(F * P[n+1] + Ep[n]) * LAIz[n] * subdt / Wmax[n]))
                        # wetness ration in layer
                        if LAIz[n] > 0 and P[n+1] + Ep[n] > 0:
                            wf[n] = (F * P[n+1] - dW[n] / (LAIz[n] * subdt)) / (F * P[n+1] + Ep[n])
                        else:
                            wf[n] = 0.0
                        # interception rate in layer
                        Ir[n] = F * (1 - wf[n]) * LAIz[n] * P[n+1]
                        # drainage rate from layer
                        P[n] = P[n+1] - Ir[n]
                        # evaporation from layer
                        Evap[n] += wf[n] * LAIz[n] * Ep[n] * subdt
                    else:  # condensation case
                        # change in storage
                        dW[n] = (Wmax[n] - W[n]) \
                                * (1.0 - np.exp(-(F * P[n+1] - Ep[n]) * LAIz[n] * subdt / Wmax[n]))
                        # wetness ration in layer
                        if LAIz[n] > 0 and P[n+1] - Ep[n] > 0:
                            wf[n] = (F * P[n+1] - Ep[n] - dW[n] / (LAIz[n] * subdt)) / (F * P[n+1] - Ep[n])
                        else:
                            wf[n] = 0.0
                        # interception rate in layer
                        Ir[n] = F * (1 - wf[n]) * LAIz[n] * P[n+1]
                        # drainage rate from layer (incl. increase by condensation drip)
                        P[n] = P[n+1] - Ir[n] - wf[n] * LAIz[n] * Ep[n]
                        # Condensation (incl. condenstation to dry leaf and drip from wet leaf)
                        Cond[n] += LAIz[n] * Ep[n] * subdt
                
                # ! condensation to dry leaf part not accounted for in energy balance here but in dry leaf module
                # Sensible heat flux [W m-2(ground)] * subdt
                Heat += wf * LAIz * Hw * subdt
                # radiative flux [W m-2(ground)] * subdt
                Fr += wf * LAIz * Frw * subdt
#                LE += wf * LAIz * Ep / MOLAR_MASS_H2O * WATER_DENSITY * L * subdt
                # update storage
                W += dW
                # interception and throughfall
                Interc += Ir * subdt
                Trfall += P[0] * subdt
                Tr += P[:-1] * subdt

        # throughfall to field layer or snowpack
        Trfall = Trfall + Unload
        Trfall_rain = fW * Trfall 
        Trfall_snow = (1 - fW) * Trfall

        # H20 source/sink per ground area due to evaporation and condensation [mol m-2 s-1]
        dqsource = (Evap + Cond) / dt / MOLAR_MASS_H2O

        if sum(W) < eps:
            W *= 0.0

        # dry canopy fraction
#        df = 1.0 - np.where(Ep >= 0, wf, 1.0)
        df = 1.0 - wf

        # update state variables
        self.W = W
        self.Tl_wet = Tl_wet
        self.df = df

        # mass-balance error [kg m-2 s-1] ! self.W is old storage
        water_closure = sum(self.W) - sum(self.oldW) - (Prec * dt - sum(Evap) - sum(Cond) - (Trfall_rain + Trfall_snow))

        fluxes = {
                  'throughfall': (Trfall_rain + Trfall_snow) / dt,
                  'throughfall_rain': Trfall_rain / dt,
                  'throughfall_snow': Trfall_snow / dt,
                  'interception': sum(Interc) / dt,
                  'evaporation': sum(Evap) / dt,
                  'condensation': sum(Cond * (1 - wf)) / dt,
                  'condensation_drip': sum(Cond * wf) / dt,
                  'water_closure': water_closure / dt,
                  'sources': {'h2o': dqsource,
                              'sensible_heat': Heat / dt,
                              'fr': Fr / dt,
                              'latent_heat': dqsource * L},
                  'evaporation_ml': (Evap + Cond * (1 - wf)) / dt,
                  'throughfall_ml': Tr / dt,
                  'condensation_drip_ml': (Cond * wf) / dt
                  }
        return fluxes

    def update(self):
        """Updates interception storage W to old W
        """

        self.oldW = self.W.copy()

def latent_heat(T):
    """
    Computes latent heat of vaporization or sublimation [J/kg]
    Args:
        T: ambient air temperature [degC]
    Returns:
        L: latent heat of vaporization or sublimation depending on 'type'[J/kg]
    """
    # latent heat of vaporizati [J/kg]
    Lv = 1e3 * (3147.5 - 2.37 * (T + DEG_TO_KELVIN))
    # latent heat sublimation [J/kg]
    Ls = Lv + 3.3e5
    L = np.where(T < 0, Ls, Lv)
    return L
