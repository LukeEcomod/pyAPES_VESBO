# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 13:21:29 2018

@author: L1656
"""

import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from micromet import e_sat, leaf_boundary_layer_conductance
from constants import WATER_DENSITY, MOLAR_MASS_H2O, SPECIFIC_HEAT_AIR, DEG_TO_KELVIN

class Interception():
    """

    """
    def __init__(self, p, LAIz):
        """
        Args:
            dt: timestep [s]
            p - parameter dict
                'wmax':  maximum interception storage capacity for rain [m per unit of LAI]
                'wmaxsnow': maximum interception storage capacity for snow [m per unit of LAI]
                'Tmin': temperature below which all is snow [degC]
                'Tmax': temperature above which all is water [degC]
                'w_ini': initial canopy storage [m]
        Returns:
            Interception -object
        """
        # parameters:
        # maximum storage capacities [m per unit of LAI]
        self.wmax = p['wmax']  # for rainfall
        self.wmaxsnow = p['wmaxsnow']  # for snowfall
        # quality of precipitation [degC]
        self.Tmin = p['Tmin']
        self.Tmax = p['Tmax']
        self.Tmelt = p['Tmelt']

        # initial state
        self.W = np.minimum(p['w_ini'], p['wmax'] * LAIz)
        self.Tl_wet = None

        self._update()

    def _multi_layer(self, dt, lt, ef, LAIz, H2O, U, T, Rabs, Prec, P, Ebal, Tl_ave, gr):
        """
        updates canopy liquid water storage by partitioning precipitation (Prec) to
        interception and throughfall (rain and snow separately) at each layer
        until soil surface or depth of snow. 
        T_leaf assumed equal to T_air.
        boundary layer conductances for heat and H2O are calculated from U(z) 
        and leaf dimensions.
        Rate of change of water stored at each layer (W) is as in 
        Watanabe & Mizutani (1996) Agric. For. Meteorol. 80, 195-214 
        Tanaka, 2002 Ecol. Mod. 147: 85-104
            (1) dW(z)/dt = F(1-W(z)/Wmax(z))I(z) - (W(z)/Wmax(z))E(z), I(z)=dPrec/dz=Fa(z)dz(1-W(z)/Wmax(z))Prec(z+dz) when E>0 (evaporation)
            (2) dW(z)/dt = (1-W(z)/Wmax(z))(FI(z)-E(z)), I(z)=dPrec/dz=Fa(z)dz(1-W(z)/Wmax(z))Prec(z+dz) + a(z)dz(W(z)/Wmax(z))E(z) when E<0(condensation)
        a(z) is one-sided plant area density (m2m-3), F fraction of leaf area (0-1)
        perpendicular to Prec and Wmax(z) maximum water storage (mm) at each layer.

        Args:
            dt: timestep [s]
            lt: leaf length scale for aerodynamic resistance [m]
                NOTE: seems that lt order of 1-10 cm is appropriate for pine evaporation
            ef: leaf emissivity
            LAIz (array): leaf area index per canopy layer [m\ :sup:`2`\ m\ :sup:`-2`\]
            H2O (array): mixing ratio for each layer [mol/mol]
            U (array): wind speed for each layer [m s-1]
            T (array): air temperature for each layer [degC]
            Rabs (array): isothermal radiation balance at each layer [W m\ :sup:`-2`\]
            Prec (float): precipitation rate above canopy [m s\ :sup:`-1`\]
            P (float): ambient pressure [Pa]
        Returns:
            updates sate self.W (array)
            LEwc: latent heat flux from wet fraction [W m-2], sum gives canopy value
            df: fraction of dry leaves per layer [-]
            Trfall_rain: throughfall as rainfall [m]
            Trfall_snow: throughfall as snowfall [m]
            Interc: interception [m]
            Evap: evaporation from canopy store [m]
            MBE: mass balance error
        
        sources: APES WetLeaf_Module2()
        """

        # number of canopy layers
        N = len(LAIz)
        ic = np.where(LAIz > 0)

        # initial guess for wet leaf temperature
        if self.Tl_wet is None or Ebal is False:
            Tl_wet = T.copy()
        else:
            Tl_wet = self.Tl_wet.copy()

        # radiative conductance mol m-2 s-1, Campbell & Norman, 1998
        gr = gr / SPECIFIC_HEAT_AIR

        # latent heat of vaporization/sublimation at temperature T [J/mol]
        L = latent_heat(T) * MOLAR_MASS_H2O

        # Leaf orientation factor with respect to incident Prec; assumed to be 1 when Prec is in vertical
        F = 0.5 # for randomdly oriented leaves

        """ --- state of precipitation (uses fW[-1] in end of code)--- """
        # fraction as water [-]
        fW = np.ones(N)
        ix = np.where(T < self.Tmin)
        fW[ix] = 0.0
        ix = np.where((T >= self.Tmin) & (T <= self.Tmax))
        fW[ix] = (T[ix] - self.Tmin) / (self.Tmax - self.Tmin)

        # maximum interception storage capacities layerwise [m]
        Wmax = (fW * self.wmax + (1 - fW) * self.wmaxsnow) * LAIz + eps

        # boundary layer conductances for H2O and heat [mol m-2 s-1]
        gb_h, _, gb_v = leaf_boundary_layer_conductance(U, lt, T, 0.0, P)  # OK to assume dt = 0.0?? convergence problems otherwise
#        gb_h, _, gb_v = leaf_boundary_layer_conductance(U, lt, T, Tl_wet - T, P)

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
            Told = Tl_wet.copy()

            if Ebal:
                # solve leaf temperature [degC]
                Tl_wet[ic] = (Rabs[ic] + SPECIFIC_HEAT_AIR*gr[ic]*Tl_ave[ic] + SPECIFIC_HEAT_AIR*gb_h[ic]*T[ic] - L[ic]*gb_v[ic]*Dleaf[ic] 
                  + L[ic]*s[ic]*gb_v[ic]*Told[ic]) / (SPECIFIC_HEAT_AIR*(gr[ic] + gb_h[ic]) + L[ic]*s[ic]*gb_v[ic])
                err = np.nanmax(abs(Tl_wet - Told))
#                Tl_wet = 0.5 * Tl_wet + 0.5 * Told
#                print ('iterNo', iterNo, 'err', err, 'Tl_wet', np.mean(Tl_wet))
                es, s = e_sat(Tl_wet)
                Dleaf = es / P - H2O  #np.maximum(0.0, es / P - H2O)  # [mol/mol]
                s = s / P  # [mol/mol/degC]
                if iterNo == itermax:
                    print 'Maximum number of iterations reached in wet leaf module'
                    print('err', err, 'Tl_wet', np.mean(Tl_wet))
            else:
                err = 0.0

        """ --- energy and water fluxes for wet leaf --- """ ##### or sublimation/deposition ????????? GITHUB SPATHY!!
        # sensible heat flux [W m-2(wet leaf)]
        Hw = SPECIFIC_HEAT_AIR * gb_h * (Tl_wet - T)
        # non-isothermal radiative flux [W m-2 (wet leaf)]???
        Frw = SPECIFIC_HEAT_AIR * gr *(Tl_wet - Tl_ave)
        # evaporation rate from wet leaf [m/s] (negative for condensation)
        Ep = gb_v * Dleaf * MOLAR_MASS_H2O / WATER_DENSITY

#        # Assume no evaporation during rain (CHANGE when energy balance added)
#        if Prec > 0.0:
#            Ep = np.zeros(N)

        """ --- canopy water storage change --- """
        W = self.oldW.copy()  # layerwise canopy storage [m]

        # Unloading in canopy, ensures also that seasonal
        # LAI development does not mess up computations
        for n in reversed(range(0, N)):  # start from highest grid point
            Unload = max(W[n] - Wmax[n], 0.0)  # unloading from layer n
            W[n] -= Unload  # update storage of layer n
            if n != 0:
                W[n-1] += Unload  # unloading added to layer below (if not lower layer)
        # Unload = unloading below canopy [m]

        # timestep subdivision to calculate change in canopy water store, no impact??
        Nsteps = 1  # number of subtimesteps
        subdt = dt / Nsteps  # [s]

        # initiate cumulative variables
        Interc = np.zeros(N)  # interception [m]
        Evap = np.zeros(N)  # evaporation [m]
        Cond = np.zeros(N)  # condesation [m]
        Heat = np.zeros(N)  # sensible heat flux [W m-2(ground)]
        Fr = np.zeros(N)  # sensible heat flux [W m-2(ground)]
        wf = np.zeros(N)  # wetness ratio
        Trfall = 0.0  # throughfall below canopy [m]

        if Prec > 0 or np.any(np.less(Ep, 0)) or np.any(np.greater(W, 0)):
            for t in range(0, Nsteps):
                Ir = np.zeros(N)  # interception rate [m/s]
                dW = np.zeros(N)  # change in storage [m]
                P = np.zeros(N+1)  # precipitation rate to layer [m/s]
                P[-1] = Prec  # above canopy equals precipitation rate [m/s]
                for n in reversed(range(0, N)):  # start from highest grid point
                    if Ep[n] >= 0:  # evaporation case
                        # change in storage [m]
                        dW[n] = (F * P[n+1] / (F * P[n+1] + Ep[n] + eps) * Wmax[n] - W[n]) \
                                * (1.0 - np.exp(-(F * P[n+1] + Ep[n]) * LAIz[n] * subdt / Wmax[n]))
                        # wetness ration in layer
                        if LAIz[n] > 0 and P[n+1] + Ep[n] > 0:
                            wf[n] = (F * P[n+1] - dW[n] / (LAIz[n] * subdt)) / (F * P[n+1] + Ep[n])
                        else:
                            wf[n] = 0.0
                        # interception rate in layer [m/s]
                        Ir[n] = F * (1 - wf[n]) * LAIz[n] * P[n+1]
                        # drainage rate from layer [m/s]
                        P[n] = P[n+1] - Ir[n]
                        # evaporation from layer [m]
                        Evap[n] += wf[n] * LAIz[n] * Ep[n] * subdt
                    else:  # condensation case
                        # change in storage [m]
                        dW[n] = (Wmax[n] - W[n]) \
                                * (1.0 - np.exp(-(F * P[n+1] - Ep[n]) * LAIz[n] * subdt / Wmax[n]))
                        # wetness ration in layer
                        if LAIz[n] > 0 and P[n+1] - Ep[n] > 0:
                            wf[n] = (F * P[n+1] - Ep[n] - dW[n] / (LAIz[n] * subdt)) / (F * P[n+1] - Ep[n])
                        else:
                            wf[n] = 0.0
                        # interception rate in layer [m/s]
                        Ir[n] = F * (1 - wf[n]) * LAIz[n] * P[n+1]
                        # drainage rate from layer [m/s] (incl. increase by condensation drip)
                        P[n] = P[n+1] - Ir[n] - wf[n] * LAIz[n] * Ep[n]
                        # Condensation [m] (incl. condenstation to dry leaf and drip from wet leaf)
                        Cond[n] += LAIz[n] * Ep[n] * subdt
                # Tl_wet represent the whole leaf in case of condensation!!
                # Sensible heat flux [W m-2(ground)] * subdt
                Heat += np.where(Ep >= 0, wf, 1.0) * LAIz * Hw * subdt
                # radiative flux [W m-2(ground)] * subdt
                Fr += np.where(Ep >= 0, wf, 1.0) * LAIz * Frw * subdt
                # update storage [m]
                W += dW
                # interception and throughfall [m]
                Interc += Ir * subdt
                Trfall += P[0] * subdt

        # throughfall to field layer or snowpack
        Trfall = Trfall + Unload
        Trfall_rain = fW[-1] * Trfall  # accoording to above canopy temperature
        Trfall_snow = (1 - fW[-1]) * Trfall

        # H20 source/sink per ground area due to evaporation and condensation [mol m-2 s-1]
        dqsource = (Evap + Cond) / dt / MOLAR_MASS_H2O * WATER_DENSITY
        # heat source [W m-2(ground)]
        dtsource = Heat / dt
        # latent heat flux [W m-2(ground)]
        LE = dqsource * L

        if sum(W) < eps:
            W *= 0.0

        # dry canopy fraction
        df = 1.0 - np.where(Ep >= 0, wf, 1.0)

        # update state variables
        self.W = W
        self.Tl_wet = Tl_wet

        # mass-balance error [m] ! self.W is old storage
        MBE = sum(self.W) - sum(self.oldW) - (Prec * dt - sum(Evap) - sum(Cond) - (Trfall_rain + Trfall_snow))

        return df, Trfall_rain, Trfall_snow, sum(Interc), sum(Evap), sum(Cond), dqsource,  dtsource, Fr/dt, LE, Tl_wet, MBE

    def _update(self):

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