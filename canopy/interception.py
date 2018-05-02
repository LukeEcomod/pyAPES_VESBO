# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 13:21:29 2018

@author: L1656
"""

import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from evapotranspiration import penman_monteith, e_sat, latent_heat
from micromet import leaf_boundary_layer_conductance

#: [kg m-3], Density of water
RHO_WATER = 1000.0
#: [kg mol\ :sup:`-1`\ ], molar mass of H\ :sub:`2`\ O
MOLAR_MASS_H2O = 18.015e-3

class Interception():
    """
    Big-leaf interception model for rain and snowfall
    ****** Kersti! ********
    We need to change this to be multi-layer; see APES Matlab-codes. 
    But for simplicity let's first start with assumption that Tleaf == Tair and neglect
    leaf energy budget.
    """
    def __init__(self, p, LAI, MLinterception, LAIz):
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
        # canopy closure [-]
        self.cf = p['cf']

        # state variables
        self.W = np.minimum(p['w_ini'], p['wmax'] * LAI) # interception storage [m]
        if MLinterception: # compute with multilayer scheme
            self.W = self.W * LAIz

        self._update()

    def _big_leaf(self, dt, LAI, T, Prec, AE, H2O, P, Ra=25.0, U=2.0):
        """
        Args:
            dt: timestep [s]
            LAI: leaf area index [m\ :sup:`2`\ m\ :sup:`-2`\]
            cf: canopy closure [-]
            T: air temperature [degC]
            Prec: precipitation rate [m s\ :sup:`-1`\]
            AE: available energy at canopy level  [W m\ :sup:`-2`\]
            H2O: mixing ratio [mol/mol]
            P: ambient pressure [Pa]
            Ra: canopy aerodynamic resistance [s m\ :sup:`-1`\]
            U: wind speed  at hc [m/s] ----------------------------------------------- ???
        Returns:
            updates sate self.W
            Trfall_rain: throughfall as rainfall [m]
            Trfall_snow: throughfall as snowfall [m]
            Interc: interception [m]
            Evap: evaporation from canopy store [m]
            MBE: mass balance error
        """
        # vapor pressure deficit [Pa]
        esat, _, _ = e_sat(T)  # [Pa]
        VPD = max(0.0, esat - H2O * P)

        Prec = Prec * dt  # [m/s] -> [m]

        """ --- state of precipitation --- """
        # fraction as water [-]
        if T >= self.Tmax:
            fW = 1.0
        elif T <= self.Tmin:
            fW = 0.0
        else:
            fW = (T - self.Tmin) / (self.Tmax - self.Tmin)

        # maximum interception storage capacity [m]
        Wmax = (fW * self.wmax + (1 - fW) * self.wmaxsnow) * LAI

        """ 'potential' evaporation and sublimation rates """
        if (Prec == 0) & (T <= self.Tmin):  # sublimation case
            Ce = 0.01 * ((self.oldW + eps) / Wmax)**(-0.4)  # exposure coeff [-]
            Sh = (1.79 + 3.0 * U**0.5)  # Sherwood numbner [-]
            gi = Sh * self.oldW * 1000 * Ce / 7.68 + eps # [m/s]
            erate = penman_monteith(AE, VPD, T, gi, 1./Ra,  units='m', type='sublimation') * dt
        elif (Prec == 0) & (T > self.Tmin):  # evaporation case
            gs = 1e6
            erate = penman_monteith(AE, VPD, T, gs, 1./Ra, units='m', type='evaporation') * dt
        else:  # negelect evaporation during precipitation events
            erate = 0.0

        """ --- canopy water storage change --- """
        W = self.oldW.copy()  # initiate

        # snow unloading from canopy, ensures also that seasonal
        # LAI development does not mess up computations
        Unload = max(0.0, W - Wmax)
        # update canopy storage [m]
        W = W - Unload

        # interception of rain or snow: asymptotic approach of saturation.
        # Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
        # Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
        Interc = (Wmax - W) * (1.0 - np.exp(-(self.cf / Wmax) * Prec))
        # update canopy storage [m]
        W = W + Interc

        # throughfall to field layer or snowpack
        Trfall = Prec + Unload - Interc
        Trfall_rain = fW * Trfall
        Trfall_snow = (1 - fW) * Trfall

        # evaporate from canopy and update storage [m]
        Evap = np.minimum(erate, W)
        W = W - Evap

        # update state variables
        self.W = W

        # mass-balance error [m] ! self.W is old storage
        MBE = (self.W - self.oldW) - (Prec - Evap - (Trfall_rain + Trfall_snow))

        return Trfall_rain, Trfall_snow, Interc, Evap, MBE

    def _multi_layer(self, dt, lt, LAIz, H2O, U, T, Prec, P):
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
            LAIz (array): leaf area index per canopy layer [m\ :sup:`2`\ m\ :sup:`-2`\]
            H2O (array): mixing ratio for each layer [mol/mol]
            U (array): wind speed for each layer [m s-1]
            T (array): air temperature for each layer [degC]
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

        # Leaf orientation factor with respect to incident Prec; assumed to be 1 when Prec is in vertical
        F = 1.0  # other than 1.0 migth not work!

        """ --- state of precipitation (uses fW[-1] in end of code)--- """
        # fraction as water [-]
        fW = np.ones(N)
        ix = np.where(T < self.Tmin)
        fW[ix] = 0.0
        ix = np.where((T >= self.Tmin) & (T <= self.Tmax))
        fW[ix] = (T[ix] - self.Tmin) / (self.Tmax - self.Tmin)

        # maximum interception storage capacities layerwise [m]
        Wmax = (fW * self.wmax + (1 - fW) * self.wmaxsnow) * LAIz + eps

        """ --- rate of evaporation/condensation --- """ ##### or sublimation/deposition ?????????
        # boundary layer conductances for H2O [mol m-2 s-1]
        _, gb_v, _, _ = leaf_boundary_layer_conductance(U, lt, T, 0.0, P)
        # vapor pressure difference between leaf and air [mol/mol]
        # assumption Tleaf = T
        es, _, _ = e_sat(T)
        Dleaf = np.maximum(0.0, es / P - H2O)  # [mol/mol]
        # latent heat of vaporization at temperature T [J/mol]
        L = latent_heat(T) * MOLAR_MASS_H2O
        # rate as latent heat flux [W m-2] per unit wet/dry leaf area
        LEw = L * gb_v * Dleaf
        # evaporation rate from wet leaf [m/s]
        Ep = gb_v * Dleaf * MOLAR_MASS_H2O / RHO_WATER

        # Assume no evaporation during rain (CHANGE when energy balance added)
        if Prec > 0.0:
            Ep = np.zeros(N)

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
                            wf = (F * P[n+1] - dW[n] / (LAIz[n] * subdt)) / (F * P[n+1] + Ep[n])
                        else:
                            wf = 0.0
                        # interception rate in layer [m/s]
                        Ir[n] = F * (1 - wf) * LAIz[n] * P[n+1]
                        # drainage rate from layer [m/s]
                        P[n] = P[n+1] - Ir[n]
                        # evaporation from layer [m]
                        Evap[n] += wf * LAIz[n] * Ep[n] * subdt
                    else:  # condensation case
                        # change in storage [m]
                        dW[n] = (Wmax[n] - W[n]) \
                                * (1.0 - np.exp(-(F * P[n+1] - Ep[n]) * LAIz[n] * subdt / Wmax[n]))
                        # wetness ration in layer
                        if LAIz[n] > 0 and P[n+1] - Ep[n] > 0:
                            wf = (F * P[n+1] - Ep[n] - dW[n] / (LAIz[n] * subdt)) / (F * P[n+1] - Ep[n])
                        else:
                            wf = 0.0
                        # interception rate in layer [m/s]
                        Ir[n] = F * (1 - wf) * LAIz[n] * P[n+1]
                        # drainage rate from layer [m/s] (incl. increase by condensation drip)
                        P[n] = P[n+1] - Ir[n] - wf[n] * LAIz[n] * Ep[n]
                        # Condensation [m] (incl. condenstation to dry leaf and drip from wet leaf)
                        Cond[n] += LAIz[n] * Ep[n] * subdt
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
        dqsource = (Evap + Cond) / dt / MOLAR_MASS_H2O * RHO_WATER

        if sum(W) < eps:
            W *= 0.0

        # dry canopy fraction
        df = 1 - W / Wmax

        # update state variables
        self.W = W

#        if Prec > 0 or np.any(np.less(Ep, 0)) or np.any(np.greater(W, 0)):
#            print " W " + str(sum(W)) + " Prec " + str(Prec * dt) +" Evap " + str(sum(Evap)) + " Tr " + str((Trfall_rain + Trfall_snow))

        # mass-balance error [m] ! self.W is old storage
        MBE = sum(self.W) - sum(self.oldW) - (Prec * dt - sum(Evap) - sum(Cond) - (Trfall_rain + Trfall_snow))

        return df, Trfall_rain, Trfall_snow, sum(Interc), sum(Evap), dqsource, MBE

    def _update(self):

        self.oldW = self.W.copy()