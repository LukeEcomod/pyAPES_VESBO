# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 11:59:23 2018

@author: L1656
"""
import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from canopy.evapotranspiration import e_sat

#: [kg mol\ :sup:`-1`\ ], molar mass of H\ :sub:`2`\ O
MOLAR_MASS_H2O = 18.015e-3

class MossLayer():
    def __init__(self, para):
        """
        Moss layer interception, evaporation and CO2 exchange model
        """
        self.f_cover = para['ground_coverage']  # fraction of moss ground coverage [-]
        self.LAI = para['LAI']  # leaf area index
        self.Amax = para['Amax']  # max photo rate [umolm-2s-1]
        self.b = self.Amax / (2.0 * para['qeff'])  # half-saturation par
        self.R10 = para['R10']  # base respiration at 10degC
        self.Q10 = para['Q10']  # temperature sensitivity [-]
        
        self.zr = para['zr']  # roughness height m
        self.Mdry = para['Mdry']
        self.Wmax = para['Mdry']*para['Wmax']
        self.Wmin = para['Mdry']*para['Wmin']

        self.W = para['Wmax']*para['Mdry']      # current water content

    def waterbalance(self, dt, Rn, Prec, U, T, H2O, P=101300.0):
        """
        Moss layer interception, evaporation and water balance.
        Args:
            dt - timestep [s]
            Prec - precipitation [mm]
            U - wind speed [m s-1]
            T - air temperature [degC]
            H2O - mixing ratio [mol mol-1]
            P - ambient pressure [Pa]
        Returns:
            Trfall - trfall rate below moss layer [mm]
            Evap - evaporation rate [mm/s]
            updates self.W
        """
        # VPD at air temperature; neglect condensation conditions
        es, _, _ = e_sat(T)
        D = np.maximum(0.0, es / P - H2O)  # mol / mol

        # initial water content
        Wo = self.W

        # interception and throughfall rate, new storage
        Ir = np.maximum(0.0, np.minimum(Prec, self.Wmax - Wo))
        Trfall = Prec - Ir  # mm

        W = Wo + Ir  # intermediate storage mm

        # evaporation from moss layer: actual conductance is boundary layer x
        # correction for internal resistance
        grel = np.minimum(0.1285 * W / self.Mdry - 0.1285, 1.0)
        gb = grel * self._boundary_layer_conductance(U)

        erate = gb * D  # mol m-2 s-1
        # rate = 1.26*eq_evap(Rn, T, units='mol')  # unrestricted rate
        Evap = np.minimum(erate * MOLAR_MASS_H2O * dt, W - self.Wmin)  # mm
        self.W = W - Evap  # mm

        Mbe = (self.W - Wo) - (Prec - Evap - Trfall)
        # print('Mbe', Mbe)

        return Evap/dt, Trfall, Mbe

    def co2_exchange(self, Par, T):
        """
        moss photosynthetic rate umolm-2s-1
        Args:
            Par (umolm-2s-1)
            T (degC)
        Returns:
            net photosynthetic rate (umolm-2s-1)
        """
        # Williams and Flanagan (1996),Oecologia 108, 38-46. Frolking et al. 1996 GCB
        a = [6.4355, -14.0605, 9.1867, -0.8720]
        b = [-4.3e-5, -8.3e-4, 0.08, 0.1]

        wn = self.W / self.Wmax

        # moisture response, always keep least 5% of capacity
        fW = np.maximum(0.05, a[3] + a[2]*wn + a[1]*wn**2.0 + a[0]*wn**3.0)

        # temperature response
        fT = b[0]*T**3.0 + b[1]*T**2.0 + b[2]*T + b[3]

        # compute photosynthetic rate [umol m-2 s-1]. Slice LAI into 10 layers, attenuate Par
        # exponentially and sum up layerwise photos. rates
        L = np.linspace(0, self.LAI, 10)
        dL = L[1] - L[0]
        Qp = Par*np.exp(-0.7*L)
        Ab = - fW * fT * np.sum(dL * (Qp / (Qp + self.b)))

        del fT, fW

        # respiration rate [umol m-2 s-1]
        if self.W <= 7.0:
            fW = -0.45 + 0.4*self.W - 0.0273*self.W**2.0
        else:
            fW = -0.04*self.W + 1.38

        fW = np.maximum(0.01, np.minimum(1.0, fW))

        Rb = self.R10 * self.Q10 ** ((T - 10.0) / 10.0) * fW

        return Ab + Rb

    def _boundary_layer_conductance(self, U):
        """
        Moss boundary layer conductance as in Rice et al. 2001 eq. 1
        Args:
            zr - roughness lenght scale [m]
            U - mean wind speed [m s-1]
        Returns:
            gb - boundary layer conductance for H2O [mol m-2 s-1]
        """

        Dv = 24e-6  # m2s-1  molecular diffusitity at 20degC
        mu = 15.1e-6  # m2s-1 viscosity of air
        Sc = mu / Dv  # 0.63  # [-] ratio of viscosity to diffusivity
        rhoa = 41.6  # molm-3, density of air

        Re = U*self.zr / mu  # [-], Reynolds number

        gb = rhoa * 10**(-3.18) * Re**1.61 * Dv / self.zr * Sc**(0.33)  # m s-1

        return gb + eps
