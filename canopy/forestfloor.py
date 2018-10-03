# -*- coding: utf-8 -*-
"""
.. module: forestfloor
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes forest floor consisting of bryophytes and/or baresoil
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Mar 13 11:59:23 2018
"""
import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from micromet import e_sat, soil_boundary_layer_conductance
from constants import *

class ForestFloor(object):
    r"""Describes forest floor consisting of bryophytes and/or baresoil.
    """
    def __init__(self, p, pp):
        r""" Initializes forestfloor object consisting on bryophytes
        and/or bare soil.
        Args:
            p (dict):
                'mossp' (list):
                    i. (dict): properties of bryotype i
                        'ground_coverage': fraction of moss ground coverage [-]
                        'Wmax'
                        'Wmin'
                        'zr': roughness height [m]
                        'Mdry'
                        'LAI': leaf area index [m2m-2]
                        'Amax': max photo rate [umolm-2s-1]
                        'Q10': temperature sensitivity [-]
                        'R10': base respiration at 10degC [umolm-2s-1]
                        'qeff'
                'soilp' (dict):  --- ALBEDO & ROUGHNESS??
                    'R10': base heterotrophic respiration rate [umolm-2s-1]
                    'Q10': temperature sensitivity [-]
                    'poros': porosity [m3m-3]
                    'limitpara' (list): Skopp respiration function param [a ,b, d, g]
                --- LITTER ???
            pp (dict):
                forestfloor albedo, emissivity, roughness height --- SHOULD CHANGE!
        Returns:
            self (object)
        """
        # Bryotypes
        brtypes = []
        f_bryo = 0.0
        for br_para in p['mossp']:
            brtypes.append(MossLayer(br_para))
            f_bryo += br_para['ground_coverage']
        self.Bryophytes = brtypes
        # soil coverage: baresoil, bryotypes (and litter?)
        self.f_baresoil = 1.0 - f_bryo
        self.f_bryo = f_bryo

        self.R10 = p['soilp']['R10']
        self.Q10 = p['soilp']['Q10']
        self.poros = p['soilp']['poros']
        self.limitpara = p['soilp']['limitpara']

        self.soil_alb = {'Par': pp['soil_Par_alb'],# soil (moss) Par-albedo [-]
                         'Nir': pp['soil_Nir_alb']}  # soil (moss) Nir-albedo [-]
        self.soil_emi = pp['soil_emi']
        self.soil_zr = 0.01  ## INPUT!!!!!!!!

    def run_water_energy_balance(self, dt, Prec, U, T, H2O, P, SWE,
                                  z_can, T_ave, T_soil, h_soil, z_soil, Kh, Kt,
                                  Par_gr, Nir_gr, LWn, Ebal):
        r"""Water and energy balance at forest floor.

        Args:
            dt: timestep [s]
            Prec: throughfall to forest floor [m s-1]
            U: wind speed from first node above ground, corresponds to z_can [m s-1]
            T: air temperature at soil surface [degC]
            H2O: mixing ratio at soil surface [mol mol-1]
            P: ambient pressure [Pa]
            SWE: snow water equivalent [m]
            z_can: heigth of first node above ground [m]
            T_ave: forest floor temperature used in LW computation [degC]
            T_soil: temperature in first soil node [degC]
            h_soil: water potential in first soil node [m]
            z_soil: depth of first soil node [m]
            Kh: hydraulic conductivity at first soil node [m s-1]
            Kt:thermal conductivity at first soil node [W m-1 K-1]
            Par_gr, Nir_gr: incident PAR and NIR at forest floor
            LWn: net longwave radiation at forferst floor
            Ebal (bool): solve energy balance
        Returns:
            Trfall
            Ebryo: evaporation from bryo (mol m-2(ground)s-1) 
            Esoil: soil evaporation (mol m-2(ground)s-1)
            Gsoil: ground heat fluxes (Wm-2)
            LE_gr: latent heat flux (Wm-2)
            Hf: forest floor sensible heat flux (Wm-2)
            Tf: forest floor temperature (degC)
            MBE: water closure
            energy_closure: energy closure
        """
        # initialize fluxes at forest floor
        Trfall = 0.0  # throughfall rate (m s-1)
        Ebryo = 0.0  # evaporation from bryo (mol m-2(ground)s-1) 
        Esoil = 0.0  # soil evaporation (mol m-2(ground)s-1)
        Hf = 0.0  # forest floor sensible heat flux (Wm-2)
        Gsoil = 0.0  # ground heat flux (Wm-2)
        LE_gr = 0.0  # latent heat flux (Wm-2)
        Tf = 0.0  # forest floor temperature (degC)
        Tbryo = 0.0  # bryophyte temperature (degC)
        Fr = 0.0  # radiative flux [W m-2]
        # water (and energy) closure
        MBE = 0.0
        energy_closure = 0.0

        if SWE > 0:  # snow on the ground
            Trfall += Prec
        else:
            if self.f_bryo > 0.0:
                for bryo in self.Bryophytes:
                    ef, trfall, mbe = bryo.waterbalance(dt, Prec * 1e3, U, T, H2O, P=P)
                    Ebryo += bryo.f_cover * ef / MOLAR_MASS_H2O  # mol m-2 s-1
                    LE_gr += bryo.f_cover * ef  / MOLAR_MASS_H2O * LATENT_HEAT
                    Tf += bryo.f_cover * T
                    Trfall += bryo.f_cover * trfall * 1e-3
                    MBE += bryo.f_cover * mbe
            if self.f_baresoil > 0.0:
                Trfall += Prec * self.f_baresoil
                # soil surface energy balance
                T_surf, Hw, Frw, Gw, Ep, LEw, closure = baresoil_energybalance(
                        z_can, U, T, H2O, P, T_ave,
                        soil_alb=self.soil_alb, soil_emi=self.soil_emi, zr=self.soil_zr,
                        T_soil=T_soil, h_soil=h_soil, z_soil=z_soil, Kh=Kh, Kt=Kt,
                        Par_gr=Par_gr, Nir_gr=Nir_gr, LWn=LWn, Ebal=Ebal)
                Hf += self.f_baresoil * Hw
                Gsoil += self.f_baresoil * Gw
                Esoil += self.f_baresoil * Ep
                LE_gr += self.f_baresoil * LEw
                Tf += self.f_baresoil * T_surf
                Fr += self.f_baresoil * Frw
                energy_closure += self.f_baresoil * closure

        return Trfall, Ebryo, Esoil, Gsoil, LE_gr, Hf, Tf, MBE, energy_closure

    def update(self):
        """Updates W of each bryotype to old W
        """
        if self.f_bryo > 0.0:
            for bryo in self.Bryophytes:
                bryo.update()

    def run_CO2(self, dt, Par, T, Ts, Ws, SWE):
        """
        run forest floor model for one timestep
        Args:
            dt - timestep [s]
            Par - incident Par [umolm-2s-1]
            T - air temperture [degC]
            Ts - soil temperature [degC]
            Ws - soil vol. moisture [m3m-3]
            SWE - snow water equivalent, >0 sets E and An to zero
        Returns:
            An - moss net CO2 exchange [umolm-2s-1]
            Rsoil - soil respiration rate [umolm-2s-1]
        """
        An = 0.0
        # moss CO2 exchange when not covered by snow (or yes???????????????)
        if SWE == 0.0:
            if self.f_bryo > 0.0:
                for bryo in self.Bryophytes:
                    An += bryo.f_cover * bryo.co2_exchange(Par, T)

        # soil respiration
        Rsoil, fm = self.soil_respiration(Ts, Ws)

        return An, Rsoil

    def soil_respiration(self, Ts, Wliq):
        """
        computes heterotrophic respiration rate (CO2-flux) based on
        Pumpanen et al. (2003) Soil.Sci.Soc.Am
        Restricts respiration by soil moisuture as in
        Skopp et al. (1990), Soil.Sci.Soc.Am
        Args:
            Ts - soil temperature [degC]
            Wliq - soil vol. moisture content [m3m-3]
        Returns:
            rsoil - soil respiration rate [umolm-2s-1]
            fm - relative modifier (Skopp et al.)
        """
        # Skopp limitparam [a,b,d,g] for two soil types
        # sp = {'Yolo':[3.83, 4.43, 1.25, 0.854], 'Valentine': [1.65,6.15,0.385,1.03]}
        Wliq = np.minimum(self.poros, Wliq)        
        afp = self.poros - Wliq + eps # air filled porosity

        p = self.limitpara

        # unrestricted respiration rate
        rs0 = self.R10 * np.power(self.Q10, (Ts - 10.0) / 10.0)

        # moisture response (substrate diffusion, oxygen limitation)
        fm = np.minimum(p[0]*Wliq**p[2], p[1]*afp**p[3])  # ]0...1]
        fm = np.minimum(fm, 1.0)
        # fm = 1.0
        rsoil = rs0 * fm

        return rsoil, fm

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
        self.Wold = self.W

    def waterbalance(self, dt, Prec, U, T, H2O, P=101300.0):
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
        es, _ = e_sat(T)
        D = np.maximum(0.0, es / P - H2O)  # mol / mol

        # initial water content
        Wo = self.Wold

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

    def update(self):
#        print('Wold',self.Wold, 'W', self.W)
        self.Wold = self.W

def baresoil_energybalance(z_can, U, T, H2O, P, T_ave, soil_alb, soil_emi, zr,
                           T_soil, h_soil, z_soil, Kh, Kt,
                           Par_gr, Nir_gr, LWn, Ebal):
    """
    Solves surface temperature from linearized energy balance
    equation for previously known soil conditions.
    """

    dz_soil = - z_soil
    # initial guess for surface temperature
    T_surf = T
    # boundary layer conductances for H2O and heat [mol m-2 s-1]
    gb_h, _, gb_v = soil_boundary_layer_conductance(u=U, z=z_can, zo=zr, Ta=T, dT=0.0, P=P)  # OK to assume dt = 0.0?
    # radiative conductance [mol m-2 s-1]
    gr = 4.0 * soil_emi * STEFAN_BOLTZMANN * T_ave**3 / SPECIFIC_HEAT_AIR

    # absorbed shortwave radiation
    SW_gr = (1 - soil_alb['Par']) * Par_gr + (1 - soil_alb['Nir']) * Nir_gr

    # Maximum LE
    # atm pressure head in equilibrium with atm. relative humidity
    es_a, _ = e_sat(T)
    RH = H2O * P / es_a  # air relative humidity above ground [-]
    h_atm = GAS_CONSTANT * (DEG_TO_KELVIN + T) * np.log(RH)/(MOLAR_MASS_H2O * GRAVITY)  # [m]
    # maximum latent heat flux constrained by h_atm
    LEmax = -LATENT_HEAT * Kh * (h_atm - h_soil - z_soil) / dz_soil * WATER_DENSITY / MOLAR_MASS_H2O  # [W/m2]

    # LE demand
    # vapor pressure deficit between leaf and air, and slope of vapor pressure curve at T
    es, s = e_sat(T_surf)
    Dsurf = es / P - H2O  # [mol/mol] - allows condensation
    s = s / P  # [mol/mol/degC]
    LE = LATENT_HEAT * gb_v * Dsurf

    if LE > LEmax:
        LE = LEmax
        s = 0.0

    """ --- solve surface temperature --- """
    itermax = 20
    err = 999.0
    iterNo = 0
    while err > 0.01 and iterNo < itermax:
        iterNo += 1
        Told = T_surf
        if Ebal:
            # solve leaf temperature [degC]
            T_surf = (SW_gr + LWn + SPECIFIC_HEAT_AIR*gr*T_ave + SPECIFIC_HEAT_AIR*gb_h*T - LE + LATENT_HEAT*s*gb_v*Told
                      + Kt / dz_soil * T_soil) / (SPECIFIC_HEAT_AIR*(gr + gb_h) + LATENT_HEAT*s*gb_v + Kt / dz_soil)
            err = np.nanmax(abs(T_surf - Told))
#            print ('iterNo', iterNo, 'err', err, 'T_surf', T_surf)
            es, s = e_sat(T_surf)
            Dsurf = es / P - H2O  # [mol/mol] - allows condensation
            s = s / P  # [mol/mol/degC]
            LE = LATENT_HEAT * gb_v * Dsurf
            if LE > LEmax:
                LE = LEmax
                s = 0.0
            if iterNo == itermax:
                print 'Maximum number of iterations reached in surface energy module'
                print('err', err, 'T_surf', T_surf)
        else:
            err = 0.0

    """ --- energy and water fluxes --- """
    # sensible heat flux [W m-2]
    Hw = SPECIFIC_HEAT_AIR * gb_h * (T_surf - T)
    # non-isothermal radiative flux [W m-2]
    Frw = SPECIFIC_HEAT_AIR * gr *(T_surf - T_ave)
    # ground heat flux [W m-2]
    Gw = Kt / dz_soil * (T_surf - T_soil)
    # evaporation rate [mol m-2 s-1]
    Ep = gb_v * Dsurf

    # energy closure
    closure = SW_gr + LWn - Hw - LE - Gw
    
    return T_surf, Hw, Frw, Gw, Ep, LE, closure