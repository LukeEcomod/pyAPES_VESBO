# -*- coding: utf-8 -*-
"""
.. module: canopy
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes H2O, CO2, energy transfer in multilayer canopy.
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Oct 02 09:04:05 2018

References:
Launiainen, S., Katul, G.G., Lauren, A. and Kolari, P., 2015. Coupling boreal
forest CO2, H2O and energy flows by a vertically structured forest canopy – 
Soil model with separate bryophyte layer. Ecological modelling, 312, pp.385-405.
"""

import numpy as np
eps = np.finfo(float).eps
from matplotlib import pyplot as plt
from constants import PAR_TO_UMOL, MOLAR_MASS_H2O, LATENT_HEAT

from radiation import Radiation
from micromet import Micromet
from interception import Interception
from planttype.planttype import PlantType
from forestfloor import ForestFloor
from snow import Snowpack

class CanopyModel(object):
    r""" Represents canopy-soil-atmosphere interactions.

    """
    def __init__(self, cpara, dz_soil):
        r""" Initializes canopy object and submodel objects using given parameters.

        Args:
            cpara (dict):
                'ctr' (dict): switches and specifications for computation
                    'Eflow' (bool): ensemble flow assumption (False solves U_normed on daily timestep)
                    'WMA' (bool): well-mixed assumption (False solves H2O, CO2, T)
                    'StomaModel' (str): stomatal model e.g. 'MEDLYN_FARQUHAR'
                    'Ebal' (bool): True solves energy balance
                    'SwModel' (str): model for shortwave radiation, e.g. 'ZHAOQUALLS'
                    'LwModel' (str): model for longwave radiation, e.g. 'ZHAOQUALLS'
                    'WaterStress' (bool): account for water stress in planttypes --- TRUE NOT SUPPORTED!
                    'seasonal_LAI' (bool): account for seasonal LAI dynamics
                    'pheno_cycle' (bool): account for phenological cycle
                'loc' (dict):
                    'lat': latitude [deg]
                    'lon': longitude [deg]
                'grid' (dict):
                    'zmax': heigth of grid from ground surface [m]
                    'Nlayers': number of layers in grid [-]
                'radi' (dict): radiation model parameters --- ALBEDO SHOULD COME FROM ffloor/planttypes?
                'aero' (dict): micromet model parameters --- ROUGHNESS HEIGHT SHOULD from ffloor?
                'interc_snow' (dict): interception and snow model parameters
                'plant_types' (list):
                    i. (dict): properties of planttype i
                'ffloor': ffloor
                    'mossp' (list):
                        i. (dict): properties of bryotype i
                    'soilp' (dict): parameters to compute soil respiration --- ???
            dz_soil (array): thickness of soilprofile layers, needed for rootzone

        Returns:
            self (object):
                .location (dict):
                    'lat': latitude [deg]
                    'lon': longitude [deg]
                .z (array): canopy model nodes, height from soil surface (= 0.0) [m]
                .dz (float): thickness of canopy layers [m]
                .ones (array): dummy of len(z)
                .Switch_Eflow (bool): ensemble flow assumption
                .Switch_WMA (bool): well-mixed assumption (H2O, CO2, T)
                .Switch_Ebal (bool): solve energy balance
                .LAI (array): total leaf area index [m2 m-2]
                .lad (array): total leaf area density [m2 m-3]
                .hc (float): canopy heigth [m]
                .rad (array): normalized total fine root density distribution [-]
                .Ptypes (list):
                    i. (object): planttype object i
                .Radi_Model (object): radiation model (SW, LW)
                .Micromet_Model (object): micromet model (U, H2O, CO2, T)
                .Interc_Model (object): interception model
                .Snow_Model (object): snow model
                .ForestFloor (object): forest floor object (bryotype/baresoil)
        """

        # --- site location ---
        self.location = cpara['loc']

        # --- grid ---
        self.z = np.linspace(0, cpara['grid']['zmax'], cpara['grid']['Nlayers'])  # grid [m] above ground
        self.dz = self.z[1] - self.z[0]  # gridsize [m]
        self.ones = np.ones(len(self.z))  # dummy

        # --- switches ---
        self.Switch_Eflow = cpara['ctr']['Eflow']
        self.Switch_WMA = cpara['ctr']['WMA']
        self.Switch_Ebal = cpara['ctr']['Ebal']

        # --- Plant types (with phenoligical models) ---
        ptypes = []
        for k in range(len(cpara['plant_types'])):
            p = cpara['plant_types'][k]
            for n in range(len(p['LAImax'])):
                pp = p.copy()
                pp['LAImax'] = p['LAImax'][n]
                pp['lad'] = p['lad'][:, n]
                ptypes.append(PlantType(self.z, pp, dz_soil, ctr=cpara['ctr']))
        self.Ptypes = ptypes

        # --- stand characteristics ---
        self.LAI = sum([pt.LAI for pt in self.Ptypes])  # total leaf area index [m2 m-2]
        self.lad = sum([pt.lad for pt in self.Ptypes])  # total leaf area density [m2 m-3]
        rad = sum([pt.Roots.rad for pt in self.Ptypes])  # total fine root density [m2 m-3]
        self.rad = rad / sum(rad)  # normalized total fine root density distribution [-]
        # canopy height [m]
        if len(np.where(self.lad > 0)[0]) > 0:
            f = np.where(self.lad > 0)[0][-1]
            self.hc = self.z[f].copy()
        else:
            self.hc = 0.0

        """ initialize submodel instance """
        # radiation
        self.Radi_Model = Radiation(cpara['radi'], cpara['ctr'])

        # micromet
        self.Micromet_Model = Micromet(self.z, self.lad, self.hc, cpara['aero'])

        # interception
        self.Interc_Model = Interception(cpara['interc_snow'], self.lad * self.dz)

        # snow
        self.Snow_Model = Snowpack(cpara['interc_snow'])

        # forest floor
        self.ForestFloor = ForestFloor(cpara['ffloor'], cpara['radi'])

    def run_daily(self, doy, Ta, PsiL=0.0):
        r""" Computatations at daily timestep.
        Updates planttypes and total canopy leaf area index and phenological state. 
        Recomputes normalize flow statistics with new leaf area density profile.

        Args:
            doy (float): day of year [days]
            Ta (float): mean daily air temperature [degC]
            PsiL (float): leaf water potential [MPa] --- CHECK??
        """

        """ update physiology and leaf area of planttypes and canopy"""
        for pt in self.Ptypes:
            pt.update_daily(doy, Ta, PsiL=PsiL)  # updates pt properties
        self.LAI = sum([pt.LAI for pt in self.Ptypes])  # total leaf area index [m2 m-2]
        self.lad = sum([pt.lad for pt in self.Ptypes])  # total leaf area density [m2 m-3]

        """ normalized flow statistics in canopy with new lad """
        if self.Switch_Eflow and self.Ptypes[0].Switch_lai:
            self.Micromet_Model.normalized_flow_stats(self.z, self.lad, self.hc)

    def run_timestep(self, dt, forcing, Tsoil, Wsoil, soil_state):
        r""" Calculates one timestep and updates state of CanopyModel object.

        Args:
            dt: timestep [s]
            forcing (dataframe): meteorological forcing data
                'Prec': precipitation rate [m s-1]
                'Tair': air temperature [degC]
                'dirPar': direct fotosynthetically active radiation [W m-2]
                'diffPar': diffuse fotosynthetically active radiation [W m-2]
                'dirNir': direct near infrared radiation [W m-2]
                'diffNir': diffuse near infrare active radiation [W m-2]
                'LWin': Downwelling long wave radiation [W/m2]
                'U': wind speed [m s-1]
                'Ustar': friction velocity [m s-1]
                'CO2': ambient CO2 mixing ratio [ppm]
                'H2O': ambient H2O mixing ratio [mol/mol]
                'P': pressure [Pa]
                'Zen': solar zenith angle [rad]
            Tsoil, Wsoil: used for soils respiration --- GOING TO CHANGE?
            soil_state (dict): properties of first soil node
                'depth': [m]
                'temperature': [\ :math:`^{\circ}`\ C]
                'water_potential': [m]
                'hydraulic_conductivity': [m s\ :sup:`-1`\ ]
                'thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ]

        Returns:
            fluxes (dict)
            states (dict)
        """

        P = forcing['P']

        """ --- flow stats --- """
        if self.Switch_Eflow is False:
            # recompute normalized flow statistics in canopy with current meteo
            self.Micromet_Model.normalized_flow_stats(
                    z=self.z,
                    lad=self.lad,
                    hc=self.hc,
                    Utop=forcing['U']/(forcing['Ustar'] + eps))
        # update U
        U = self.Micromet_Model.update_state(ustaro=forcing['Ustar'])

        """ --- SW profiles within canopy --- """
        # --- PAR ---
        Q_sl1, Q_sh1, q_sl1, q_sh1, q_soil1, f_sl, Par_gr = self.Radi_Model.SW_profiles(
                LAIz=self.lad * self.dz,
                zen=forcing['Zen'],
                dirSW=forcing['dirPar'],
                diffSW=forcing['diffPar'],
                radtype='Par')

        if self.Switch_Ebal:
            # --- NIR ---
            Q_sl2, Q_sh2, q_sl2, q_sh2, q_soil2, _, Nir_gr = self.Radi_Model.SW_profiles(
                LAIz=self.lad * self.dz,
                zen=forcing['Zen'],
                dirSW=forcing['dirNir'],
                diffSW=forcing['diffNir'],
                radtype='Nir')
            # absorbed radiation by sunlit and shaded leafs [W m-2(leaf)]
            SWabs_sl = q_sl1 + q_sl2
            SWabs_sh = q_sh1 + q_sh2
        else: # values need to be given but are not used
            SWabs_sl = 0.0
            SWabs_sh = 0.0
            Nir_gr = 0.0

        """ --- iterative solution of H2O, CO2, T, Tleaf and Tsurf --- """
        max_err = 0.01  # maximum relative error
        max_iter = 20  # maximum iterations
        gam = 0.5  # weight for old value in iterations
        err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
        Switch_WMA = self.Switch_WMA

        # initialize canopy micrometeo, Tleaf and Tsurf
        T, H2O, CO2, Tleaf, Tsurf = self._init_state(forcing)

        iter_no = 0

        while (err_t > max_err or err_h2o > max_err or
               err_co2 > max_err or err_Tl > max_err or
               err_Ts > max_err) and iter_no <= max_iter:

            iter_no += 1
            Tleaf_prev = Tleaf.copy()
            Tsurf_prev = Tsurf

            if self.Switch_Ebal:
                """ --- compute LW profiles and net LW at leaf --- """
                LWl, LWd, LWu, gr = self.Radi_Model.LW_profiles(
                        LAIz=self.lad*self.dz,
                        Tleaf=Tleaf,
                        Tsurf=Tsurf,
                        LWdn0=forcing['LWin'])
            else: # values need to be given but are not used
                LWl = 0.0 * self.ones
                gr = 0.0
                LWd, LWu = np.array([0.0]), np.array([0.0])

            # --- T, h2o and co2 sink/source terms
            tsource = 0.0 * self.ones
            qsource = 0.0 * self.ones
            csource = 0.0 * self.ones
            frsource = 0.0 * self.ones
            lsource = 0.0 * self.ones
            # dark respiration
            Rstand = 0.0
            # create empty list to append planttype results
            pt_stats = []

            """ --- interception and interception evaporation --- """
            # isothermal radiation balance at each layer,
            # sunlit and shaded leaves together [W m-2]
            Rabs = SWabs_sl*f_sl + SWabs_sh*(1 - f_sl) + LWl  # zero if Ebal not solve
            df, Trfall_rain, Trfall_snow, Interc, Evap, Cond, Ew, Hw, dFr, LEw, Tleaf_w, MBE_interc = self.Interc_Model.run(
                    dt=dt,
                    lt=0.1,  ############ to inputs!!!
                    LAIz=self.lad * self.dz,
                    H2O=H2O,
                    U=U,
                    T=T,
                    Rabs=Rabs,
                    Prec=forcing['Prec'],
                    P=P,
                    Ebal=self.Switch_Ebal,
                    Tl_ave=Tleaf_prev,
                    gr=gr)
            # --- update ---
            # heat and h2o sink/source terms
            tsource += Hw / self.dz  # [W m-3]
            qsource += Ew / self.dz  # [mol m-3 s-1]
            lsource += LEw / self.dz
            frsource += dFr / self.dz  # [W m-3]
            # canopy leaf temperature
            Tleaf = Tleaf_w * (1 - df) * self.lad

            """ --- leaf gas-exchange at each layer and for each PlantType --- """
            for pt in self.Ptypes:
                # --- sunlit and shaded leaves
                pt_stats_i, dtsource, dqsource, dcsource, dRstand, dFr = pt.leaf_gasexchange(
                        f_sl, H2O, CO2, T, U, P, Q_sl1*PAR_TO_UMOL, Q_sh1*PAR_TO_UMOL,
                        SWabs_sl=SWabs_sl, SWabs_sh=SWabs_sh, LWl=LWl, # only used if Ebal=True
                        df=df, Ebal=self.Switch_Ebal, Tl_ave=Tleaf_prev, gr=gr) # only used if Ebal=True
                # --- update ---
                # heat, h2o and co2 sink/source terms
                tsource += dtsource  # W m-3
                qsource += dqsource  # mol m-3 s-1
                csource += dcsource  # umol m-3 s-1
                lsource += dqsource * LATENT_HEAT
                frsource += dFr  # [W m-3]
                # dark respiration umol m-2 s-1
                Rstand +=  dRstand
                # PlantType results
                pt_stats.append(pt_stats_i)
                # canopy leaf temperature
                Tleaf += pt_stats_i['Tleaf'] * df * pt.lad
            # canopy leaf temperature as weighted average
            Tleaf = Tleaf / (self.lad + eps)

            err_Tl = max(abs(Tleaf - Tleaf_prev))

            """ --- snowpack ---"""
            PotInf, MBE_snow = self.Snow_Model._run(
                    dt=dt,
                    T=T[0],
                    Trfall_rain=Trfall_rain,
                    Trfall_snow=Trfall_snow)

            """ --- forest floor --- """
            # Water balance at forest floor (interception and evaporation of moss layer)
            Trfall_gr, Ebryo, Esoil, Gsoil, LE_gr, H_gr, Tsurf, MBE_ff, ene_ff = self.ForestFloor.run_water_energy_balance(
                    dt=dt,
                    Prec=PotInf,
                    U=U[1],  # from first node above ground, corresponds to z_can
                    T=T[0],
                    H2O=H2O[0],
                    P=P,
                    SWE=self.Snow_Model.swe,
                    z_can=self.z[1], T_ave=Tsurf_prev, 
                    T_soil=soil_state['temperature'],
                    h_soil=soil_state['water_potential'],
                    z_soil=soil_state['depth'],
                    Kh=soil_state['hydraulic_conductivity'],
                    Kt=soil_state['thermal_conductivity'],
                    Par_gr=Par_gr, Nir_gr=Nir_gr, LWn=LWd[0]-LWu[0], Ebal=self.Switch_Ebal)
            err_Ts = abs(Tsurf - Tsurf_prev)
            # CO2 flux at forest floor
            An_gr, R_gr = self.ForestFloor.run_CO2(
                    dt=dt,
                    Par=Par_gr * PAR_TO_UMOL,
                    T=T[1],  ## 1 or 0 ??
                    Ts=Tsoil,
                    Ws=Wsoil,
                    SWE=self.Snow_Model.swe)
            Fc_gr = An_gr + R_gr

            """  --- solve scalar profiles (H2O, CO2, T) """
            if Switch_WMA is False:
                H2O, CO2, T, err_h2o, err_co2, err_t = self.Micromet_Model.scalar_profiles(
                        gam, H2O, CO2, T, P,
                        source={'H2O': qsource, 'CO2': csource, 'T': tsource},
                        lbc={'H2O': Ebryo + Esoil, 'CO2': Fc_gr, 'T': H_gr})
                if (iter_no > max_iter or any(np.isnan(T)) or err_t > 50.0 or
                    any(np.isnan(H2O)) or any(np.isnan(CO2))):
                    Switch_WMA = True  # if no convergence, re-compute with WMA -assumption
                    # reset values
                    iter_no = 0
                    err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
                    T, H2O, CO2, Tleaf, Tsurf = self._init_state(forcing)
                    self.Interc_Model.Tl_wet = None
            else:
                err_h2o, err_co2, err_t = 0.0, 0.0, 0.0

        """ --- update state variables --- """
        self.Snow_Model.update()  # snow water equivalent SWE, SWEi, SWEl
        self.Interc_Model.update()  # interception storage
        self.ForestFloor.update()  # bryo storage

        # ecosystem fluxes
        flux_co2 = np.cumsum(csource)*self.dz + Fc_gr  # umolm-2s-1
        flux_latent_heat = np.cumsum(lsource)*self.dz + LE_gr  # Wm-2
        flux_sensible_heat = np.cumsum(tsource)*self.dz + H_gr  # Wm-2

        # net ecosystem exchange umolm-2s-1
        NEE = flux_co2[-1]
        # ecosystem respiration umolm-2s-1
        Reco = Rstand + R_gr
        # ecosystem GPP umolm-2s-1
        GPP = - NEE + Reco
        # stand transpiration [m/dt]
        Tr = sum([pt_st['E'] * MOLAR_MASS_H2O * 1e-3 * dt for pt_st in pt_stats])

        # energy closure of canopy  -- THIS IS EQUAL TO frsource (the error caused by linearizing sigma*ef*T^4)
        energy_closure =  sum((SWabs_sl*f_sl + SWabs_sh*(1 - f_sl) + LWl) * self.lad * self.dz) - (  # absorbed radiation
                          sum(tsource*self.dz) + sum(lsource*self.dz))  # sensible and latent heat flux

        # evaporation from moss layer [m/dt]
        Efloor = Ebryo * MOLAR_MASS_H2O * dt * 1e-3
        Esoil = Esoil * MOLAR_MASS_H2O * dt * 1e-3

        # return state and fluxes in dictionary
        state = {'interception_storage': sum(self.Interc_Model.W),
                 'snow_water_equivalent': self.Snow_Model.swe,
                 'LAI': self.LAI,
                 'lad': self.lad,
                 'sunlit_fraction': f_sl,
                 'phenostate': sum([pt.LAI * pt.pheno_state for pt in self.Ptypes])/self.LAI,
                 'IterWMA': iter_no
                 }
        fluxes = {'potential_infiltration': Trfall_gr / dt,
                  'throughfall': (Trfall_rain + Trfall_snow) / dt,
                  'interception': Interc / dt,
                  'evaporation': Evap / dt,
                  'condensation': Cond / dt,
                  'transpiration': Tr / dt,
                  'moss_evaporation': Efloor / dt,
                  'baresoil_evaporation': Esoil / dt,
                  'NEE': NEE,
                  'GPP': GPP,
                  'respiration': Reco,
                  'LE': flux_latent_heat[-1],
                  'co2_flux': flux_co2,  # [umol m-2 s-1]
                  'latent_heat_flux': flux_latent_heat,  # [W m-2]
                  'pt_transpiration': np.array([pt_st['E'] * MOLAR_MASS_H2O * 1e-3 for pt_st in pt_stats]),
                  'pt_gpp': np.array([pt_st['An'] + pt_st['Rd'] for pt_st in pt_stats]),
                  'pt_respiration': np.array([pt_st['Rd'] for pt_st in pt_stats]),
                  'water_closure_canopy': MBE_interc,
                  'water_closure_snow': MBE_snow,
                  'water_closure_forestfloor': MBE_ff
                  }

        if self.Switch_WMA is False:
            state.update({'h2o': H2O,
                          'co2': CO2,
                          'temperature': T,
                          'WMA_assumption': 1.0*Switch_WMA})

        if self.Switch_Ebal:
            state.update({'Tleaf_wet': np.where(self.lad > 0.0, Tleaf_w, np.nan),
                          'Tleaf_sl': np.where(self.lad > 0.0, sum([pt_st['Tleaf_sl'] for pt_st in pt_stats])/(self.lad+eps), np.nan),
                          'Tleaf_sh':np.where(self.lad > 0.0, sum([pt_st['Tleaf_sh'] for pt_st in pt_stats])/(self.lad+eps), np.nan),
                          'Tleaf': np.where(self.lad > 0.0, Tleaf, np.nan),
                          'Tsurf': Tsurf})
            fluxes.update({'net_leaf_radiation': Rabs,
                           'net_leaf_LW': LWl,
                           'ground_heat_flux': Gsoil,
                           'soil_sensible_heat_flux': H_gr,
                           'sensible_heat_flux': flux_sensible_heat,  # [W m-2]
                           'energy_closure_canopy': energy_closure,
                           'fr_source': sum(frsource*self.dz)})

        return fluxes, state

    def _init_state(self, forcing):
        """ initialize state variables """

        T = self.ones * ([forcing['Tair']])
        H2O = self.ones * ([forcing['H2O']])
        CO2 = self.ones * ([forcing['CO2']])
        Tleaf = T.copy() * self.lad / (self.lad + eps)
        Tsurf = T[0].copy()

        return T, H2O, CO2, Tleaf, Tsurf

