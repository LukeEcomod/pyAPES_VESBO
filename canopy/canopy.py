# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:01:50 2017

@author: slauniai

******************************************************************************
CanopyModel:

Gridded canopy and snow hydrology model for SpatHy -integration
Uses simple schemes for computing water flows and storages within vegetation
canopy and snowpack at daily or sub-daily timesteps.

(C) Samuli Launiainen, 2016-2017
last edit: 1.11.2017: Added CO2-response to dry_canopy_et
******************************************************************************

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

class CanopyModel():
    def __init__(self, cpara, dz_soil):
        """
        initializes CanopyModel -object

        Args:
            cpara (dict): see canopy_parameters.py
        Returns:
            self (object):
                parameters/state variables:
                    Switch_MLM (boolean): control for multilayer computation
                    location (dict): 
                        'lat'(float): latitude
                        'lon'(float): longitude
                    pheno_state (float): LAI-weigthed average phenological state of canopy (0..1)
                grid (only for multilayer computation):
                    z, dz, Nlayers (floats)
                objects:
                    Ptypes (.Pheno_Model, .LAI_model), Forestfloor (.Evap_Model, ...)
                    Radi_Model, Micromet_Model, Interc_model, Canopy_Tr, Snow_Model
        """

        # --- site location ---
        self.location = cpara['loc']

        # --- grid ---
        self.Nlayers = cpara['grid']['Nlayers']  # number of layers
        self.z = np.linspace(0, cpara['grid']['zmax'], self.Nlayers)  # grid [m] above ground
        self.dz = self.z[1] - self.z[0]  # gridsize [m]

        # --- switches ---
        self.Switch_Eflow = cpara['ctr']['Eflow']
        self.Switch_WMA = cpara['ctr']['WMA']
        self.Switch_Ebal = cpara['ctr']['Ebal']

        self.ones = np.ones(len(self.z))  # dummy

        # --- Plant types (with phenoligical models) ---
        ptypes = []
        stand_lai, stand_lad = 0.0, 0.0
        # para = [pinep, shrubp]
        for k in range(len(cpara['plant_types'])):
            p = cpara['plant_types'][k]
            for n in range(len(p['LAImax'])):
                pp = p.copy()
                pp['LAImax'] = p['LAImax'][n]
                pp['lad'] = p['lad'][:, n]
                ptypes.append(PlantType(self.z, pp, dz_soil,ctr=cpara['ctr']))
                stand_lai += pp['LAImax']
                stand_lad += pp['lad']
        self.Ptypes = ptypes
        del p, pp, k, ptypes

        # --- stand characteristics ---
        self.LAI = stand_lai  # canopy total 1-sided leaf area index [m2m-2]
        self.lad = stand_lad  # canopy total 1-sided leaf area density [m2m-3]
        rad = sum([pt.Roots.rad for pt in self.Ptypes])
        self.rad = rad / sum(rad)
        # canopy height [m]
        if len(np.where(self.lad > 0)[0]) > 0:
            f = np.where(self.lad > 0)[0][-1]
            self.hc = self.z[f].copy()
        else:
            self.hc = 0.0

        """ initialize submodels """
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

    def _run_daily(self, doy, Ta, PsiL=0.0):  ### PsiL???
        """
        Computatations done at daily timestep. Updates 
        LAI and phenology
        Args:
            doy: day of year [days]
            Ta: mean daily air temperature [degC]
            PsiL: leaf water potential [MPa]
        Returns:
            None
        """

        """ update physiology and leaf area of planttypes & canopy"""
        stand_lai, stand_lad, pheno_state = 0.0, 0.0, 0.0
        for pt in self.Ptypes:
            pt._update_daily(doy, Ta, PsiL=PsiL)  # updates pt properties
            stand_lad += pt.lad
            stand_lai += pt.LAI
            pheno_state += pt.LAI * pt.pheno_state
        self.lad = stand_lad
        self.LAI = stand_lai
        self.pheno_state = pheno_state / stand_lai

        """ normalized flow statistics in canopy with new lad """
        if self.Switch_Eflow and self.Ptypes[0].Switch_lai:
            self.Micromet_Model.normalized_flow_stats(self.z, self.lad, self.hc)

    def _run_timestep(self, dt, forcing, Rew, beta, Tsoil, Wsoil,
                      Ts, hs, zs, Kh, Kt):  # properties of first soil node
        """
        Runs CanopyModel instance for one timestep
        Args:
            dt: timestep [s]
            forcing (dataframe): meteorological forcing data
                'Prec': precipitation rate [m s-1]
                'Tair': air temperature [degC]
                'Rg': global radiation [W m-2]
                'vpd': vapor pressure deficit [kPa]
                'Par': fotosynthetically active radiation [W m-2]
                'U': wind speed [m s-1]
                'CO2': atm. CO2 mixing ratio [ppm]
                'P': pressure [Pa]
            Rew - relative extractable water [-], scalar or matrix
            beta - term for soil evaporation resistance (Wliq/FC) [-]
        Returns:
            fluxes (dict):
            state (dict):
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

        """ start iterative solution of scalar profiles, Tleaf and Tsurf """
        max_err = 0.01  # maximum relative error
        max_iter = 20  # maximum iterations
        gam = 0.5  # weight for old value in iterations
        err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
        Switch_WMA = self.Switch_WMA

        # initialize canopy micrometeo, Tleaf and Tsurf
        T = self.ones * ([forcing['Tair']])
        H2O = self.ones * ([forcing['H2O']])
        CO2 = self.ones * ([forcing['CO2']])
        Tleaf = T.copy() * self.lad / (self.lad + eps)
        Tsurf = T[0].copy()

        iter_no = 0

        while (err_t > max_err or
               err_h2o > max_err or
               err_co2 > max_err or
               err_Tl > max_err or
               err_Ts > max_err) and iter_no <= max_iter:

            iter_no += 1
            Tleaf_prev = Tleaf.copy()
            Tsurf_prev = Tsurf

            if self.Switch_Ebal:
                """ --- compute LW profiles and net isothermal LW at leaf --- """
                LWl, LWd, LWu, gr = self.Radi_Model.LW_profiles(
                        LAIz=self.lad*self.dz,
                        Tleaf=Tleaf,
                        Tsurf=Tsurf,
                        LWdn0=forcing['LWin'])
            else: # values need to be given but are not used
                LWl = np.zeros(self.Nlayers)
                gr = 0.0
                LWd, LWu = np.array([0.0]),np.array([0.0])

            # --- T, h2o and co2 sink/source terms
            tsource = np.zeros(self.Nlayers)
            qsource = np.zeros(self.Nlayers)
            csource = np.zeros(self.Nlayers)
            frsource = np.zeros(self.Nlayers)
            lsource = np.zeros(self.Nlayers)
            # dark respiration
            Rstand = 0.0
            # PlantType results, create empty list to append results
            pt_stats = []

            """ --- multilayer interception and interception evaporation --- """
            # isothermal radiation balance at each layer,
            # sunlit and shaded leaves together [W m-2]
            Rabs = SWabs_sl*f_sl + SWabs_sh*(1 - f_sl) + LWl  # zero if Ebal not solve
            df, Trfall_rain, Trfall_snow, Interc, Evap, Cond, Ew, Hw, dFr, LEw, Tleaf_w, MBE_interc = self.Interc_Model._multi_layer(
                    dt=dt,
                    lt=0.1,  ############ to inputs!!!
                    ef=self.Radi_Model.leaf_emi,
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
            Tleafff = Tleaf_w * (1 - df) * self.lad

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
                Tleafff += pt_stats_i['Tleaf'] * df * pt.lad
            # canopy leaf temperature as weighted average
            Tleafff = Tleafff / (self.lad + eps)
            Tleaf = Tleafff.copy()

            err_Tl = max(abs(Tleaf - Tleaf_prev))

            """ --- snowpack ---"""
            PotInf, MBE_snow = self.Snow_Model._run(
                    dt=dt,
                    T=T[0],
                    Trfall_rain=Trfall_rain,
                    Trfall_snow=Trfall_snow)

            """ --- forest floor --- """
            # Water balance at forest floor (interception and evaporation of moss layer)
            Trfall_gr, Ebryo, Esoil, Gsoil, LE_gr, H_gr, Tsurf, MBE_ff, ene_ff = self.ForestFloor._run_water_energy_balance(
                    dt=dt,
                    Prec=PotInf,
                    U=U[1],  # from first node above ground, corresponds to z_can
                    T=T[0],
                    H2O=H2O[0],
                    P=P,
                    SWE=self.Snow_Model.swe,
                    z_can=self.z[1], T_ave=Tsurf_prev, 
                    T_soil=Ts, h_soil=hs, z_soil=zs, Kh=Kh, Kt=Kt,  # from soil model
                    Par_gr=Par_gr, Nir_gr=Nir_gr, LWn=LWd[0]-LWu[0], Ebal=self.Switch_Ebal)
            err_Ts = abs(Tsurf - Tsurf_prev)
            # CO2 flux at forest floor
            An_gr, R_gr = self.ForestFloor._run_CO2(
                    dt=dt,
                    Par=Par_gr*PAR_TO_UMOL,
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
#                print('iterNo', iter_no, 'err_h2o', err_h2o, 'err_co2', err_co2, 'err_t', err_t)
                if (iter_no > max_iter or any(np.isnan(T)) or err_t > 50.0 or
                    any(np.isnan(H2O)) or any(np.isnan(CO2))):  # if no convergence, re-compute with WMA -assumption
                    Switch_WMA = True
#                    print 'Maximum number of iterations reached, WMA assumed'
#                    print('iter_no', iter_no, 'err_h2o', err_h2o, 'err_co2', err_co2, 'err_t', err_t)
                    iter_no = 0
                    err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
                    T = self.ones * ([forcing['Tair']])
                    H2O = self.ones * ([forcing['H2O']])
                    CO2 = self.ones * ([forcing['CO2']])
                    Tleaf = T.copy() * self.lad / (self.lad + eps)
                    self.Interc_Model.Tl_wet = None
                    Tsurf = T[0].copy()
            else:
                err_h2o, err_co2, err_t = 0.0, 0.0, 0.0

#        print('iterNo', iter_no, 'err_h2o', err_h2o, 'err_co2', err_co2, 'err_t', err_t, 'err_Tl', err_Tl)

        """ --- update state variables --- """
        self.Snow_Model._update()  # snow water equivalent SWE, SWEi, SWEl
        self.Interc_Model._update()  # interception storage
        self.ForestFloor._update()  # bryo storage

        # ecosystem fluxes
        Fc = np.cumsum(csource)*self.dz + Fc_gr  # umolm-2s-1
        LE = np.cumsum(lsource)*self.dz + LE_gr  # Wm-2

        # net ecosystem exchange umolm-2s-1
        NEE = Fc[-1]
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
        PotInf = Trfall_gr

        # return state and fluxes in dictionary
        state = {'snow_water_equivalent': self.Snow_Model.swe,
                 'LAI': self.LAI,
                 'phenostate': self.pheno_state
                 }
        fluxes = {'potential_infiltration': PotInf / dt,
                  'throughfall': (Trfall_rain + Trfall_snow) / dt,
                  'interception': Interc / dt,
                  'evaporation': Evap / dt,
                  'condensation': Cond / dt,
                  'transpiration': Tr / dt,
                  'moss_evaporation': Efloor / dt,
                  'baresoil_evaporation': Esoil / dt,
                  'ground_heat_flux': Gsoil,
                  'U_ground': U[0],
                  'MBE1': MBE_interc,
                  'MBE2': MBE_snow,
                  'MBE3': MBE_ff
                  }
        fluxes.update({'NEE': NEE,
                       'GPP': GPP,
                       'Reco': Reco,
                       'Rsoil': R_gr,
                       'LE': LE[-1],
                       'LEgr': LE_gr,
                       'sensible_heat_flux':H_gr,
                       'GPPgr': -An_gr,
                       'pt_transpiration': np.array([pt_st['E'] * MOLAR_MASS_H2O * 1e-3 for pt_st in pt_stats]),
                       'pt_An': np.array([pt_st['An']+pt_st['Rd'] for pt_st in pt_stats]),
                       'pt_Rd': np.array([pt_st['Rd'] for pt_st in pt_stats]),
                       'Tleaf': np.where(self.lad > 0.0, Tleafff, np.nan),
                       'Tsurf': Tsurf,
                       'Tleaf_wet':np.where(self.lad > 0.0, Tleaf_w, np.nan),
                       'Tleaf_sl': np.where(self.lad > 0.0, sum([pt_st['Tleaf_sl'] for pt_st in pt_stats])/(self.lad+eps), np.nan),
                       'Tleaf_sh':np.where(self.lad > 0.0, sum([pt_st['Tleaf_sh'] for pt_st in pt_stats])/(self.lad+eps), np.nan),
                       'Rabs': Rabs,
                       'LWleaf': LWl,
                       'IterWMA': iter_no,
                       'energy_closure': energy_closure,
                       'Frsource': sum(frsource*self.dz)})
        state.update({'wind_speed': U,
                      'PAR_sunlit': Q_sl1,
                      'PAR_shaded': Q_sh1,
                      'lad': self.lad,
                      'sunlit_fraction': f_sl})
        if self.Switch_WMA is False:
            state.update({'h2o': H2O,
                          'co2': CO2,
                          'T': T,
                          'WMA': 1.0*Switch_WMA})
        state.update({'interception_storage': sum(self.Interc_Model.W),
                      'Tleaf_wet': Tleaf_w})

        return fluxes, state

