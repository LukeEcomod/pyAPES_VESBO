#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module: canopy
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes H2O, CO2, energy transfer in multilayer canopy.
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Oct 02 09:04:05 2018

Note:
    Migrated to python3:
        - absolute imports
        - for each loop: usage of enumerate() instead of list(range(len(foo)))
        - if dictionary IS modified in a for-loop dict.items()/.keys()/.values() wrapped in a list
        (this is not necessary if only modifying and not adding/deleting)
        - if dictionary IS NOT modified in a for-loop, wrapped anyways (for caution).

References:
Launiainen, S., Katul, G.G., Lauren, A. and Kolari, P., 2015. Coupling boreal
forest CO2, H2O and energy flows by a vertically structured forest canopy – 
Soil model with separate bryophyte layer. Ecological modelling, 312, pp.385-405.
"""

import numpy as np
from copy import deepcopy
from matplotlib import pyplot as plt
from .constants import MOLAR_MASS_H2O, EPS

from .radiation import Radiation
from .micromet import Micromet
from .interception import Interception
from .planttype.planttype import PlantType
from .forestfloor.forestfloor import ForestFloor

import logging
logger = logging.getLogger(__name__)

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
                'grid' (dict):
                    'zmax': heigth of grid from ground surface [m]
                    'Nlayers': number of layers in grid [-]
                'radiation' (dict): radiation model parameters
                'micromet' (dict): micromet model parameters --- ROUGHNESS HEIGHT SHOULD from ffloor?
                'interception' (dict): interception and snow model parameters
                'planttypes' (list):
                    i. (dict): properties of planttype i
                'forestfloor': forestfloor parameters
                    'bryophytes'(list):
                        i. (dict): properties of bryotype i
                    'baresoil' (dict): baresoil parameters
                    'snowpack' (dict): smow model parameters
                    'initial_conditions': initial conditions for forest floor
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
                .planttypes (list):
                    i. (object): planttype object i
                .radiation (object): radiation model (SW, LW)
                .micromet (object): micromet model (U, H2O, CO2, T)
                .intercetion (object): interception model
                .forestfloor (object): forest floor object (bryotype/baresoil/snow)
        """

#         --- site location ---
#        self.location = cpara['loc']

        # --- grid ---
        self.z = np.linspace(0, cpara['grid']['zmax'], cpara['grid']['Nlayers'])  # grid [m] above ground
        self.dz = self.z[1] - self.z[0]  # gridsize [m]
        self.ones = np.ones(len(self.z))  # dummy

        # --- switches ---
        self.Switch_Eflow = cpara['ctr']['Eflow']
        self.Switch_WMA = cpara['ctr']['WMA']
        self.Switch_Ebal = cpara['ctr']['Ebal']

        #logger.info('Eflow: %s, WMA: %s, Ebal: %s',
        #            self.Switch_Eflow,
        #            self.Switch_WMA,
        #            self.Switch_Ebal)

        # --- Plant types (with phenoligical models) ---
        ptypes = []

        # Dictionary IS NOT modified in loop. Wrapped in to a list just in case.
        for p in cpara['planttypes'].values():

            for idx, lai_max in enumerate(p['LAImax']):

                pp = p.copy()
                pp['LAImax'] = lai_max
                pp['lad'] = p['lad'][:, idx]
                ptypes.append(PlantType(self.z, pp, dz_soil, ctr=cpara['ctr']))

        self.planttypes = ptypes

        # --- stand characteristics ---
        self.LAI = sum([pt.LAI for pt in self.planttypes])  # total leaf area index [m2 m-2]
        self.lad = sum([pt.lad for pt in self.planttypes])  # total leaf area density [m2 m-3]
        rad = sum([pt.Roots.rad for pt in self.planttypes])  # total fine root density [m2 m-3]
        self.rad = rad / sum(rad)  # normalized total fine root density distribution [-]
        # canopy height [m]
        if len(np.where(self.lad > 0)[0]) > 0:
            f = np.where(self.lad > 0)[0][-1]
            self.hc = self.z[f].copy()
        else:
            self.hc = 0.0

        # --- radiation, micromet, interception, and forestfloor instances
        self.radiation = Radiation(cpara['radiation'], cpara['ctr'])

        self.micromet = Micromet(self.z, self.lad, self.hc, cpara['micromet'])

        self.interception = Interception(cpara['interception'], self.lad * self.dz)

        self.forestfloor = ForestFloor(cpara['forestfloor'])

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
        for pt in self.planttypes:
            pt.update_daily(doy, Ta, PsiL=PsiL)  # updates pt properties
        self.LAI = sum([pt.LAI for pt in self.planttypes])  # total leaf area index [m2 m-2]
        self.lad = sum([pt.lad for pt in self.planttypes])  # total leaf area density [m2 m-3]

        """ normalized flow statistics in canopy with new lad """
        if self.Switch_Eflow and self.planttypes[0].Switch_lai:
            self.micromet.normalized_flow_stats(self.z, self.lad, self.hc)

    def run_timestep(self, dt, forcing):
        r""" Calculates one timestep and updates state of CanopyModel object.

        Args:
            dt: timestep [s]
            forcing (dataframe): meteorological and soil forcing data  !! NOT UP TO DATE
                'precipitation': precipitation rate [m s-1]
                'air_temperature': air temperature [\ :math:`^{\circ}`\ C]
                'dir_par': direct fotosynthetically active radiation [W m-2]
                'dif_par': diffuse fotosynthetically active radiation [W m-2]
                'dir_nir': direct near infrared radiation [W m-2]
                'dif_nir': diffuse near infrare active radiation [W m-2]
                'lw_in': Downwelling long wave radiation [W m-2]
                'wind_speed': wind speed [m s-1]
                'friction_velocity': friction velocity [m s-1]
                'co2': ambient CO2 mixing ratio [ppm]
                'h2o': ambient H2O mixing ratio [mol mol-1]
                'air_pressure': pressure [Pa]
                'zenith_angle': solar zenith angle [rad]
                'depth': [m] properties of first soil node
                'soil_temperature': [\ :math:`^{\circ}`\ C] properties of first soil node
                'soil_water_potential': [m] properties of first soil node
                'soil_volumetric_water': [m m\ :sup:`-3`\ ] properties of first soil node
                'hydraulic_conductivity': [m s\ :sup:`-1`\ ] properties of first soil node
                'thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ] properties of first soil node

        Returns:
            fluxes (dict)
            states (dict)
        """
        logger = logging.getLogger(__name__)

        forcing.update({'Ebal': self.Switch_Ebal,
                        'radiation': {}})
        ff_forcing = deepcopy(forcing)
        ff_forcing.update({'height': self.z[1],  # height to first calculation node
                           'nsteps': 20})

        """ --- flow stats --- """

        if self.Switch_Eflow is False:
            # recompute normalized flow statistics in canopy with current meteo
            self.micromet.normalized_flow_stats(
                    z=self.z,
                    lad=self.lad,
                    hc=self.hc,
                    Utop=forcing['wind_speed'] / (forcing['friction_velocity'] + EPS))
        # update U
        U = self.micromet.update_state(ustaro=forcing['friction_velocity'])

        """ --- SW profiles within canopy --- """

        ff_albedo = self.forestfloor.shortwave_albedo()
        forcing.update({'ff_albedo': ff_albedo,
                        'LAIz': self.lad * self.dz})

        # --- PAR ---
        forcing['radiation']['PAR'] = self.radiation.SW_profiles(
                forcing=forcing,
                radtype='PAR')
        f_sl = forcing['radiation']['PAR']['sunlit']['fraction']
        ff_forcing.update({'par': forcing['radiation']['PAR']['ground']})

        if self.Switch_Ebal:
            # --- NIR ---
            forcing['radiation']['NIR']= self.radiation.SW_profiles(
                forcing=forcing,
                radtype='NIR')
            ff_forcing.update({'nir': forcing['radiation']['NIR']['ground']})

            # absorbed radiation by leafs [W m-2(leaf)]
            forcing['radiation']['SW_absorbed'] = (
                        forcing['radiation']['PAR']['sunlit']['absorbed'] * f_sl
                        + forcing['radiation']['NIR']['sunlit']['absorbed'] * f_sl
                        + forcing['radiation']['PAR']['shaded']['absorbed'] * (1 - f_sl)
                        + forcing['radiation']['NIR']['shaded']['absorbed'] * (1 - f_sl))

        """ --- iterative solution of H2O, CO2, T, Tleaf and Tsurf --- """

        max_err = 0.01  # maximum relative error
        max_iter = 25  # maximum iterations
        gam = 0.5  # weight for new value in iterations
        err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
        Switch_WMA = self.Switch_WMA

        # initialize state variables
        T, H2O, CO2, Tleaf = self._restore(forcing)
        sources = {'h2o': None,  # [mol m-3 s-1]
                   'co2': None,  # [umol m-3 s-1]
                   'sensible_heat': None,  # [W m-3]
                   'latent_heat': None,  # [W m-3]
                   'fr': None}  # [W m-3]

        iter_no = 0
        while (err_t > max_err or err_h2o > max_err or
               err_co2 > max_err or err_Tl > max_err or
               err_Ts > max_err) and iter_no <= max_iter:

            iter_no += 1
            Tleaf_prev = Tleaf.copy()
            Tsurf_prev = self.forestfloor.temperature

            forcing.update({'leaf_temperature': Tleaf_prev,
                            'iteration': iter_no})

            if self.Switch_Ebal:
                """ ---  LW profiles within canopy --- """
                ff_longwave = self.forestfloor.longwave_radiation()

                forcing.update({'ff_longwave': ff_longwave})

                forcing['radiation']['LW'] = self.radiation.LW_profiles(
                        forcing=forcing)

                ff_forcing.update({'lw_dn': forcing['radiation']['LW']['down'][0],
                                   'lw_up': forcing['radiation']['LW']['up'][0]})

            # --- heat, h2o and co2 source terms
            for key in sources.keys():
                sources[key] = 0.0 * self.ones

            """ --- wet leaf water and energy balance --- """
            wetleaf_fluxes = self.interception.run(
                    dt=dt,
                    H2O=H2O, U=U, T=T,
                    forcing=forcing)
            # dry leaf fraction
            df = self.interception.df

            # update source terms
            for key in wetleaf_fluxes['sources'].keys():
                sources[key] += wetleaf_fluxes['sources'][key] / self.dz

            # canopy leaf temperature
            Tleaf = self.interception.Tl_wet * (1 - df) * self.lad

            """---  dry leaf gas-exchange --- """
            pt_stats = []
            for pt in self.planttypes:

                # --- sunlit and shaded leaves
                pt_stats_i, pt_sources = pt.leaf_gasexchange(
                        H2O, CO2, T, U, df=df, forcing=forcing)

                # update source terms
                # Dictionary IS modiefied in a loop. Wrapped in a list.
                for key in pt_sources.keys():
                    sources[key] += pt_sources[key]

                # append results
                pt_stats.append(pt_stats_i)

                # canopy leaf temperature
                Tleaf += pt_stats_i['Tleaf'] * df * pt.lad

            # canopy leaf temperature as weighted average
            Tleaf = Tleaf / (self.lad + EPS)

            err_Tl = max(abs(Tleaf - Tleaf_prev))

            """ --- forest floor --- """

            ff_forcing.update({'throughfall_rain': wetleaf_fluxes['throughfall_rain'],
                               'throughfall_snow': wetleaf_fluxes['throughfall_snow'],
                               'air_temperature': T[1],
                               'wind_speed': U[1],  # windspeed above forestfloor ca. 30 cm
                               'h2o': H2O[1],
                               'iteration': iter_no})

            fluxes_ffloor, states_ffloor = self.forestfloor.run(
                    dt=dt,
                    forcing=ff_forcing)

            err_Ts = abs(Tsurf_prev - self.forestfloor.temperature)

# check the sign of photosynthesis
            Fc_gr = -fluxes_ffloor['photosynthesis_bryo'] + fluxes_ffloor['respiration']

            """  --- solve scalar profiles (H2O, CO2, T) """
            if Switch_WMA is False:
                # to recognize oscillation
                if iter_no > 1:
                    T_prev2 = T_prev.copy()
                T_prev = T.copy()

                H2O, CO2, T, err_h2o, err_co2, err_t = self.micromet.scalar_profiles(
                        gam, H2O, CO2, T, forcing['air_pressure'],
                        source=sources,
                        lbc={'H2O': fluxes_ffloor['evaporation'],
                             'CO2': Fc_gr,
                             'T': fluxes_ffloor['sensible_heat_flux']},
                        Ebal=self.Switch_Ebal)

                # to recognize oscillation
                if iter_no > 5 and np.mean((T_prev - T)**2) > np.mean((T_prev2 - T)**2):
                    T = (T_prev + T) / 2
                    gam = max(gam / 2, 0.25)

                if (iter_no == max_iter or any(np.isnan(T)) or
                    any(np.isnan(H2O)) or any(np.isnan(CO2))):

                    if (any(np.isnan(T)) or any(np.isnan(H2O)) or any(np.isnan(CO2))):
                        logger.debug('%s Solution of profiles blowing up, T nan %s, H2O nan %s, CO2 nan %s',
                                         forcing['date'],
                                         any(np.isnan(T)), any(np.isnan(H2O)), any(np.isnan(CO2)))
                    elif max(err_t, err_h2o, err_co2, err_Tl, err_Ts) < 0.05:
                        if max(err_t, err_h2o, err_co2, err_Tl, err_Ts) > 0.01:
                            logger.debug('%s Maximum iterations reached but error tolerable < 0.05',
                                         forcing['date'])
                        break

                    Switch_WMA = True  # if no convergence, re-compute with WMA -assumption

                    logger.debug('%s Switched to WMA assumption: err_T %.4f, err_H2O %.4f, err_CO2 %.4f, err_Tl %.4f, err_Ts %.4f',
                                 forcing['date'],
                                 err_t, err_h2o, err_co2, err_Tl, err_Ts)

                    # reset values
                    iter_no = 0
                    err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
                    T, H2O, CO2, Tleaf = self._restore(forcing)
            else:
                err_h2o, err_co2, err_t = 0.0, 0.0, 0.0

        """ --- update state variables --- """
        self.interception.update()  # interception storage
        self.forestfloor.update()  # update old states to new states

        # --- ecosystem fluxes ---

        flux_co2 = (np.cumsum(sources['co2']) * self.dz 
                    + Fc_gr)  # [umol m-2 s-1]
        flux_latent_heat = (np.cumsum(sources['latent_heat']) * self.dz
                            + fluxes_ffloor['latent_heat_flux'])  # [W m-2]
        flux_sensible_heat = (np.cumsum(sources['sensible_heat']) * self.dz
                              + fluxes_ffloor['sensible_heat_flux'])  # [W m-2]

        # net ecosystem exchange [umol m-2 s-1]
        NEE = flux_co2[-1]
        # ecosystem respiration [umol m-2 s-1]
        Reco = sum([pt_st['dark_respiration'] for pt_st in pt_stats]) + fluxes_ffloor['respiration']
        # ecosystem GPP [umol m-2 s-1]
        GPP = - NEE + Reco
        # stand transpiration [m s-1]
        Tr = sum([pt_st['transpiration'] * MOLAR_MASS_H2O * 1e-3 for pt_st in pt_stats])

        if self.Switch_Ebal:
            # energy closure of canopy  -- THIS IS EQUAL TO frsource (the error caused by linearizing sigma*ef*T^4)
            energy_closure =  sum((forcing['radiation']['SW_absorbed'] + 
                                   forcing['radiation']['LW']['net_leaf']) * self.lad * self.dz) - (  # absorbed radiation
                              sum(sources['sensible_heat'] * self.dz)  # sensible heat
                              + sum(sources['latent_heat'] * self.dz))  # latent heat

        # --- RESULTS ---

        fluxes_ffloor.update({
                'potential_infiltration': fluxes_ffloor['potential_infiltration'],
                'evaporation_bryo': fluxes_ffloor['evaporation_bryo'] * MOLAR_MASS_H2O * 1e-3,  # [m s-1]
                'evaporation_litter': fluxes_ffloor['evaporation_litter'] * MOLAR_MASS_H2O * 1e-3,  # [m s-1]
                'evaporation_soil': fluxes_ffloor['evaporation_soil'] * MOLAR_MASS_H2O * 1e-3,  # [m s-1]
                'evaporation': fluxes_ffloor['evaporation'] * MOLAR_MASS_H2O * 1e-3  # [m s-1]
                })

        # return state and fluxes in dictionary
        state_canopy = {
                'interception_storage': sum(self.interception.W),
                'LAI': self.LAI,
                'lad': self.lad,
                'sunlit_fraction': f_sl,
                'phenostate': sum([pt.LAI * pt.pheno_state for pt in self.planttypes])/self.LAI,
                'IterWMA': iter_no
                }

        fluxes_canopy = {
                'throughfall': wetleaf_fluxes['throughfall'],
                'interception': wetleaf_fluxes['interception'],
                'evaporation': wetleaf_fluxes['evaporation'],
                'condensation': wetleaf_fluxes['condensation'],
                'transpiration': Tr,
                'NEE': NEE,
                'GPP': GPP,
                'respiration': Reco,
                'LE': flux_latent_heat[-1],
                'co2_flux': flux_co2,  # [umol m-2 s-1]
                'latent_heat_flux': flux_latent_heat,  # [W m-2]
                'pt_transpiration': np.array([pt_st['transpiration'] * MOLAR_MASS_H2O * 1e-3 for pt_st in pt_stats]),
                'pt_gpp': np.array([pt_st['net_co2'] + pt_st['dark_respiration'] for pt_st in pt_stats]),
                'pt_respiration': np.array([pt_st['dark_respiration'] for pt_st in pt_stats]),
                'pt_stomatal_conductance_h2o':  np.array([pt_st['stomatal_conductance'] for pt_st in pt_stats]),
                'pt_boundary_conductance_h2o':  np.array([pt_st['boundary_conductance'] for pt_st in pt_stats]),
                'pt_leaf_internal_co2':  np.array([pt_st['leaf_internal_co2'] for pt_st in pt_stats]),
                'pt_leaf_surface_co2':  np.array([pt_st['leaf_surface_co2'] for pt_st in pt_stats]),
                'water_closure': wetleaf_fluxes['water_closure'],
                }

        if self.Switch_WMA is False:
            state_canopy.update({'h2o': H2O,
                          'co2': CO2,
                          'temperature': T,
                          'WMA_assumption': 1.0*Switch_WMA})

        if self.Switch_Ebal:

            Tleaf_sl = np.where(self.lad > 0.0,
                                sum([pt_st['Tleaf_sl'] for pt_st in pt_stats]) / (self.lad + EPS),
                                np.nan)
            Tleaf_sh = np.where(self.lad > 0.0,
                                sum([pt_st['Tleaf_sh'] for pt_st in pt_stats]) / (self.lad + EPS),
                                np.nan)
            Tleaf_wet = np.where(self.lad > 0.0,
                                 self.interception.Tl_wet,
                                 np.nan)

            state_canopy.update({
                    'Tleaf_wet': Tleaf_wet,
                    'Tleaf_sl': Tleaf_sl,
                    'Tleaf_sh': Tleaf_sh,
                    'Tleaf': np.where(self.lad > 0.0, Tleaf, np.nan)
                    })

            fluxes_canopy.update({
                    'leaf_SW_absorbed': forcing['radiation']['SW_absorbed'],
                    'leaf_net_LW': forcing['radiation']['LW']['net_leaf'],
                    'sensible_heat_flux': flux_sensible_heat,  # [W m-2]
                    'energy_closure': energy_closure,
                    'fr_source': sum(sources['fr'] * self.dz)})

        return fluxes_canopy, state_canopy, fluxes_ffloor, states_ffloor

    def _restore(self, forcing):
        """ initialize state variables """

        T = self.ones * ([forcing['air_temperature']])
        H2O = self.ones * ([forcing['h2o']])
        CO2 = self.ones * ([forcing['co2']])
        Tleaf = T.copy() * self.lad / (self.lad + EPS)
        self.forestfloor.restore()
        if self.Switch_Ebal is False:
            self.forestfloor.temperature = (forcing['soil_temperature']
                                            + T[0]) / 2.0
        self.interception.Tl_wet = T.copy()
        for pt in self.planttypes:
            pt.Tl_sh = T.copy()
            pt.Tl_sl = T.copy()

        return T, H2O, CO2, Tleaf

# EOF
