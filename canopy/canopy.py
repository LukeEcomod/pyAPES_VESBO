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

LAST EDITS: Samuli Launiainen 15.01.2020
    - new forestfloor: in canopy changes in:
        * ff_forcing
        * forestfloor optical properties taken from previous timestep
        * ff_longwave computed in canopy
        * _restore() returns also Tsurf = forestfloor temperature in previous timestep
    - units of water fluxes [m s-1] -->  [kg m-2 s-1]
        -> * forcing['precipitation'] now in [kg m-2 s-1]
         * interception updated accordingly
    - cleaned forcing
    - returns altered
TODO:
    - clean interception -code

"""

import logging
import numpy as np
from .constants import MOLAR_MASS_H2O, EPS, STEFAN_BOLTZMANN, DEG_TO_KELVIN, WATER_DENSITY
from .radiation import Radiation
from .micromet import Micromet
from .interception import Interception
from .planttype.planttype import PlantType
from .forestfloor.forestfloor import ForestFloor

logger = logging.getLogger(__name__)

class CanopyModel(object):
    r""" Radiation, momentum, CO2 and energy exchange inside plant canopies.
    """

    def __init__(self, cpara, dz_soil):
        r""" Initializes canopy object and submodel objects using given parameters.

        Args:
            cpara (dict):
                'ctr' (dict): switches and specifications for computation
                    'Eflow' (bool): ensemble flow assumption (False solves U_normed on daily timestep)
                    'WMA' (bool): well-mixed assumption (False solves H2O, CO2, T)
                    'Ebal' (bool): True solves energy balance
                    'WaterStress' (str): account for water stress in planttypes via 'Rew', 'PsiL' or None omits
                    'seasonal_LAI' (bool): account for seasonal LAI dynamics
                    'pheno_cycle' (bool): account for phenological cycle
                'grid' (dict):
                    'zmax': heigth of grid from ground surface [m]
                    'Nlayers': number of layers in grid [-]
                'radiation' (dict): radiation model parameters
                'micromet' (dict): micromet model parameters
                'interception' (dict): interception model parameters
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
                .interception (object): interception model
                .forestfloor (object): forest floor object (bryotype/baresoil/snow)
        """

        # --- grid ---
        self.z = np.linspace(0, cpara['grid']['zmax'], cpara['grid']['Nlayers'])  # grid [m] above ground
        self.dz = self.z[1] - self.z[0]  # gridsize [m]
        self.ones = np.ones(len(self.z))  # dummy

        # --- switches ---
        self.Switch_Eflow = cpara['ctr']['Eflow'] # True assumes constant U/ustar at upper boundary
        self.Switch_WMA = cpara['ctr']['WMA']   # True solves scalar profiles
        self.Switch_Ebal = cpara['ctr']['Ebal'] # True solves leaf energy balance

        logger.info('Eflow: %s, WMA: %s, Ebal: %s',
                    self.Switch_Eflow,
                    self.Switch_WMA,
                    self.Switch_Ebal)

        # --- PlantTypes ---
        ptypes = []
        ptnames = list(cpara['planttypes'].keys())
        ptnames.sort()
        for pt in ptnames:
            ptypes.append(PlantType(self.z, cpara['planttypes'][pt], dz_soil, ctr=cpara['ctr'], loc=cpara['loc']))
        self.planttypes = ptypes

        # --- stand characteristics: sum over planttypes---

        # total leaf area index [m2 m-2]
        self.LAI = sum([pt.LAI for pt in self.planttypes])
        # total leaf area density [m2 m-3]
        self.lad = sum([pt.lad for pt in self.planttypes])

         # layerwise mean leaf characteristic dimension [m] for interception model
        self.leaf_length = sum([pt.leafp['lt'] * pt.lad for pt in self.planttypes]) / (self.lad + EPS)

        # root area density [m2 m-3]
        rad = np.zeros(np.shape(dz_soil))
        imax = 1
        for pt in self.planttypes:
            rad[:len(pt.Roots.rad)] += pt.Roots.rad
            imax = max(imax, len(pt.Roots.rad))

        self.ix_roots = np.array(range(imax)) # soil model layers corresponding to root zone
        self.rad = rad[self.ix_roots]

        # total root area index [m2 m-2]
        self.RAI = sum([pt.Roots.RAI for pt in self.planttypes])
        # distribution of roots across soil model layers [-]
        self.root_distr = self.rad * dz_soil[self.ix_roots] / (self.RAI + EPS)

        # canopy height [m]
        if len(np.where(self.lad > 0)[0]) > 0:
            f = np.where(self.lad > 0)[0][-1]
            self.hc = self.z[f].copy()
        else:
            self.hc = 0.0

        # --- create radiation, micromet, interception, and forestfloor model instances
        self.radiation = Radiation(cpara['radiation'], self.Switch_Ebal)

        self.micromet = Micromet(self.z, self.lad, self.hc, cpara['micromet'])

        self.interception = Interception(cpara['interception'], self.lad * self.dz)

        # forestfloor. Now soil respiration profile is as root_distr
        self.forestfloor = ForestFloor(cpara['forestfloor'],
                                       respiration_profile=self.root_distr)

    def run_daily(self, doy, Ta, Rew=1.0):
        r""" Computatations at daily timestep.
        Updates planttypes and total canopy leaf area index and phenological state.
        Recomputes normalize flow statistics with new leaf area density profile.

        Args:
            doy (float): day of year [days]
            Ta (float): mean daily air temperature [degC]
            PsiL (float): leaf water potential [MPa] --- CHECK??
            Rew (float): relatively extractable water (-)
        """
        """ update physiology and leaf area of planttypes and canopy"""
        for pt in self.planttypes:
            if pt.LAImax > 0.0:
                PsiL = (pt.Roots.h_root - self.z) / 100.0  # MPa
                pt.update_daily(doy, Ta, PsiL=PsiL, Rew=Rew)  # updates pt properties

        # total leaf area index [m2 m-2]
        self.LAI = sum([pt.LAI for pt in self.planttypes])
        # total leaf area density [m2 m-3]
        self.lad = sum([pt.lad for pt in self.planttypes])
         # layerwise mean leaf characteristic dimension [m]
        self.leaf_length = sum([pt.leafp['lt'] * pt.lad for pt in self.planttypes]) / (self.lad + EPS)

        """ normalized flow statistics in canopy with new lad """
        if self.Switch_Eflow and self.planttypes[0].Switch_lai:
            self.micromet.normalized_flow_stats(self.z, self.lad, self.hc)

    def run(self, dt, forcing, parameters):
        r""" Calculates one timestep and updates state of CanopyModel object.
        Args:
            dt: timestep [s]
            forcing (dict): meteorological forcing at ubc and soil state at uppermost soil layer.
                'precipitation': precipitation rate [kg m-2 s-1]
                'air_temperature': air temperature [\ :math:`^{\circ}`\ C]
                'wind_speed': wind speed [m s-1]
                'friction_velocity': friction velocity [m s-1]
                'co2': ambient CO2 mixing ratio [ppm]
                'h2o': ambient H2O mixing ratio [mol mol-1]
                'air_pressure': pressure [Pa]
                'zenith_angle': solar zenith angle [rad]
                'PAR' (dict): with keys 'direct', 'diffuse' [Wm-2]
                'NIR' (dict): with keys 'direct', 'diffuse' [Wm-2]

                'soil_temperature': [\ :math:`^{\circ}`\ C] properties of first soil node
                'soil_water_potential': [m] properties of first soil node
                'soil_volumetric_water': [m m\ :sup:`-3`\ ] properties of first soil node
                'soil_volumetric_air': [m m\ :sup:`-3`\ ] properties of first soil node
                'soil_pond_storage': [kg m-2] properties of first soil node
            'parameters':
                'date'
                'thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ] properties of first soil node
                'hydraulic_conductivity': [m s\ :sup:`-1`\ ] properties of first soil node
                'depth': [m] properties of first soil node

        Returns:
            fluxes (dict)
            states (dict)
        """
        logger = logging.getLogger(__name__)

        # --- solve Canopy flow ---

        if self.Switch_Eflow is False:
            # recompute normalized flow statistics in canopy with current meteo
            self.micromet.normalized_flow_stats(
                z=self.z,
                lad=self.lad,
                hc=self.hc,
                Utop=forcing['wind_speed'] / (forcing['friction_velocity'] + EPS))

        # get U and ustar profiles [m s-1]
        U, ustar = self.micromet.update_state(ustaro=forcing['friction_velocity'])

        # --- solve SW profiles within canopy ---

        radiation_profiles = {}

        radiation_params = {
            'ff_albedo': self.forestfloor.albedo,
            'LAIz': self.lad * self.dz,
        }

        # --- PAR ---
        radiation_params.update({
            'radiation_type': 'par',
        })

        radiation_profiles['par'] = self.radiation.shortwave_profiles(
            forcing=forcing,
            parameters=radiation_params
        )

        sunlit_fraction = radiation_profiles['par']['sunlit']['fraction']

        if self.Switch_Ebal:
            # --- NIR ---
            radiation_params['radiation_type'] = 'nir'

            radiation_profiles['nir'] = self.radiation.shortwave_profiles(
                forcing=forcing,
                parameters=radiation_params
            )

            # absorbed radiation by leafs [W m-2(leaf)]
            radiation_profiles['sw_absorbed'] = (
                radiation_profiles['par']['sunlit']['absorbed'] * sunlit_fraction
                + radiation_profiles['nir']['sunlit']['absorbed'] * sunlit_fraction
                + radiation_profiles['par']['shaded']['absorbed'] * (1. - sunlit_fraction)
                + radiation_profiles['nir']['shaded']['absorbed'] * (1. - sunlit_fraction)
            )

        # --- start iterative solution of H2O, CO2, T, Tleaf and Tsurf ---

        max_err = 0.01  # maximum relative error
        max_iter = 25  # maximum iterations
        gam = 0.5  # weight for new value in iterations
        err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
        Switch_WMA = self.Switch_WMA

        # initialize state variables
        T, H2O, CO2, Tleaf, Tsurf = self._restore(forcing)
        sources = {
            'h2o': None,  # [mol m-3 s-1]
            'co2': None,  # [umol m-3 s-1]
            'sensible_heat': None,  # [W m-3]
            'latent_heat': None,  # [W m-3]
            'fr': None  # [W m-3]
        }

        iter_no = 0
        while (err_t > max_err or err_h2o > max_err or
               err_co2 > max_err or err_Tl > max_err or
               err_Ts > max_err) and iter_no <= max_iter:

            iter_no += 1
            Tleaf_prev = Tleaf.copy()
            Tsurf_prev = Tsurf

            if self.Switch_Ebal:
                # ---  LW profiles within canopy ---

                # emitted longwave from forest floor
                lw_surf = self.forestfloor.emissivity * STEFAN_BOLTZMANN *\
                                np.power(Tsurf + DEG_TO_KELVIN, 4)
                lw_forcing = {
                    'lw_in': forcing['lw_in'],
                    'lw_up': lw_surf,
                    'leaf_temperature': Tleaf_prev,
                }

                lw_params = {
                    'LAIz': self.lad * self.dz,
                    'ff_emissivity': self.forestfloor.emissivity
                }

                # solve LW profiles
                radiation_profiles['lw'] = self.radiation.longwave_profiles(
                    forcing=lw_forcing,
                    parameters=lw_params
                )

            # --- heat, h2o and co2 source terms
            for key in sources.keys():
                sources[key] = 0.0 * self.ones

            # --- interception of rainfall & wet leaf water and energy balance ---
            interception_forcing = {
                'h2o': H2O,
                'wind_speed': U,
                'air_temperature': T,
                'air_pressure': forcing['air_pressure'],
                'leaf_temperature': Tleaf_prev,
                'precipitation': forcing['precipitation'],
            }

            if self.Switch_Ebal:
                interception_forcing.update({
                    'sw_absorbed': radiation_profiles['sw_absorbed'],
                    'lw_radiative_conductance': radiation_profiles['lw']['radiative_conductance'],
                    'net_lw_leaf': radiation_profiles['lw']['net_leaf'],
                })

            interception_params = {
                'LAIz': self.lad * self.dz,
                'leaf_length': self.leaf_length
            }

            interception_controls = {
                'energy_balance': self.Switch_Ebal,
                'logger_info': 'date: {} iteration: {}'.format(
                    parameters['date'],
                    iter_no
                )
            }

            if self.Switch_Ebal:
                interception_forcing.update({
                    'sw_absorbed': radiation_profiles['sw_absorbed'],
                    'lw_radiative_conductance': radiation_profiles['lw']['radiative_conductance'],
                    'net_lw_leaf': radiation_profiles['lw']['net_leaf'],
                })

            # --- solve interception model
            wetleaf_fluxes = self.interception.run(
                dt=dt,
                forcing=interception_forcing,
                parameters=interception_params,
                controls=interception_controls
            )

            # dry leaf fraction
            df = self.interception.df

            # update source terms
            for key in wetleaf_fluxes['sources'].keys():
                sources[key] += wetleaf_fluxes['sources'][key] / self.dz

            # canopy leaf temperature
            Tleaf = self.interception.Tl_wet * (1 - df) * self.lad

            # --- dry leaf gas-exchange: solve for each planttype ---
            forcing_pt = {
                'h2o': H2O,
                'co2': CO2,
                'air_temperature': T,
                'air_pressure': forcing['air_pressure'],
                'wind_speed': U,
                'par': radiation_profiles['par'],
                'average_leaf_temperature': Tleaf_prev, # sl changed 25.11.
                'wet_leaf_temperature': self.interception.Tl_wet # kh for wet leaf rd
            }

            if self.Switch_Ebal:
                forcing_pt.update({
                    'nir': radiation_profiles['nir'],
                    'lw': radiation_profiles['lw']
                })

            parameters_pt = {
                'dry_leaf_fraction': self.interception.df,
                'sunlit_fraction': sunlit_fraction
            }

            controls_pt = {
                'energy_balance': self.Switch_Ebal,
                'logger_info': 'date: {} iteration: {}'.format(
                    parameters['date'],
                    iter_no
                )
            }

            pt_stats = []
            pt_layerwise = []

            for pt in self.planttypes:

                # --- sunlit and shaded leaves gas-exchange and optional energy balance
                pt_stats_i, layer_stats_i = pt.run(
                    forcing=forcing_pt,
                    parameters=parameters_pt,
                    controls=controls_pt
                )

                # update source terms: now do explicitly (SL 17.11.19)
                sources['co2'] -= layer_stats_i['net_co2'] # in planttype, net uptake >0, here we want sink negative
                sources['h2o'] += layer_stats_i['transpiration']
                sources['sensible_heat'] += layer_stats_i['sensible_heat']
                sources['latent_heat'] += layer_stats_i['latent_heat']
                sources['fr'] += layer_stats_i['fr']

                # update canopy leaf temperature; it must be lad-weighted average
                Tleaf += layer_stats_i['leaf_temperature'] * df * pt.lad

                # append results
                pt_stats.append(pt_stats_i)

                layer_stats_i['leaf_temperature'] *= pt.mask # set pt.lad=0 to np.NaN
                pt_layerwise.append(layer_stats_i)
            del pt, pt_stats_i, layer_stats_i

            # --- end of planttype -loop

            # mean leaf temperature of canopy layer
            Tleaf = Tleaf / (self.lad + EPS)

            err_Tl = max(abs(Tleaf - Tleaf_prev))

            # --- solve forest floor water & heat balance & carbon exchange ---

            # compile control and forcing dicts
            ff_controls = {
                'energy_balance': self.Switch_Ebal,
                'logger_info': 'date: {} iteration: {}'.format(
                    parameters['date'],
                    iter_no
                )
            }

            ff_params = {
                'reference_height': self.z[1],  # height of first node above ground
                'soil_depth': parameters['soil_depth'],
                'soil_hydraulic_conductivity': parameters['soil_hydraulic_conductivity'][0],
                'soil_thermal_conductivity': parameters['soil_thermal_conductivity'],
                'iteration': iter_no,
            }

            ff_forcing = {
                'precipitation_rain': wetleaf_fluxes['throughfall_rain'],
                'precipitation_snow': wetleaf_fluxes['throughfall_snow'],
                'par': radiation_profiles['par']['ground'],
                'air_temperature': T[1],
                'h2o': H2O[1],
                'co2': CO2[1],
                'air_pressure': forcing['air_pressure'],
                'wind_speed': U[1],
                'friction_velocity': ustar[1],
                'soil_temperature': forcing['soil_temperature'],
                'soil_water_potential': forcing['soil_water_potential'][0],
                'soil_volumetric_water': forcing['soil_volumetric_water'],
                'soil_volumetric_air': forcing['soil_volumetric_air'],
                'soil_pond_storage': forcing['soil_pond_storage']
            }

            if self.Switch_Ebal:

                ff_forcing.update({
                    'nir': radiation_profiles['nir']['ground'],
                    'lw_dn': radiation_profiles['lw']['down'][0],
                })

            # --- solve forestfloor
            ff_fluxes, ff_states, gt_results = self.forestfloor.run(
                        dt=dt,
                        forcing=ff_forcing,
                        parameters=ff_params,
                        controls=ff_controls
                        )

            Tsurf = ff_states['surface_temperature']
            err_Ts = abs(Tsurf_prev - Tsurf)

            # --- compute scalar profiles (H2O, CO2, T) in canopy air-space
            if Switch_WMA is False:
                # to recognize oscillation
                if iter_no > 1:
                    T_prev2 = T_prev.copy()
                T_prev = T.copy()

                # --- solve scalar profiles
                H2O, CO2, T, err_h2o, err_co2, err_t = self.micromet.scalar_profiles(
                        gam, H2O, CO2, T, forcing['air_pressure'],
                        source=sources,
                        lbc={'H2O': ff_fluxes['evaporation'] / MOLAR_MASS_H2O,
                             'CO2': ff_fluxes['net_co2'],
                             'T': ff_fluxes['sensible_heat']},
                        Ebal=self.Switch_Ebal)

                # to recognize oscillation
                if iter_no > 5 and np.mean((T_prev - T)**2) > np.mean((T_prev2 - T)**2):
                    T = (T_prev + T) / 2
                    gam = max(gam / 2, 0.25)

                if (iter_no == max_iter or any(np.isnan(T)) or
                    any(np.isnan(H2O)) or any(np.isnan(CO2))):

                    if (any(np.isnan(T)) or any(np.isnan(H2O)) or any(np.isnan(CO2))):
                        logger.debug('%s Solution of profiles blowing up, T nan %s, H2O nan %s, CO2 nan %s',
                                         parameters['date'],
                                         any(np.isnan(T)), any(np.isnan(H2O)), any(np.isnan(CO2)))
                    elif max(err_t, err_h2o, err_co2, err_Tl, err_Ts) < 0.05:
                        if max(err_t, err_h2o, err_co2, err_Tl, err_Ts) > 0.01:
                            logger.debug('%s Maximum iterations reached but error tolerable < 0.05',
                                         parameters['date'])
                        break

                    Switch_WMA = True  # if no convergence, re-compute with WMA -assumption

                    logger.debug('%s Switched to WMA assumption: err_T %.4f, err_H2O %.4f, err_CO2 %.4f, err_Tl %.4f, err_Ts %.4f',
                                 parameters['date'],
                                 err_t, err_h2o, err_co2, err_Tl, err_Ts)

                    # reset values
                    iter_no = 0
                    err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
                    T, H2O, CO2, Tleaf, Tsurf = self._restore(forcing)
            else:
                err_h2o, err_co2, err_t = 0.0, 0.0, 0.0

        # --- end of iterative solution of timestep

        # --- update state variables
        self.interception.update()
        self.forestfloor.update()

        # --- Compile outputs

        # --- integrate to ecosystem fluxes (per m-2 ground) ---
        flux_co2 = (np.cumsum(sources['co2']) * self.dz + ff_fluxes['net_co2'])  # [umol m-2 s-1]
        flux_latent_heat = (np.cumsum(sources['latent_heat']) * self.dz + ff_fluxes['latent_heat'])  # [W m-2]
        flux_sensible_heat = (np.cumsum(sources['sensible_heat']) * self.dz + ff_fluxes['sensible_heat'])  # [W m-2]

        # net ecosystem exchange [umol m-2 s-1]
        NEE = flux_co2[-1]
        # ecosystem respiration [umol m-2 s-1]
        Reco = sum([pt_st['dark_respiration'] for pt_st in pt_stats]) + ff_fluxes['respiration']
        # ecosystem GPP [umol m-2 s-1]
        GPP = - NEE + Reco
        # stand transpiration [m s-1]
        Tr = sum([pt_st['transpiration'] * MOLAR_MASS_H2O * 1e-3 for pt_st in pt_stats])

        # divide total canopy root water uptake into soil model layers [m s-1]
        rootsink = np.zeros(np.shape(self.rad))
        pt_index = 0
        for pt in self.planttypes:
            if pt.LAImax > 0.0:
                Tr_pt = pt_stats[pt_index]['transpiration'] * MOLAR_MASS_H2O / WATER_DENSITY
                rootsink[pt.Roots.ix] += pt.Roots.wateruptake(
                        transpiration=Tr_pt,
                        h_soil=forcing['soil_water_potential'],
                        kh_soil=parameters['soil_hydraulic_conductivity'])
            pt_index += 1
        del pt_index, pt

        if self.Switch_Ebal:
            # check energy closure of canopy
            # -- THIS IS EQUAL TO frsource (the error caused by linearizing sigma*ef*T^4)
            energy_closure =  sum((radiation_profiles['sw_absorbed'] +
                                   radiation_profiles['lw']['net_leaf']) * self.lad * self.dz) - (# absorbed radiation
                              sum(sources['sensible_heat'] * self.dz)  # sensible heat
                              + sum(sources['latent_heat'] * self.dz))  # latent heat

        # return state and fluxes in dictionaries: to save into results, name convention and dimensions must
        # follow those in 'general parameters'

        outputs_canopy = {

                # canopy state
                'LAI': self.LAI, # m2m-2
                'lad': self.lad, # m2 m-3
                'IterWMA': iter_no,
                'WMA_assumption': 1.0*Switch_WMA,
                'phenostate': sum([pt.LAI * pt.pheno_state for pt in self.planttypes])/(self.LAI + EPS),

                # micromet profiles
                'wind_speed': U,    # [m s-1]
                'friction_velocity': ustar, # [m s-1]
                'h2o': H2O, # [mol mol-1]
                'co2': CO2, # [ppm]
                'temperature': T, # [degC]

                # radiation profiles
                'sunlit_fraction': sunlit_fraction, # [-]
                'par_down': radiation_profiles['par']['down'], # [W m-2]
                'par_up': radiation_profiles['par']['up'], # [W m-2]
                'par_absorbed_sunlit': radiation_profiles['par']['sunlit']['absorbed'], #  [W m-2 (leaf)]
                'par_absorbed_shaded': radiation_profiles['par']['shaded']['absorbed'],

                # total fluxes from interception model [kg m-2 s-1 = mm s-1], divide with WATER_DENSITY to get [kg m-2 s-1 = mm s-1]
                'interception_storage': sum(self.interception.W),
                'throughfall': wetleaf_fluxes['throughfall'],
                'interception': wetleaf_fluxes['interception'],
                'evaporation': wetleaf_fluxes['evaporation'],
                'condensation': wetleaf_fluxes['condensation'],
                'condensation_drip': wetleaf_fluxes['condensation_drip'],
                'water_closure': wetleaf_fluxes['water_closure'],

                # vertical profiles from interception model [kg m-2 s-1 = mm s-1]
                'evaporation_ml': wetleaf_fluxes['evaporation_ml'],
                'throughfall_ml': wetleaf_fluxes['throughfall_ml'],
                'condensation_drip_ml': wetleaf_fluxes['condensation_drip_ml'],

                # ecosystem fluxes per m2 ground
                'SH': flux_sensible_heat[-1], # W m-2
                'LE': flux_latent_heat[-1], # Wm-2
                'NEE': NEE, # net ecosystem exchange [umol m-2 s-1]
                'GPP': GPP, # gross-primary productivity [umol m-2 s-1]
                'Reco': Reco, # ecosystem respiration [umol m-2 s-1]
                'transpiration': Tr, # transpiration of all planttypes [m s-1]

                # flux profiles
                'co2_flux': flux_co2,  # [umol m-2 s-1]
                'latent_heat_flux': flux_latent_heat,  # [W m-2]
                'sensible_heat_flux': flux_sensible_heat,  # [W m-2]

                # root sink profile for Soil-model
                'root_sink' : rootsink, # [m s-1]
                }

        if self.Switch_Ebal:
            # layer - averaged leaf temperature; average over all plant-types
            Tleaf_sl = np.where(self.lad > 0.0, sum([k['leaf_temperature_sunlit'] for k in pt_layerwise]), np.nan)
            Tleaf_sh = np.where(self.lad > 0.0, sum([k['leaf_temperature_shaded'] for k in pt_layerwise]), np.nan)
            Tleaf_wet = np.where(self.lad > 0.0, self.interception.Tl_wet, np.nan)

            outputs_canopy.update({
                    'Tleaf_wet': Tleaf_wet,
                    'Tleaf_sl': Tleaf_sl,
                    'Tleaf_sh': Tleaf_sh,
                    'Tleaf': np.where(self.lad > 0.0, Tleaf, np.nan)
                    })

            # net SW and LW radiation at uppermost layer.
            SWnet = (radiation_profiles['nir']['down'][-1] - radiation_profiles['nir']['up'][-1] +
                     radiation_profiles['par']['down'][-1] - radiation_profiles['par']['up'][-1])
            LWnet = (radiation_profiles['lw']['down'][-1] - radiation_profiles['lw']['up'][-1])

            outputs_canopy.update({
                    'sensible_heat_flux': flux_sensible_heat,  # [W m-2]
                    'energy_closure': energy_closure, # energy balance closure [Wm-2]
                    'SWnet': SWnet, # W m-2 (ground)
                    'LWnet': LWnet, # W m-2 (ground)
                    'Rnet': SWnet + LWnet,
                    'fr_source': sum(sources['fr'] * self.dz),

                    # outputs from radiation-model
                    'leaf_net_SW': radiation_profiles['sw_absorbed'], # [Wm-2 leaf]
                    'leaf_net_LW': radiation_profiles['lw']['net_leaf'], # [Wm-2 leaf]
                    'nir_down': radiation_profiles['nir']['down'], # [Wm-2 ground]
                    'nir_up': radiation_profiles['nir']['up'],
                    'lw_down': radiation_profiles['lw']['down'],
                    'lw_up': radiation_profiles['lw']['up']
                    })

        # plant-type specific results: this is dictionary of lists where each planttype value is element in list
        # integrated values over planttype
        pt_results = {
                'root_water_potential': np.array([pt.Roots.h_root for pt in self.planttypes]),
                'total_transpiration': np.array([pt_st['transpiration'] * MOLAR_MASS_H2O for pt_st in pt_stats]), #[kg m-2]
                'total_gpp': np.array([pt_st['net_co2'] + pt_st['dark_respiration'] for pt_st in pt_stats]),
                'total_dark_respiration': np.array([pt_st['dark_respiration'] for pt_st in pt_stats]),
                'total_stomatal_conductance_h2o':  np.array([pt_st['stomatal_conductance'] for pt_st in pt_stats]),
                'total_boundary_conductance_h2o':  np.array([pt_st['boundary_conductance'] for pt_st in pt_stats])
                }

        # add vertical profiles: convert list of dicts to dict of lists and append. Love Python!
        pt_profs = {}
        for k,v in pt_layerwise[0].items():
            pt_profs[k] = [x[k] for x in pt_layerwise]

        pt_results.update(pt_profs)
        del pt_profs

        # --- outputs from forestfloor: ff_fluxes and groudtype -specific resutls
        ff_fluxes.update(ff_states)
        ff_fluxes.update({'par_albedo': self.forestfloor.albedo['PAR'],
                          'nir_albedo': self.forestfloor.albedo['NIR']}
                        )
        return outputs_canopy, pt_results, ff_fluxes, gt_results

    def _restore(self, forcing):
        """ initialize state variables """

        T = self.ones * ([forcing['air_temperature']])
        H2O = self.ones * ([forcing['h2o']])
        CO2 = self.ones * ([forcing['co2']])
        Tleaf = T.copy() * self.lad / (self.lad + EPS)
        Tground = self.forestfloor.temperature
        self.interception.Tl_wet = T.copy()

        for pt in self.planttypes:
            pt.Tl_sh = T.copy()
            pt.Tl_sl = T.copy()

        return T, H2O, CO2, Tleaf, Tground

# EOF
