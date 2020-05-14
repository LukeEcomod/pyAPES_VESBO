#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module: forestfloor
    :synopsis: APES-model component
.. moduleauthor:: Antti-Jussi Kieloaho

Describes forest floor consisting of bottomlayer types, baresoil and snowpack.

canopy.forestfloor is interface that handles tiling of bottomlayer types.

DEVELOPMENT VERSION: Last edit SL 26.3.2020

Todo:
    - snowpack: energy balance and snowdepth changes with a process model
    - bare soil energy balance implementation
    - organic layer freezing/thawing
    - add heat advection (infiltration & capillary rise) to soil.heat!

Created on Tue Mar 13 11:59:23 2018
"""

from canopy.constants import EPS, MOLAR_MASS_H2O, LATENT_HEAT

from .organiclayer import OrganicLayer
from .snowpack import DegreeDaySnow
from .carbon import SoilRespiration

import logging
logger = logging.getLogger(__name__)

class ForestFloor(object):
    r"""Describes forest floor consisting of organic layers and possible snowpack.
    """

    def __init__(self, para, respiration_profile=None):
        r""" Initializes forestfloor object
        Args
            para (dict):
             org_types (dict):
                'name' (str)
                'layer_type': 'bryophyte' or 'litter'
                'coverage':  [-]
                'height': [m]
                'roughness_height': [m]
                #'leaf_area_index': [m\ :sup:`2` m :sup:`-2`\ ]
                #'specific_leaf_area': [m\ :sup:`3` m :sup:`-3`\ ]
                'dry_mass': [kg m\ :sup:`-2`]
                'bulk_density': [kg m\ :sup:`-3`]
                'max_water_content': [g g\ :sup:`-1`\ ]
                'max_symplast_water_content': [g g\ :sup:`-1`\ ]
                'min_water_content': [g g\ :sup:`-1`\ ]
                'water_retention' (dict):
                    #'theta_s': saturated water content [m\ :sup:`3` m :sup:`-3`\ ]
                    #'theta_r': residual water content [m\ :sup:`3` m :sup:`-3`\ ]
                    'alpha': air entry suction [cm\ :sup:`-1`]
                    'n': pore size distribution [-]
                    'saturated conductivity': [m s\ :sup:`-1`]
                    'pore connectivity': (l) [-]
                'porosity': [m\ :sup:`3` m\ :sup:`-3`\ ]
                'photosynthesis' (dict): : only if layer_type == 'bryophyte'
                    if farquhar-model as now:
                        'Vcmax', 'Jmax', 'Rd', 'alpha', 'theta', 'beta',
                        'gmax', 'wopt', 'a0', 'a1', 'CAP_desic', 'tresp'
                'respiration' (dict): only if layer_type == 'litter'
                    'q10' [-]
                    'r10' [\ :math:`\mu`\ mol m\ :sup:`-1`\ :sub:`leaf` s\ :sup:`-1`]
                'optical_properties' (dict):
                    'albedo_par': [-] photosynthetically active radiation (PAR)
                    'albedo_nir': [-] near infrared radiation (NIR)
                    'emissivity': [-]
                'initial_conditions' (dict)
            snowpack (dict):
                'kmelt' Melting coefficient [m degC-1 s-1] (=2.0 mm/C/d)
                'kfreeze': Freezing  coefficient [m degC-1 s-1] (=0.5 mm/C/d)
                'retention': max fraction of liquid water in snow [-]
                'Tmelt': temperature when melting starts [degC]
                'optical_properties':
                        'emissivity':
                        'albedo_PAR':
                        'albedo_NIR':
                'swe_initial': [kg m-2]
            soil_respiration (dict):
                    'r10': base respiration rate [umolm-2s-1]
                    'q10': temperature sensitivity [-]
                    'moisture_coeff' (list): moisture response parameters
        Returns:
            self (object)
        """

        # -- forest floor tiled surface of organic layers. snowpack can overly ground
        self.snowpack = DegreeDaySnow(para['snowpack'])

        self.soilrespiration = SoilRespiration(para['soil_respiration'], weights=respiration_profile)

        # Organiclayers; includes both bryophytes and litter. Append to list and
        # Compute area-weighted forest floor temperature and optical properties
        bltypes = []
        bltnames = list(para['bottom_layer_types'].keys())
        bltnames.sort()

        for bt in bltnames:
            if para['bottom_layer_types'][bt]['coverage'] > 0:
                # case coverage > 0:
                bltypes.append(OrganicLayer(para['bottom_layer_types'][bt]))
            else:
                logger.info('Forestfloor: %s, coverage 0.0 omitted!', bt)

        f_organic = sum([bt.coverage for bt in bltypes])
        if abs(f_organic - 1.0) > EPS:
            raise ValueError('The sum of organic type coverage must = 1, ' +
                             'now %.2f' % f_organic)

        self.bottomlayer_types = bltypes
        logger.info('Forestfloor has %s bottomlayer types', len(self.bottomlayer_types))

        if self.snowpack.swe > 0:
            self.temperature = self.snowpack.temperature
            self.surface_temperature = self.snowpack.temperature
            self.albedo = self.snowpack.optical_properties['albedo']
            self.emissivity = self.snowpack.optical_properties['emissivity']
        else:
            self.temperature = sum([bt.coverage * bt.temperature
                                    for bt in self.bottomlayer_types])
            self.surface_temperature = sum([bt.coverage * bt.surface_temperature
                                    for bt in self.bottomlayer_types])
            self.albedo = {'PAR': sum([bt.coverage * bt.albedo['PAR']
                                       for bt in self.bottomlayer_types]),
                           'NIR': sum([bt.coverage * bt.albedo['NIR']
                                       for bt in self.bottomlayer_types])}
            self.emissivity = sum([bt.coverage * bt.emissivity
                                   for bt in self.bottomlayer_types])
        self.water_storage = sum([bt.coverage * bt.water_storage
                                  for bt in self.bottomlayer_types])


    def update(self):
        """ Updates forestfloor states
        """
        self.snowpack.update()

        for bt in self.bottomlayer_types:
            bt.update_state()

        if self.snowpack.swe > 0:
            self.temperature = self.snowpack.temperature
            self.surface_temperature = self.snowpack.temperature
            self.albedo = self.snowpack.optical_properties['albedo']
            self.emissivity = self.snowpack.optical_properties['emissivity']
        else: # NOTE! forestfloor temperature in snow-free conditions is weighted average of moss surface temperature!
            self.temperature = sum([bt.coverage * bt.temperature
                                    for bt in self.bottomlayer_types])
            self.surface_temperature = sum([bt.coverage * bt.surface_temperature
                                            for bt in self.bottomlayer_types])
            self.albedo['PAR'] = sum([bt.coverage * bt.albedo['PAR']
                                      for bt in self.bottomlayer_types])
            self.albedo['NIR'] = sum([bt.coverage * bt.albedo['NIR']
                                      for bt in self.bottomlayer_types])
            self.emissivity = sum([bt.coverage * bt.emissivity
                                   for bt in self.bottomlayer_types])
            self.water_storage = sum([bt.coverage * bt.water_storage
                                      for bt in self.bottomlayer_types])

    def run(self, dt, forcing, parameters, controls):
        r"""Water and energy balance at the forestfloor; handles 'tiled' surfaces
        at forest floor and aggregates average fluxes and state.

        Args:
            dt: timestep [s]
            forcing (dict): states of microclimate
                'precipitation_rain': [kg m-2 s\ :sup:`-1`\ ]
                'precipitation_snow': [kg m-2 s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ] if energy_balance is True
                'lw_dn': [W m\ :sup:`-2`\ ] if energy_balance is True
                'h2o': [mol mol\ :sup:`-1`\ ]
                'co2': [ppm]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [Pa]
                'soil_volumetric_water': [m\ :sup:`3`\  m\`-3`\ ]
                'soil_volumetric_air': [m\ :sup:`3`\  m\`-3`\ ]
                'soil_pond_storage': [kg m-2]
            parameters (dict):
                'soil_thermal_conductivity': [] if energy_balance is True
                'soil_hydraulic_conductivity': []
                'depth': [m] first soil calculation node
                'reference_height': [m] first canopy calculation node
            controls (dict):
                'energy_balance': boolean
                'logger_info': str
        Returns:
            fluxes (dict): forestfloor aggregated fluxes
                'net_radiation' [W m-2]
                'sensible_heat'
                'latent_heat'
                'ground_heat'
                'energy_closure'

                'evaporation': tiles + soil below [kg m-2 s-1]
                'soil_evaporation': from soil [kg m-2 s-1]
                'throughfall'
                'capillar_rise'
                'pond_recharge'
                'water_closure'

                'co2_flux' [umolm m-2 (ground) s-1]
                'photosynthesis'
                'respiration'
                'soil_respiration' - from belowground

            state (dict): forestfloor aggregated state
                'temperature' [degC]
                'surface_temperature' [degC]
                'water_storage' [kg m-2]
                'snow_water_equivalent' [kg m-2]

            blt_outputs (dict of lists): bottomlayer_type fluxes and state
                'net_radiation' [W m\ :sup:`-2`\ ]
                'latent_heat' [W m\ :sup:`-2`\ ]
                'sensible_heat' [W m\ :sup:`-2`\ ]
                'ground_heat' [W m\ :sup:`-2`\ ] (negative towards soil)
                'heat_advection' [W m\ :sup:`-2`\ ]
                'water_closure' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'energy_closure' [W m\ :sup:`-2`\ ]
                'evaporation' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'interception' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'pond_recharge' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'capillary_rise' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'throughfall' [kg m\ :sup:`-2`\ s\ :sup`-1`\]

                'temperature': [\ :math:`^{\circ}`\ C]
                'volumetric_water': [m\ :sup:`3` m\ :sup:`-3`\ ]
                'water_potential': [m]
                'water_content': [g g\ :sup:`-1`\ ]
                'water_storage': [kg m\ :sup:`-2`\ ]
                'hydraulic_conductivity': [m s\ :sup:`-1`\]
                'thermal_conductivity': [W m-1 K-1]
        """
        # initialize fluxes and states
        fluxes = {
            'net_radiation': 0.0, # [W m-2]
            'sensible_heat': 0.0, # [W m-2]
            'latent_heat': 0.0, # [W m-2]
            'ground_heat': 0.0, # [W m-2]
            'energy_closure': 0.0, # [W m-2]

            'evaporation': 0.0,  # [kg m-2 s-1]
            'soil_evaporation': 0.0, # [kg m-2 s-1]
            'throughfall': 0.0,  # [kg m-2 s-1]
            'capillary_rise': 0.0,  # [kg m-2 s-1]
            'pond_recharge': 0.0, # [kg m-2 s-1]
            'water_closure': 0.0, # [kg m-2 s-1]

            'net_co2': 0.0, # [umol m-2(ground) s-1]
            'photosynthesis': 0.0,  # [umol m-2(ground) s-1]
            'respiration': 0.0,  # [umol m-2(ground) s-1]
            'soil_respiration': 0.0,  # [umol m-2(ground) s-1]
        }

        state = {
            'temperature': 0.0,  # [degC]
            'surface_temperature': 0.0,  # [degC]
            'water_storage': 0.0, # [kg m-2]
            'snow_water_equivalent': 0.0, # [kg m-2]
            # not needed as we take optical properties from previous dt
            #'albedo': None,
            #'emissivity': None
         }

        # --- Soil respiration
        fluxes['soil_respiration'] = self.soilrespiration.respiration(
                                        forcing['soil_temperature'],
                                        forcing['soil_volumetric_water'],
                                        forcing['soil_volumetric_air'])

        fluxes['respiration'] += fluxes['soil_respiration']
        fluxes['net_co2'] += fluxes['soil_respiration']

        # --- Snow: now simple degree-day model ---
        snow_forcing = {
            'precipitation_rain': forcing['precipitation_rain'],
            'precipitation_snow': forcing['precipitation_snow'],
            'air_temperature': forcing['air_temperature'],
        }

        fluxes_snow, states_snow = self.snowpack.run(dt=dt, forcing=snow_forcing)

        # --- solve bottomlayer types and aggregate forest floor fluxes & state
        org_forcing = forcing.copy()
        del org_forcing['precipitation_rain'], org_forcing['precipitation_snow']

        org_forcing.update(
                {'precipitation': fluxes_snow['potential_infiltration'],
                'soil_temperature': forcing['soil_temperature'][0],
                'snow_water_equivalent': states_snow['snow_water_equivalent']}
                )

        # bottomlayer-type specific fluxes and state for output: list of dicts
        bt_results = []

        for bt in self.bottomlayer_types:
            bt_flx, bt_state = bt.run(dt, org_forcing, parameters, controls)

            # effective forest floor fluxes and state
            for key in fluxes.keys():
                if key in bt_flx.keys():
                    fluxes[key] += bt.coverage * bt_flx[key]

            state['temperature'] += bt.coverage * bt_state['temperature']
            state['surface_temperature'] += bt.coverage * bt_state['surface_temperature']
            state['water_storage'] += bt.coverage * bt_state['water_storage']

            # merge dicts and append to gt_results
            bt_flx.update(bt_state)
            bt_results.append(bt_flx)
            del bt_flx, bt_state

        fluxes['evaporation'] += fluxes['soil_evaporation']
        fluxes['latent_heat'] += LATENT_HEAT / MOLAR_MASS_H2O * fluxes['soil_evaporation']

#            # this creates problem if forcing['air_temperature'] <0 and self.snowpack.swe >0
#        if self.snowpack.swe > 0:
#            fluxes['throughfall'] += fluxes_snow['potential_infiltration']
#
#            # some groundheat flux to keep soil temperatures reasonable
#            # this creates problem if forcing['air_temperature'] <0 and self.snowpack.swe >0
#            fluxes['ground_heat'] += (
#                parameters['soil_thermal_conductivity']
#                / abs(parameters['soil_depth'])
#                * (min(forcing['air_temperature'],0.0) - forcing['soil_temperature'][0])
#            )

        state['snow_water_equivalent'] = states_snow['snow_water_equivalent']
        state['temperature'] = states_snow['temperature']
        state['surface_temperature'] = states_snow['temperature']

        # bottomlayer_type specific results (fluxes & state): convert list of dicts to dict of lists
        blt_outputs = {}
        for k,v in bt_results[0].items():
            blt_outputs[k] = [x[k] for x in bt_results]

        return fluxes, state, blt_outputs

# EOF
