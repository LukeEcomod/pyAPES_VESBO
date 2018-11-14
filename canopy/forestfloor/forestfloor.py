#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module: forestfloor
    :synopsis: APES-model component
.. moduleauthor:: Antti-Jussi Kieloaho

Describes forest floor consisting of bryophytes, baresoil and snowpack.

Based on MatLab implementation by Samuli Launiainen.

Note:
    migrated to python3
    - absolute imports
    - dictionary.items() in for-loops

Note: roughness length z0 is [1/15 - 1/30] * x where x is height of element.
        given at the moment as a parameter, but for mosses this is problematic.
        There should be species specific 'effective roughness' (and also
        something that takes into account roughness of foresfloor)


Created on Tue Mar 13 11:59:23 2018
"""

from canopy.constants import EPS, WATER_DENSITY, MOLAR_MASS_H2O

from .bryophyte import Bryophyte
from .baresoil import Baresoil
from .snowpack import Snowpack
from .litter import Litter

from .carbon import soil_respiration
from .heat_and_water import bryophyte_shortwave_albedo, emitted_longwave_radiation


class ForestFloor(object):
    r"""Describes forest floor consisting of bryophytes and/or baresoil.
    """

    def __init__(self, properties, initial_temperature=None):
        r""" Initializes forestfloor object consisting on bryophytes
        and/or bare soil.
        Args:
            p (dict):
                'mossp' (list):
                    i. (dict): properties of bryotype i
                        'ground_coverage': fraction of moss ground coverage [-]
                        'Wmax'
                        'Wmin'
                        'zr': roughness length [m]
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

        self.snowpack = Snowpack(properties['snowpack'],
                                 properties['initial_conditions']['snowpack'])

        self.baresoil = Baresoil(properties['baresoil'],
                                 properties['initial_conditions']['baresoil'])
        self.f_baresoil = self.baresoil.coverage

        self.litter = Litter(properties['litter'],
                             properties['initial_conditions']['litter'])
        self.f_litter = self.litter.coverage

        # Bryotypes
        bryotypes = []

        f_bryo = 0.0
        for key, bryo in properties['bryophytes'].items():
            bryotypes.append(Bryophyte(bryo,
                                       properties['initial_conditions']['bryophytes'][key]))
            f_bryo += bryo['ground_coverage']

        self.bryotypes = bryotypes

        if abs(1.0 - (f_bryo + self.f_baresoil + self.f_litter)) > EPS:
            raise ValueError("The sum of bryophytes, litter and baresoil coverages "
                             + "should be one! Now %.2f" %
                             (f_bryo + self.f_baresoil + self.f_litter))

        self.f_bryo = f_bryo

        if initial_temperature is not None:
            self.temperature = initial_temperature

        else:
            self.temperature = 10.

        self.old_temperature = self.temperature

    def update(self):
        """ Updates new states to the forestfloor.
        """

        self.old_temperature = self.temperature

        if self.f_bryo > 0.0:

            for bryo in self.bryotypes:
                bryo.update()

        if self.f_litter > 0.0:

            self.litter.update()

        if self.f_baresoil > 0.0:

            self.baresoil.update()

        self.snowpack.update()

    def restore(self):
        """ Restores new states back to states before iteration.
        """

        self.temperature = self.old_temperature

        if self.f_bryo > 0.0:

            for bryo in self.bryotypes:
                bryo.restore()

        if self.f_litter > 0.0:

            self.litter.restore()

        if self.f_baresoil > 0.0:

            self.baresoil.restore()

        self.snowpack.restore()

    def run(self, dt, forcing, parameters, controls):
        r"""Water and energy balance at the forestfloor.

        Args:
            dt: timestep [s]
            forcing (dict): states of microclimate
                'throughfall_rain': [mm s\ :sup:`-1`\ ]
                'throughfall_snow': [mm s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ] if energy_balance is True
                'lw_dn': [W m\ :sup:`-2`\ ] if energy_balance is True
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'precipitation_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [Pa]
                'soil_volumetric_water': [m\ :sup:`3`\  m\`-3`\ ]
                'soil_volumetric_air': [m\ :sup:`3`\  m\`-3`\ ]
                'soil_pond_storage': [m]
                'nsteps' number of steps in odesolver
            parameters (dict):
                'soil_thermal_conductivity': [] if energy_balance is True
                'soil_hydraulic_conductivity': []
                'depth': [m] first soil calculation node
                'height': [m] first canopy calculation node
                'nsteps': number of steps in odesolver if energy_balance is True
            controls (dict):
                'energy_balance': boolean
        Returns:
            fluxes (dict)
            states (dict)

        """
        # initialize fluxes at forest floor

        fluxes = {
            'sensible_heat': 0.0,  # [W m-2]
            'ground_heat': 0.0,  # [W m-2]
            'latent_heat': 0.0,  # [W m-2]
            'potential_infiltration': 0.0,  # [m s-1]
            'respiration': 0.0,  # [umol m-2(ground) s-1]
            'soil_evaporation': 0.0,  # [mol m-2 s-1]
            'soil_energy_closure': 0.0,  # [W m-2]
            'radiative_flux': 0.0,  # [W m-2]
            'litter_evaporation': 0.0,  # [mol m-2 s-1]
            'litter_respiration': 0.0,  # [umol m-2 s-1]
            'litter_water_closure': 0.0,  # [W m-2]
            'litter_energy_closure': 0.0,  # [W m-2]
            'bryo_evaporation': 0.0,  # [mol m-2 s-1]
            'capillar_rise': 0.0,  # water movement due to capillar forces [m s-1]
            'pond_recharge': 0.0,  # water moved from pond to moss [m s-1]
            'bryo_water_closure': 0.0,  # [W m-2]
            'bryo_energy_closure': 0.0,  # [W m-2]
            'bryo_photosynthesis': 0.0,  # [umol m-2 s-1]
            'bryo_respiration': 0.0,  # [umol m-2 s-1]
        }

        states = {
            'temperature': 0.0,  # [degC]
            'soil_temperature': 0.0,  # [degC]
            'litter_temperature': 0.0,  # [degC]
            'litter_carbon_pool': 0.0,  # [g C m-2]
            'litter_water_storage': 0.0,  # [kg m-2]
            'bryo_temperature': 0.0,  # [degC]
            'bryo_water_storage': 0.0,  # [kg m-2]
            'bryo_carbon_pool': 0.0,  # [g C m-2]
        }

        # --- Snow ---

        snow_forcing = {
            'throughfall_rain': forcing['throughfall_rain'],
            'throughfall_snow': forcing['throughfall_snow'],
            'air_temperature': forcing['air_temperature'],
        }

        fluxes_snow, states_snow = self.snowpack.run(dt=dt, forcing=snow_forcing)

        forcing.update({
            'throughfall_ffloor': fluxes_snow['potential_infiltration']
        })

        # --- Soil respiration ---

        fluxes['respiration'] = soil_respiration(
            self.baresoil.properties['respiration'],
            forcing['soil_temperature'],
            forcing['soil_volumetric_water'],
            forcing['soil_volumetric_air'])

        if self.snowpack.snowcover():  # snow on the ground

            fluxes['potential_infiltration'] += fluxes_snow['potential_infiltration']
            # some groundheat flux to keep soil temperatures reasonable
            fluxes['ground_heat'] += (
                0.01 * forcing['soil_thermal_conductivity']
                / abs(forcing['depth'])
                * (forcing['air_temperature'] - forcing['soil_temperature'])
            )

            for bryo in self.bryotypes:
                states['bryo.temperature'] = 0.0
                states['bryo_water_storage'] += bryo.coverage * bryo.old_water_storage / WATER_DENSITY
                states['bryo_carbon_pool'] += bryo.old_carbon_pool

        else:
            # if there is bryophyte cover on the ground
            if self.f_bryo > 0.0:

                bryo_forcing = {
                    'throughfall': fluxes_snow['potential_infiltration'],
                    'h2o': forcing['h2o'],
                    'air_pressure': forcing['air_pressure'],
                    'par': forcing['par'],
                    'air_temperature': forcing['air_temperature'],
                    'wind_speed': forcing['wind_speed'],
                    'soil_temperature': forcing['soil_temperature'],
                    'pond_storage': forcing['soil_pond_storage'],
                }

                bryo_params = {
                    'soil_hydraulic_conductivity': parameters['soil_hydraulic_conductivity'],
                    'soil_depth': parameters['depth'],
                }

                bryo_controls = {
                    'energy_balance': controls['energy_balance'],
                    'solver': 'forward_euler',
                }

                if controls['energy_balance']:
                    bryo_forcing.update({
                        'lw_dn': forcing['lw_dn'],
                        'nir': forcing['nir']
                    })

                    bryo_params.update({
                        'soil_thermal_conductivity': parameters['soil_thermal_conductivity']
                    })

                    bryo_controls.update({
                        'nsteps': 20
                    })

                for bryo in self.bryotypes:

                    if bryo.coverage == 0.0:
                        continue

                    # bryophyte's heat, water and carbon balance
                    fluxes_bryo, states_bryo = bryo.run(
                        dt=dt,
                        forcing=bryo_forcing,
                        parameters=bryo_params,
                        controls=bryo_controls
                    )

                    fluxes['bryo_evaporation'] += bryo.coverage * fluxes_bryo['evaporation'] / MOLAR_MASS_H2O
                    fluxes['soil_evaporation'] += bryo.coverage * fluxes_bryo['soil_evaporation'] / MOLAR_MASS_H2O

                    fluxes['potential_infiltration'] += bryo.coverage * fluxes_bryo['throughfall'] / WATER_DENSITY
                    fluxes['capillar_rise'] += bryo.coverage * fluxes_bryo['capillar_rise'] / WATER_DENSITY
                    fluxes['pond_recharge'] += bryo.coverage * fluxes_bryo['pond_recharge'] / WATER_DENSITY

                    fluxes['latent_heat'] += bryo.coverage * fluxes_bryo['latent_heat']
                    fluxes['sensible_heat'] += bryo.coverage * fluxes_bryo['sensible_heat']
                    fluxes['ground_heat'] += bryo.coverage * fluxes_bryo['ground_heat']

                    fluxes['bryo_water_closure'] += bryo.coverage * fluxes_bryo['water_closure']
                    fluxes['bryo_energy_closure'] += bryo.coverage * fluxes_bryo['energy_closure']

                    states['bryo_temperature'] += bryo.coverage * states_bryo['temperature']
                    states['bryo_water_storage'] += bryo.coverage * states_bryo['water_storage'] / WATER_DENSITY

                    # In carbon calculations are per bryophyte's coverage
                    fluxes['bryo_photosynthesis'] += fluxes_bryo['photosynthesis_rate']
                    fluxes['bryo_respiration'] += fluxes_bryo['respiration_rate']
                    states['bryo_carbon_pool'] += states_bryo['carbon_pool']

                fluxes['respiration'] += fluxes['bryo_respiration']
                states['temperature'] += states['bryo_temperature']

                states['bryo_temperature'] = states['bryo_temperature'] / self.f_bryo

            if self.f_litter > 0.0:

                litter_forcing = {
                    'throughfall': fluxes_snow['potential_infiltration'],
                    'h2o': forcing['h2o'],
                    'air_pressure': forcing['air_pressure'],
                    'par': forcing['par'],
                    'air_temperature': forcing['air_temperature'],
                    'wind_speed': forcing['wind_speed'],
                    'soil_temperature': forcing['soil_temperature'],
                    'pond_storage': forcing['soil_pond_storage'],
                }

                litter_params = {
                    'soil_hydraulic_conductivity': 0.0,
                    'soil_depth': parameters['depth'],
                }

                litter_controls = {
                    'energy_balance': controls['energy_balance'],
                    'solver': 'forward_euler',
                }

                if controls['energy_balance']:
                    litter_forcing.update({
                        'lw_dn': forcing['lw_dn'],
                        'nir': forcing['nir']
                    })

                    litter_params.update({
                        'soil_thermal_conductivity': parameters['soil_thermal_conductivity']
                    })

                    litter_controls.update({
                        'nsteps': 20
                    })

                # litters's heat, water and respiration
                fluxes_litter, states_litter = self.litter.run(
                    dt=dt,
                    forcing=litter_forcing,
                    parameters=litter_params,
                    controls=litter_controls
                )

                fluxes['litter_evaporation'] += self.f_litter * fluxes_litter['evaporation'] / MOLAR_MASS_H2O
                fluxes['soil_evaporation'] += self.f_litter * fluxes_litter['soil_evaporation'] / MOLAR_MASS_H2O

                fluxes['potential_infiltration'] += self.f_litter * fluxes_litter['throughfall'] / WATER_DENSITY
                fluxes['capillar_rise'] += self.f_litter * fluxes_litter['capillar_rise'] / WATER_DENSITY
                fluxes['pond_recharge'] += self.f_litter * fluxes_litter['pond_recharge'] / WATER_DENSITY

                fluxes['latent_heat'] += self.f_litter * fluxes_litter['latent_heat']
                fluxes['sensible_heat'] += self.f_litter * fluxes_litter['sensible_heat']
                fluxes['ground_heat'] += self.f_litter * fluxes_litter['ground_heat']

                fluxes['litter_water_closure'] += self.f_litter * fluxes_litter['water_closure']
                fluxes['litter_energy_closure'] += self.f_litter * fluxes_litter['energy_closure']

                states['litter_temperature'] = states_litter['temperature']
                states['temperature'] += self.f_litter * states['litter_temperature']
                states['litter_water_storage'] = self.f_litter * states_litter['water_storage'] / WATER_DENSITY

                # In carbon calculations are per litter's coverage
                fluxes['litter_respiration'] += fluxes_litter['respiration_rate']
                states['litter_carbon_pool'] += states_litter['carbon_pool']
                fluxes['respiration'] += fluxes['litter_respiration']

            if self.f_baresoil > 0.0:

                bare_forcing = {
                    'wind_speed': forcing['wind_speed'],
                    'air_temperature': forcing['air_temperature'],
                    'h2o': forcing['h2o'],
                    'air_pressure': forcing['air_pressure'],
                    'forestfloor_temperature': self.temperature,
                    'soil_temperature': forcing['soil_temperature'],
                    'soil_water_potential': forcing['soil_water_potential'],
                    'par': forcing['par'],
                    'nir': forcing['nir'],
                    'lw_dn': forcing['lw_dn'],
                    'lw_up': forcing['lw_up'],
                }

                bare_params = {
                    'soil_hydraulic_conductivity': parameters['soil_hydraulic_conductivity'],
                    'soil_thermal_conductivity': parameters['soil_thermal_conductivity'],
                    'depth': parameters['depth'],
                    'height': parameters['height']
                }

                bare_controls = {
                    'energy_balance': controls['energy_balance']
                }

                fluxes_soil, states_soil = self.baresoil.run(
                    dt=dt,
                    forcing=bare_forcing,
                    parameters=bare_params,
                    controls=bare_controls
                )

                fluxes['soil_evaporation'] += self.f_baresoil * fluxes_soil['evaporation']

                fluxes['latent_heat'] += self.f_baresoil * fluxes_soil['latent_heat']
                fluxes['sensible_heat'] += self.f_baresoil * fluxes_soil['sensible_heat']
                fluxes['ground_heat'] += self.f_baresoil * fluxes_soil['ground_heat']

                fluxes['radiative_flux'] += self.f_baresoil * fluxes_soil['radiative_flux']
                fluxes['soil_energy_closure'] += self.f_baresoil * fluxes_soil['energy_closure']

                fluxes['potential_infiltration'] += self.f_baresoil * forcing['throughfall_ffloor']

                states['soil_temperature'] = states_soil['temperature']
                states['temperature'] += self.f_baresoil * states['soil_temperature']

        self.temperature = states['temperature']
        fluxes['evaporation'] = fluxes['soil_evaporation'] + fluxes['bryo_evaporation'] + fluxes['litter_evaporation']
        fluxes['water_closure_snow'] = fluxes_snow['water_closure']

        states.update(states_snow)

        return fluxes, states

    def shortwave_albedo(self):
        """ Forestfloor albdedo for shortwave radiation.

        Returns:
            (dict):
                'PAR_albedo': [-]
                'NIR_albedo': [-]
        """

        if self.snowpack.snowcover():
            par_albedo = self.snowpack.properties['optical_properties']['albedo_PAR']
            nir_albedo = self.snowpack.properties['optical_properties']['albedo_NIR']

        else:
            par_albedo = 0.0
            nir_albedo = 0.0

            if self.f_bryo > 0.0:

                for bryo in self.bryotypes:
                    if bryo.coverage == 0.0:
                        continue

                    bryo_albedo = bryophyte_shortwave_albedo(
                        water_content=bryo.water_content,
                        properties=bryo.properties)

                    par_albedo += bryo.coverage * bryo_albedo['PAR']
                    nir_albedo += bryo.coverage * bryo_albedo['NIR']

            if self.f_baresoil > 0.0:

                bare_optical = self.baresoil.properties['optical_properties']
                par_albedo += self.baresoil.coverage * bare_optical['albedo_PAR']
                nir_albedo += self.baresoil.coverage * bare_optical['albedo_NIR']

        return {'PAR': par_albedo, 'NIR': nir_albedo}

    def longwave_radiation(self):
        """ Forestfloor longwave radiation.

        Args:
            snow_water_equivalent: [m]
            air_temperature: [degC] in lowest calculation node
        Returns:
            (float): [W m-2]
        """

        if self.snowpack.snowcover():

            lw_radiation = emitted_longwave_radiation(self.temperature,
                                                      self.snowpack.properties)

            emissivity = self.snowpack.properties['optical_properties']['emissivity']

        else:
            lw_radiation = 0.0
            emissivity = 0.0

            if self.f_bryo > 0.0:

                for bryo in self.bryotypes:

                    if bryo.coverage == 0.0:
                        continue

                    lw_radiation += (
                            bryo.coverage
                            * emitted_longwave_radiation(
                                bryo.temperature,
                                bryo.properties))

                    emissivity += (bryo.coverage
                                   * bryo.properties['optical_properties']['emissivity'])

            if self.f_baresoil > 0.0:

                lw_radiation += (self.baresoil.coverage
                                 * emitted_longwave_radiation(
                                         self.baresoil.temperature,
                                         self.baresoil.properties))

                emissivity += (self.baresoil.coverage
                               * self.baresoil.properties['optical_properties']['emissivity'])

        return {'radiation': lw_radiation, 'emissivity': emissivity}


# EOF
