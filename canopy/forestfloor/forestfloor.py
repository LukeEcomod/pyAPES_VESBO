# -*- coding: utf-8 -*-
"""
.. module: forestfloor
    :synopsis: APES-model component
.. moduleauthor:: Antti-Jussi Kieloaho

Describes forest floor consisting of bryophytes, baresoil and snowpack.

Based on MatLab implementation by Samuli Launiainen.


Note: roughness length z0 is [1/15 - 1/30] * x where x is height of element.
        given at the moment as a parameter, but for mosses this is problematic.
        There should be species specific 'effective roughness' (and also
        something that takes into account roughness of foresfloor)


Created on Tue Mar 13 11:59:23 2018
"""
import numpy as np

from canopy.micromet import e_sat
from canopy.constants import EPS, WATER_DENSITY, MOLAR_MASS_H2O

from bryophyte import Bryophyte
from baresoil import Baresoil
from snowpack import Snowpack

from carbon import soil_respiration
from heat_and_water import soil_boundary_layer_conductance
from heat_and_water import bryophyte_shortwave_albedo, emitted_longwave_radiation

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

        self.baresoil = Baresoil(properties['baresoil'],
                                 properties['initial_conditions']['baresoil'])

        self.snowpack = Snowpack(properties['snowpack'],
                                 properties['initial_conditions']['snowpack'])

        # Bryotypes
        bryotypes = []

        f_bryo = 0.0
        for bryo in properties['bryophytes']:
            bryotypes.append(Bryophyte(bryo,
                                       properties['initial_conditions']['bryophytes'][bryo['species']]))
            f_bryo += bryo['ground_coverage']

        self.bryotypes = bryotypes
        # soil coverage: baresoil, bryotypes (and litter?)

        f_baresoil = 1.0 - f_bryo

        if f_bryo + f_baresoil > 1.0:
            raise ValueError("The sum of bryophytes and baresoil coverages "
                             + "more than one!")

        self.f_baresoil = f_baresoil
        self.baresoil.coverage = f_baresoil

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

        if self.f_baresoil > 0.0:

            self.baresoil.restore()

        self.snowpack.restore()

    def run(self, dt, forcing):
        r"""Water and energy balance at the forestfloor.

        Args:
            dt: timestep [s]
            forcing (dict): states of microclimate
                'throughfall_rain': [mm s\ :sup:`-1`\ ]
                'throughfall_snow': [mm s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ]
                'lwdn': [W m\ :sup:`-2`\ ]
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'precipitation_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'depth': [m]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [Pa]
                'soil_hydraulic_conductivity': [m s\ :sup:`-1`\ ]
                'soil_thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ]
                'nsteps' number of steps in odesolver
        Returns:
            fluxes (dict)
            states (dict)

        """
        # initialize fluxes at forest floor

        fluxes = {}
        states = {}

        # --- Forestfloor ---

        sensible_heat = 0.0  # [W m-2]
        ground_heat = 0.0  # [W m-2]
        latent_heat = 0.0  # [W m-2]
        temperature = 0.0  # [degC]
        evaporation = 0.0  # [mol m-2(ground) s-1]

        potential_infiltration = 0.0  # [m s-1]
        respiration = 0.0  # [umol m-2(ground) s-1]

        # --- Soil ---

        soil_evaporation = 0.0  # [mol m-2 s-1]

        # baresoil
        soil_temperature = 0.0  # [degC]
        soil_energy_closure = 0.0  # [W m-2]
        radiative_flux = 0.0  # [W m-2]

        # --- Bryphytes ---

        bryo_evaporation = 0.0  # [mol m-2 s-1]
        capillar_rise = 0.0  # water movement due to capillar forces [m s-1]
        pond_recharge = 0.0  # water moved from pond to moss [m s-1]

        bryo_temperature = 0.0  # [degC]
        bryo_water_storage = 0.0  # [kg m-2]
        bryo_water_closure = 0.0  # [W m-2]
        bryo_energy_closure = 0.0  # [W m-2]

        bryo_carbon_pool = 0.0  # [g C m-2]
        bryo_photosynthesis = 0.0  # [umol m-2 s-1]
        bryo_respiration = 0.0  # [umol m-2 s-1]

        # --- Snow ---

        fluxes_snow, states_snow = self.snowpack.run(dt=dt, forcing=forcing)
        forcing.update({'throughfall_ffloor': fluxes_snow['potential_infiltration']})

        # --- Soil respiration ---

        respiration, _ = soil_respiration(self.baresoil.properties,
                                          forcing['soil_temperature'],
                                          forcing['soil_volumetric_water'])

        if self.snowpack.snowcover():  # snow on the ground

            potential_infiltration += fluxes_snow['potential_infiltration']
            # some groundheat flux to keep soil temperatures reasonable
            ground_heat += 0.01 * forcing['soil_thermal_conductivity'] / abs(forcing['depth']) * (
                           forcing['air_temperature'] - forcing['soil_temperature'])

            for bryo in self.bryotypes:
                # if something goes wrong in snow melt check this!
                bryo.temperature = 0.0
                bryo_water_storage += bryo.water_storage
                bryo_carbon_pool += bryo.carbon_pool

        else:

            # if there is bryophyte cover on the ground
            if self.f_bryo > 0.0:

                for bryo in self.bryotypes:

                    if bryo.coverage == 0.0:
                        continue

                    # bryophyte's heat, water and carbon balance
                    fluxes_bryo, states_bryo = bryo.run(dt, forcing)

                    bryo_evaporation += bryo.coverage * fluxes_bryo['evaporation'] / MOLAR_MASS_H2O
                    soil_evaporation += bryo.coverage * fluxes_bryo['soil_evaporation'] / MOLAR_MASS_H2O

                    potential_infiltration += bryo.coverage * fluxes_bryo['throughfall'] / WATER_DENSITY
                    capillar_rise += bryo.coverage * fluxes_bryo['capillar_rise'] / WATER_DENSITY
                    pond_recharge += bryo.coverage * fluxes_bryo['pond_recharge'] / WATER_DENSITY

                    latent_heat += bryo.coverage * fluxes_bryo['latent_heat']
                    sensible_heat += bryo.coverage * fluxes_bryo['sensible_heat']
                    ground_heat += bryo.coverage * fluxes_bryo['ground_heat']

                    bryo_water_closure += bryo.coverage * fluxes_bryo['water_closure']
                    bryo_energy_closure += bryo.coverage * fluxes_bryo['energy_closure']

                    bryo_temperature += bryo.coverage * states_bryo['temperature']
                    bryo_water_storage += bryo.coverage * states_bryo['water_storage'] / WATER_DENSITY

                    # In carbon calculations are per bryophyte's coverage
                    bryo_photosynthesis += fluxes_bryo['photosynthesis_rate']
                    bryo_respiration += fluxes_bryo['respiration_rate']
                    bryo_carbon_pool += states_bryo['carbon_pool']

                respiration += bryo_respiration
                temperature += bryo_temperature
                evaporation += bryo_evaporation

            if self.f_baresoil > 0.0:

                # baresoil surface energy balance
                forcing.update({'forestfloor_temperature': self.temperature})
                fluxes_soil, states_soil = self.baresoil.run(dt, forcing)

                soil_evaporation += self.baresoil.coverage * fluxes_soil['evaporation']
                evaporation += soil_evaporation

                latent_heat += self.baresoil.coverage * fluxes_soil['latent_heat']
                sensible_heat += self.baresoil.coverage * fluxes_soil['sensible_heat']
                ground_heat += self.baresoil.coverage * fluxes_soil['ground_heat']

                radiative_flux += self.baresoil.coverage * fluxes_soil['radiative_flux']
                soil_energy_closure += self.baresoil.coverage * fluxes_soil['energy_closure']

                potential_infiltration += self.baresoil.coverage * forcing['throughfall_ffloor']

                soil_temperature += states_soil['temperature']
                temperature += self.baresoil.coverage * states_soil['temperature']

        self.temperature = temperature

        fluxes.update({
                'evaporation': evaporation,
                'soil_evaporation': soil_evaporation,
                'bryo_evaporation': bryo_evaporation,
                'latent_heat_flux': latent_heat,
                'sensible_heat_flux': sensible_heat,
                'ground_heat_flux': ground_heat,
                'water_closure_bryo': bryo_water_closure,
                'water_closure_snow': fluxes_snow['water_closure'],
                'energy_closure_bryo': bryo_energy_closure,
                'energy_closure_soil': soil_energy_closure,
                'radiative_flux': radiative_flux,
                'potential_infiltration': potential_infiltration,
                'capillar_rise': capillar_rise,
                'pond_recharge': pond_recharge,
                'photosynthesis_bryo': bryo_photosynthesis,
                'respiration_bryo': bryo_respiration,
                'respiration': respiration,
                })

        states.update(states_snow)
        states.update({
                'temperature': temperature,
                'bryo_temperature': bryo_temperature,
                'soil_temperature': soil_temperature,
                'bryo_water_storage': bryo_water_storage,
                'bryo_carbon_pool': bryo_carbon_pool
                })

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
