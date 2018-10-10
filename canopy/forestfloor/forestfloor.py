# -*- coding: utf-8 -*-
"""
.. module: forestfloor
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes forest floor consisting of bryophytes and/or baresoil
Based on MatLab implementation by Samuli Launiainen.


Note: roughness length z0 is [1/15 - 1/30] * x where x is height of element.
        given at the moment as a parameter, but for mosses this is problematic.
        There should be species specific 'effective roughness' (and also
        something that takes into account roughness of foresfloor)


Created on Tue Mar 13 11:59:23 2018
"""
import numpy as np
eps = np.finfo(float).eps  # machine epsilon
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

        self.baresoil = Baresoil(properties['baresoil'])

        self.snowpack = Snowpack(properties['snowpack'])

        # Bryotypes
        bryotypes = []

        f_bryo = 0.0
        for br_para in properties['bryophytes']:
            bryotypes.append(Bryophyte(br_para))
            f_bryo += br_para['ground_coverage']

        self.bryotypes = bryotypes
        # soil coverage: baresoil, bryotypes (and litter?)

        f_baresoil = 1.0 - f_bryo

        assert (f_bryo + f_baresoil <= 1.0), "Forest floor coverage more than 1.0"

        self.f_baresoil = f_baresoil
        self.baresoil.coverage = f_baresoil

        self.f_bryo = f_bryo

        if initial_temperature is not None:
            self.temperature = initial_temperature

        else:
            self.temperature = 10.

        self.old_temperature = self.temperature

    def update(self):
        """ Updates new states to the forestfloor
        """

        if self.f_bryo > 0.0:

            for bryo in self.bryotypes:
                bryo.update()

        if self.f_baresoil > 0.0:
            self.baresoil.update()

        self.old_temperature = self.temperature

    def restore(self):
        """ Restores new states back to states before iteration.
        """

        self.temperature = self.old_temperature

        for bryo in self.bryotypes:
            bryo.restore()

        self.baresoil

    def run(self, dt, forcing):
        r"""Water and energy balance at the forestfloor.

        Args:
            dt: timestep [s]
            forcing (dict): states of microclimate
                'throughfall': [mm s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ]
                'lwdn': [W m\ :sup:`-2`\ ]
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'precipitation_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_depth': [m]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [Pa]
                'soil_hydraulic_conductivity': [m s\ :sup:`-1`\ ]
                'soil_thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ]
                'snow_water_equivalent': [m]
                'nsteps' number of steps in odesolver
        Returns:
            fluxes (dict)
            states (dict)

        """
        # initialize fluxes at forest floor

        fluxes = {}
        states = {}

        soil_evaporation = 0.0  # soil evaporation (mol m-2(ground)s-1)

        sensible_heat = 0.0  # forest floor sensible heat flux (Wm-2)
        ground_heat = 0.0  # ground heat flux (Wm-2)
        latent_heat = 0.0  # latent heat flux (Wm-2)
        temperature = 0.0  # forest floor temperature (degC)

        throughfall = 0.0  # throughfall rate (m s-1)

        respiration = 0.0  # [umol m-2(ground) s-1]

        respiration, _ = soil_respiration(self.baresoil.properties,
                                          forcing['soil_temperature'],
                                          forcing['soil_volumetric_water'])

        fluxes.update({'respiration': respiration})

        # if there is snow cover forest floor is passed
        if forcing['snow_water_equivalent'] > 0:  # snow on the ground
            throughfall = forcing['throughfall']

        # calculation of forest floor heat and water balance
        else:

            # if there is bryophyte cover on the ground
            if self.f_bryo > 0.0:

                bryo_evaporation = 0.0  # evaporation from bryo [mol m-2(ground)s-1]
                capillar_water = 0.0  # water risen due to capillar forces [m s-1]
                pond_recharge = 0.0  # water moved from pond to moss [m s-1]

                bryo_temperature = 0.0  # bryophytes temperature [degC]
                bryo_water_storage = 0.0  # bryophytes water storage [kg m-2]
                bryo_water_closure = 0.0  # [W m-2]
                bryo_energy_closure = 0.0  # [W m-2]

                bryo_carbon_pool = 0.0  # [g C m-2]
                bryo_photosynthesis = 0.0  # [umol m-2 s-1]
                bryo_respiration = 0.0  # [umol m-2 s-1]

                capillar_rise = 0.0  # [m s-1]

                for bryo in self.bryotypes:

                    # bryophyte's heat, water and carbon balance
                    flxs_bryo, stts_bryo = bryo.run(dt, forcing)

                    bryo_evaporation += bryo.coverage * flxs_bryo['evaporation'] / MOLAR_MASS_H2O
                    soil_evaporation += bryo.coverage * flxs_bryo['soil_evaporation'] / MOLAR_MASS_H2O

                    throughfall += bryo.coverage * flxs_bryo['throughfall'] / WATER_DENSITY
                    capillar_rise += bryo.coverage * flxs_bryo['capillar_rise'] / WATER_DENSITY
                    pond_recharge += bryo.coverage * flxs_bryo['pond_recharge'] / WATER_DENSITY

                    latent_heat += bryo.coverage * flxs_bryo['latent_heat']
                    sensible_heat += bryo.coverage * flxs_bryo['sensible_heat']
                    ground_heat += bryo.coverage * flxs_bryo['ground_heat']

                    bryo_water_closure += bryo.coverage * flxs_bryo['water_closure']
                    bryo_energy_closure += bryo.coverage * flxs_bryo['energy_closure']

                    bryo_temperature += bryo.coverage * stts_bryo['temperature']
                    bryo_water_storage += bryo.coverage * stts_bryo['water_storage'] / WATER_DENSITY

                    # In carbon calculations are per bryophyte coverage
                    bryo_photosynthesis += flxs_bryo['photosynthesis_rate']
                    bryo_respiration += flxs_bryo['respiration_rate']
                    bryo_carbon_pool += stts_bryo['carbon_pool']

                respiration += bryo_respiration
                temperature += bryo_temperature

                fluxes.update({
                   'soil_evaporation': soil_evaporation,
                   'bryo_evaporation': bryo_evaporation,
                   'latent_heat': latent_heat,
                   'sensible_heat': sensible_heat,
                   'ground_heat': ground_heat,
                   'water_closure_bryo': bryo_water_closure,
                   'energy_closure_bryo': bryo_energy_closure,
                   'throughfall': throughfall,
                   'capillar_rise': capillar_rise,
                   'pond_recharge': pond_recharge,
                   'photosynthesis_bryo': bryo_photosynthesis,
                   'respiration_bryo': bryo_respiration,
                   'respiration': respiration
                   })

                states.update({'temperature': temperature,
                               'bryo_temperature': bryo_temperature,
                               'bryo_water_storage': bryo_water_storage,
                               'bryo_carbon_pool': bryo_carbon_pool
                               })

            if self.f_baresoil > 0.0:

                # water and energy closure
                soil_energy_closure = 0.0
                radiative_flux = 0.0  # radiative flux [W m-2]

                # soil surface energy balance
                flxs_soil, stts_soil = self.baresoil.run(dt, forcing)

                soil_evaporation += self.baresoil.coverage * flxs_soil['evaporation']
                throughfall += forcing['throughfall'] * self.baresoil.coverage

                latent_heat += self.baresoil.coverage * flxs_soil['latent_heat']
                sensible_heat += self.baresoil.coverage * flxs_soil['sensible_heat']
                ground_heat += self.baresoil.coverage * flxs_soil['ground_heat']

                soil_energy_closure += self.baresoil.coverage * flxs_soil['energy_closure']

                temperature += self.baresoil.coverage * stts_soil['temperature']
                radiative_flux += self.baresoil.coverage * flxs_soil['radiative_flux']

                fluxes.update({
                        'soil_evaporation': soil_evaporation,
                        'latent_heat': latent_heat,
                        'sensible_heat': sensible_heat,
                        'ground_heat': ground_heat,
                        'energy_closure_soil': soil_energy_closure,
                        'radiative_flux': radiative_flux,
                        'throughfall': throughfall,
                        'respiration': respiration
                        })

                states.update({'temperature': temperature,
                               'soil_temperature': stts_soil['temperature']})

        self.temperature = temperature

        return fluxes, states

    def shortwave_albedo(self, snow_water_equivalent):
        """ Forestfloor albdedo for shortwave radiation.

        Returns:
            (dict):
                'PAR_albedo': [-]
                'NIR_albedo': [-]
        """

        if snow_water_equivalent > 0.0:
            return {'PAR': 0.8,
                    'NIR': 0.8}
        else:
            par_albedo = 0.0
            nir_albedo = 0.0

            if self.f_bryo > 0.0:

                for bryo in self.bryotypes:
                    bryo_albedo = bryophyte_shortwave_albedo(
                            water_content=bryo.water_content,
                            properties=bryo.properties)

                    par_albedo += bryo.coverage * bryo_albedo['PAR']
                    nir_albedo += bryo.coverage * bryo_albedo['NIR']

            if self.f_baresoil > 0.0:
                bare_albedo = self.baresoil.albedo
                par_albedo += bare_albedo['PAR']
                nir_albedo += bare_albedo['NIR']

            return {'PAR': par_albedo, 'NIR': nir_albedo}

    def longwave_radiation(self, snow_water_equivalent, air_temperature):
        """ Forestfloor longwave radiation.

        Args:
            snow_water_equivalent: [m]
            air_temperature: [degC] in lowest calculation node
        Returns:
            (float): [W m-2]
        """

        if snow_water_equivalent > 0.0:
            return emitted_longwave_radiation(air_temperature, 0.97)

        else:
            lw_radiation = 0.0
            emissivity = 0.0

            if self.f_bryo > 0.0:

                for bryo in self.bryotypes:

                    lw_radiation += (bryo.coverage
                                     * emitted_longwave_radiation(
                                             bryo.temperature,
                                             bryo.properties))

                    emissivity += bryo.properties['optical_properties']['emissivity']

            if self.f_baresoil > 0.0:

                lw_radiation += emitted_longwave_radiation(
                        self.baresoil.temperature,
                        self.baresoil.emissivity)

                emissivity += self.f_baresoil * self.baresoil.emissivity

            return {'radiation': lw_radiation, 'emissivity': emissivity}


# EOF
