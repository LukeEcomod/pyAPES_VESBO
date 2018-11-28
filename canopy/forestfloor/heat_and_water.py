#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:23:08 2018

Note:
    migrated to python3
    - absolute import

@author: ajkieloaho
"""

import numpy as np
from scipy.integrate import odeint

from canopy.constants import EPS, MOLAR_MASS_H2O, SPECIFIC_HEAT_H2O
from canopy.constants import SPECIFIC_HEAT_ORGANIC_MATTER, LATENT_HEAT
from canopy.constants import DEG_TO_KELVIN, STEFAN_BOLTZMANN
from canopy.constants import SPECIFIC_HEAT_AIR, WATER_DENSITY, GAS_CONSTANT
from canopy.constants import MOLECULAR_DIFFUSIVITY_CO2, MOLECULAR_DIFFUSIVITY_H2O
from canopy.constants import THERMAL_DIFFUSIVITY_AIR, GRAVITY
from canopy.constants import AIR_VISCOSITY, AIR_DENSITY
from .odesolver import solver_array, ForwardEuler_array
#from canopy.constants import *


class heat_and_water:
    """ Defines conservation of heat and water balance of bryotype object

    Note:
        Class is callable so it can be used with scipy odeint and odesolver
    """

    def __init__(self, properties, states):
        r""" Inititalises a heat and water balance calculation object

        Args:
            properties (dict): BryoType instance properties

            states (dict): forcing states of bryophyte heat and water balance
                'soil_thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ]
                'partial_pressure_h2o': [mol mol\ :sup:`-1`\ ]
                'air_pressure': [Pa]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ]
                'lwdn': [W m\ :sup:`-2`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'wind_speed': [m s\ :sup:`-1`\ ]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [m]
                'soil_hydraulic_conductivity': [m s\ :sup:`-1`\ ]
                'soil_depth': [m]
                'precipitation': [mm s\ :sup:`-1`\ ]
                'precipitation_temperature': [\ :math:`^{\circ}`\ C]
                'pond_storage': [mm s-1]
        """

        self.dudt = np.zeros(12)

        self.properties = properties
        self.states = states

        self.min_water = (
            properties['min_water_content'] * properties['dry_mass'])

        self.max_water = (
            properties['max_water_content'] * properties['dry_mass'])

    def __call__(self, y, dt):
        """
        Args:
            y (array): initial values
                0. temperature (du/dt)
                1. water content (du/dt)
                2. pond_storage
                3. capillar_water
                4. throughfall
            dt (float): time step
        """

        if dt == 0.0:
            dt = dt + EPS

        # [kg m-2] or [mm]

        water = min(self.max_water, y[1])
        water = max(self.min_water, water)

        max_recharge = max(self.max_water - water, 0.0)
        max_recharge = min(self.max_water, max_recharge)

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max_recharge / dt

        # [kg m-2 s-1] or [mm s-1]
        max_evaporation_rate = (y[1] - (self.min_water + EPS)) / dt

        if np.isinf(max_evaporation_rate) or max_evaporation_rate < 0.0:
            max_evaporation_rate = 0.0

        # [kg m-2 s-1] or [mm s-1]
        max_condensation_rate = -((self.max_water - EPS) - y[1]) / dt

        if np.isinf(max_condensation_rate) or max_condensation_rate > 0.0:
            max_condensation_rate = 0.0

        # boundary layer conductances for H2O, heat and CO2

        temp_difference = y[0] - self.states['air_temperature']

        # [mol m-2 s-1]
        conductance_to_air = moss_atm_conductance(self.states['wind_speed'],
                                                  self.properties['roughness_height'],
                                                  dT=temp_difference)

        # water vapor conductance from moss to air
        # Relative conductance is from Williams and Flanagan (1996), Oecologia.
        # Crosses 1 at appr. when water content/dry mass is 8.8
        relative_conductance = min(1.0,
                                   (0.1285 * y[1]
                                    / self.properties['dry_mass'] - 0.1285))

        # [mol m-2 s-1]
        conductance_to_air_h2o = (
            conductance_to_air['h2o'] * relative_conductance)

        # [kg m-2 s-1]
        evaporation_rate = (conductance_to_air_h2o
                            * (611.0 / self.states['air_pressure']
                                * np.exp(17.502 * y[0] / (y[0] + 240.0))
                                - self.states['h2o'])
                            * MOLAR_MASS_H2O)

        evaporation_rate = min(evaporation_rate, max_evaporation_rate)
        evaporation_rate = max(evaporation_rate, max_condensation_rate)

        max_recharge_rate = max(max_recharge_rate + evaporation_rate, 0.0)
        max_recharge = max_recharge_rate * dt

        # take into account that maximum water capacity is not fully achieved

        # [mm] or [kg m-2]
        interception = (max_recharge *
                        (1.0 - np.exp(-(1.0 / self.max_water)
                         * self.states['precipitation'] * dt)))

        # [kg m-2 s-1] or [mm s-1]
        interception_rate = interception / dt

        # [kg m-2 s-1] or [mm s-1]
        max_recharge = max(max_recharge - interception,
                           0.0)

        pond_recharge = (max_recharge *
                         (1.0 - np.exp(-(1.0 / self.max_water)
                          * self.states['soil_pond_storage'] * dt)))

        # [kg m-2 s-1] or [mm s-1]
        pond_recharge_rate = pond_recharge / dt

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max(max_recharge - pond_recharge, 0.0) / dt

        # [g g-1]
        water_content = (water / self.properties['dry_mass'])

        # [m m-3]
        volumetric_water = (water_content / WATER_DENSITY
                            * self.properties['bulk_density'])

        # [m]
        water_potential = convert_hydraulic_parameters(
                volumetric_water,
                self.properties['water_retention'],
                'volumetric_water')

        # [m s-1]
        hydraulic_conductivity = hydraulic_conduction(
                water_potential,
                self.properties['water_retention'])

        # [kg m-2 s-1] or [mm s-1]
        capillary_rise = capillarity(dt=dt,
                                     properties=self.properties,
                                     hydraulic_conductivity=hydraulic_conductivity,
                                     water_potential=water_potential,
                                     water_content=water_content,
                                     soil_hydraulic_conductivity=self.states['soil_hydraulic_conductivity'],
                                     soil_water_potential=self.states['soil_water_potential'],
                                     soil_depth=self.states['soil_depth'])

        # [kg m-2 s-1] or [mm s-1]
        capillary_rise = min(capillary_rise, max_recharge_rate)

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max(max_recharge_rate - capillary_rise, 0.0)

        # calculate mass balance of water

        # [kg m-2 s-1] or [mm s-1]
        dy_water = (
                interception_rate
                + pond_recharge_rate
                + capillary_rise
                - evaporation_rate
                )

        albedo = bryophyte_shortwave_albedo(
                water_content,
                self.properties)

        emissivity = self.properties['optical_properties']['emissivity']

        # [J m-2 s-1] or [W m-2]
        total_absorbed_radiation = (
                self.states['par'] * (1.0 - albedo['PAR'])
                + self.states['nir'] * (1.0 - albedo['NIR'])
                + emissivity * self.states['lw_dn'])

        # [J m-2 s-1] or [W m-2]
        net_radiation = (
            total_absorbed_radiation
            - emissivity * STEFAN_BOLTZMANN
            * (y[0] + DEG_TO_KELVIN)**4.0)

        # [J m-2 s-1] or [W m-2]
        sensible_heat = (
            SPECIFIC_HEAT_AIR
            * conductance_to_air['heat']
            * (y[0] - self.states['air_temperature']))

        # [J m-2 s-1] or [W m-2]
        latent_heat = LATENT_HEAT / MOLAR_MASS_H2O * evaporation_rate + EPS

        # thermal conductance between moss and soil

        # [W m-2 K-1]
        moss_thermal_conductivity = thermal_conduction(volumetric_water)

        # geometric mean
        # thermal conductivity x depth. depth = soildepth + 0.5 x moss height
        thermal_cond_to_soil = (
                np.power(moss_thermal_conductivity
                         * self.states['soil_thermal_conductivity'], 0.5)
                / (abs(self.states['soil_depth']) + 0.5 * self.properties['height']))

#        # two resistors in series
#        moss_resistance = moss_thermal_conductivity / (0.5 * self.properties['height'])
#        soil_resistance = self.states['soil_thermal_conductivity'] / abs(self.states['soil_depth'])
#        thermal_cond_to_soil = (moss_resistance * soil_resistance
#                                / (moss_resistance + soil_resistance))

#        # thermal conductivity as arithmetic average
#        k = 0.5 * (moss_thermal_conductivity + self.states['soil_thermal_conductivity'])
#        thermal_conductance3 = k / (abs(self.states['soil_depth']) + 0.5 * self.properties['height'])

        # [J m-2 s-1] or [W m-2]
        heat_conduction = (
                thermal_cond_to_soil
                * (y[0] - self.states['soil_temperature']))

        # internal heat lost or gained with water removing/entering

        # [J m-2 s-1] or [W m-2]
        heat_advection = SPECIFIC_HEAT_H2O * (
                        interception_rate * self.states['precipitation_temperature']
                        + capillary_rise * self.states['soil_temperature']
                        + pond_recharge_rate * self.states['soil_temperature']
                        )

        # heat capacities

        # [J K-1]
        heat_capacity_old = (SPECIFIC_HEAT_ORGANIC_MATTER
                         * self.properties['dry_mass']
                         + SPECIFIC_HEAT_H2O * y[1])

        heat_capacity_new = (SPECIFIC_HEAT_ORGANIC_MATTER
                             * self.properties['dry_mass']
                             + SPECIFIC_HEAT_H2O * (y[1] + dy_water * dt))


        # calculate heat balance

        heat_fluxes = (
                net_radiation
                - sensible_heat
                - latent_heat
                - heat_conduction
                + heat_advection
                )

        new_temperature = (heat_fluxes * dt + heat_capacity_old * y[0]) / heat_capacity_new
        dy_temperature = (new_temperature - y[0]) / dt

        # water components
        # [kg m-2 s-1] or [mm s-1]
        self.dudt[2] = pond_recharge_rate
        self.dudt[3] = capillary_rise
        self.dudt[4] = interception_rate
        self.dudt[5] = evaporation_rate

        # energy components
        # [J m-2 s-1] or [W m-2]
        self.dudt[6] = net_radiation
        self.dudt[7] = sensible_heat
        self.dudt[8] = latent_heat
        self.dudt[9] = heat_conduction
        self.dudt[10] = heat_advection
        self.dudt[11] = heat_fluxes

        # [K m-2 s-1]
        self.dudt[0] = dy_temperature
        # [kg m-2 s-1] or [mm s-1]
        self.dudt[1] = dy_water

        return self.dudt

def heat_and_water_exchange(properties,
                            temperature,
                            water_content,
                            dt,
                            forcing,
                            parameters,
                            solver='forward_euler',
                            nsteps=20):
    r""" Estimates energy and water balances of a bulk bryosphyte layer.

    Takes into account for soil-moss-air interactions.

    Args:
        properties (dict): characteristics of Bryophyte instance
        moss_temperature: [\ :math:`^{\circ}`\ C]
        moss_water_content: [g g\ :sup:`-1`\ ]
        dt: [s], timestep
        forcing (dict):
            'throughfall': [mm s\ :sup:`-1`\ ]
            'par': [W m\ :sup:`-2`\ ]
            'nir': [W m\ :sup:`-2`\ ]
            'lw_dn': [W m\ :sup:`-2`\ ]
            'h2o': [mol mol\ :sup:`-1`\ ]
            'air_temperature': [\ :math:`^{\circ}`\ C]
            'precipitation_temperature': [\ :math:`^{\circ}`\ C]
            'air_pressure': [Pa]
            'soil_depth': [m]
            'soil_temperature': [\ :math:`^{\circ}`\ C]
            'soil_water_potential': [Pa]
        parameters (dict):
            'soil_hydraulic_conductivity'
            'soil_thermal_conductivity'
        solver (str):
            'forward_euler': default
            'odeint': scipy.odeint
        nsteps (int): number of steps in odesolver

    Returns:
        fluxes (dict):
            'net_radiation_balance': [W m\ :sup:`-2`\ ]
            'latent_heat_flux': [W m\ :sup:`-2`\ ]
            'sensible_heat_flux': [W m\ :sup:`-2`\ ]
            'ground_heat_flux': [W m\ :sup:`-2`\ ] (negative towards soil)
            'emitted_longwave_radiation': [W m\ :sup:`-2`\ ]
            'heat_storage_change': [W m\ :sup:`-2`\ ]
            'interception': [kg m\ :sup:`-2`\ ]
            'throughfall_rate': [kg m\ :sup:`-2` s\ :sup:`-1`\ ]
            'capillary_rise': [mm s\ :sup:`-1`\ ]
        states (dict):
            'temperature': [\ :math:`^{\circ}`\ C]
            'volumetric_water': [m\ :sup:`3` m\ :sup:`-3`\ ]
            'water_potential': [m]
            'water_content': [g g\ :sup:`-1`\ ]
            'water_storage_change': [kg m\ :sup:`-2`\ ]
            'hydraulic_conductivity': [ms-1]
            'thermal_conductivity': [Wm-1K-1]
            'pond recharge': [m]
    """
    solver = solver.lower()

    precipitation = forcing['throughfall'] * WATER_DENSITY  # [m s-1] -> [mm s-1]

    # store previous dt value for computing water closure;
    # noticed later this is as input moss_water_content
    # water content from previous time step

    states = {
        'h2o': forcing['h2o'],
        'air_pressure': forcing['air_pressure'],
        'par': forcing['par'],
        'nir': forcing['nir'],
        'lw_dn': forcing['lw_dn'],
        'air_temperature': forcing['air_temperature'],  # [deg C]
        'wind_speed': forcing['wind_speed'],
        'soil_temperature': forcing['soil_temperature'],  # [deg C]
        'soil_water_potential': forcing['soil_water_potential'],
        'precipitation': precipitation,
        'precipitation_temperature': forcing['air_temperature'],
        'soil_pond_storage': forcing['soil_pond_storage'] * 1000.0 / dt,  # [mm s-1]
        'soil_hydraulic_conductivity': parameters['soil_hydraulic_conductivity'],
        'soil_thermal_conductivity': parameters['soil_hydraulic_conductivity'],
        'soil_depth': parameters['soil_depth'],
    }

    # initializing state equation
    state_equation = heat_and_water(properties, states)

    # initial values for solver

    # [deg C]
    initial_temperature = temperature
    # [kg m-2] or [mm]
    initial_water_content = water_content * properties['dry_mass']

    # initial values for water components

    # [mm]
    initial_pond_recharge = 0.0
    # [kg m-2] or [mm]
    initial_capillar_water = 0.0
    # [mm]
    initial_interception = 0.0
    # [mm]
    initial_evaporation = 0.0

    # initial values for energy components

    # [J m-2]
    initial_radiative_heat = 0.0
    # [J m-2]
    initial_sensible_heat = 0.0
    # [J m-2]
    initial_latent_heat = 0.0
    # [J m-2]
    initial_conducted_heat = 0.0
    # [J m-2]
    initial_advected_heat = 0.0
    # [J m-2]
    initial_energy_closure = 0.0

    tt = np.linspace(0.0, dt, nsteps)

    if solver == 'forward_euler':

        new_states, time_points = solver_array(
                func=state_equation,
                init=[initial_temperature,
                      initial_water_content,
                      initial_pond_recharge,
                      initial_capillar_water,
                      initial_interception,
                      initial_evaporation,
                      initial_radiative_heat,
                      initial_sensible_heat,
                      initial_latent_heat,
                      initial_conducted_heat,
                      initial_advected_heat,
                      initial_energy_closure],
                time=tt,
                method=ForwardEuler_array)

    elif solver == 'odeint':

        new_states = odeint(
                func=state_equation,
                y0=[initial_temperature,
                    initial_water_content,
                    initial_pond_recharge,
                    initial_capillar_water,
                    initial_interception,
                    initial_evaporation,
                    initial_radiative_heat,
                    initial_sensible_heat,
                    initial_latent_heat,
                    initial_conducted_heat,
                    initial_advected_heat,
                    initial_energy_closure],
                t=tt,
                full_output=False,
                mxstep=3600)

    del tt

    if np.any(np.isnan(new_states)):
        print('AAAAARGH!!! temp: {}, waterc: {}, evap: {}, h2o{}'.format(
                new_states[0][-1],
                new_states[1][-1],
                new_states[5][-1],
                forcing['h2o']
                ))

    # new states in every solver timestep
    # [deg C]
    new_temperature = new_states[:, 0]
    # [kg m-2]
    new_water_content = new_states[:, 1]

    # heat capacity in solver timesteps
    # [J m-2]
    heat_storage = ((SPECIFIC_HEAT_ORGANIC_MATTER * properties['dry_mass']
                     + SPECIFIC_HEAT_H2O * new_water_content)
                    * new_temperature)

    # heat storage change in solver timesteps
    # calculated change from previous timestep
    # [W m-2] or [J m-2 s-1]

    heat_storage_change = (heat_storage[-1] - heat_storage[0]) / dt

    net_radiation = new_states[-1, 6] / dt
    sensible_heat_flux = new_states[-1, 7] / dt
    latent_heat_flux = new_states[-1, 8] / dt
    ground_heat_flux = new_states[-1, 9] / dt
    heat_advection = new_states[-1, 10] / dt

    # [W m-2] or [J m-2 s-1]
    heat_fluxes = (net_radiation
                   - sensible_heat_flux
                   - latent_heat_flux
                   - ground_heat_flux
                   + heat_advection)

    energy_closure = -heat_storage_change + heat_fluxes

    # water storage change in solver timesteps
    # calculated change from previous timestep
    # [kg m-2 s-1] or [mm s-1]
    water_storage_change = (new_water_content[-1] - new_water_content[0]) / dt

    # water fluxes
    # [kg m-2 s-1] or [mm s-1]
    pond_recharge_rate = new_states[-1, 2] / dt
    capillary_rise = new_states[-1, 3] / dt
    interception_rate = new_states[-1, 4] / dt
    evaporation_rate = new_states[-1, 5] / dt

    # [kg m-2 s-1] or [mm s-1]
    water_fluxes = (pond_recharge_rate
                    + capillary_rise
                    - evaporation_rate
                    + interception_rate)

    water_closure = -water_storage_change + water_fluxes

    # last solution is new state
    new_temperature = new_temperature[-1]
    new_water_content = new_water_content[-1]
    pond_recharge = new_states[-1, 2]

    # [mm]
    new_capillar_water = new_states[-1, 3]

    # [kg m-2 s-1)] of [mm s-1]
    throughfall = (precipitation - interception_rate) * dt

    # State variables of bryophyte layer
    # [g g-1]
    new_gravimetric_water = (
        new_water_content / properties['dry_mass'])

    # [m3 m-3]
    new_volumetric_water = (
        new_water_content
        / WATER_DENSITY
        * properties['bulk_density'])

    # [m]
    new_matrix_potential = convert_hydraulic_parameters(
        new_volumetric_water,
        properties['water_retention'],
        'volumetric_water')

    # [m s-1]
    hydraulic_conductivity = hydraulic_conduction(
        new_matrix_potential,
        properties['water_retention'])

    # [W m-1 K-1]
    thermal_conductivity = thermal_conduction(new_volumetric_water)

    # return fluxes and state variables
    fluxes = {
        'net_radiation': net_radiation,  # [W m-2]
        'latent_heat': latent_heat_flux,  # [W m-2]
        'sensible_heat': sensible_heat_flux,  # [W m-2]
        'ground_heat': ground_heat_flux,  # [W m-2]
        'heat_advection': heat_advection,  # [W m-2]
        'water_closure': water_closure,  # [mm s-1]
        'energy_closure': energy_closure,  # [W m-2]
        'heat_fluxes': heat_fluxes,  # [W m-2]
        'evaporation': evaporation_rate,  # [mm s-1]
        'pond_recharge': pond_recharge / dt,  # [mm s-1]
        'capillar_rise': new_capillar_water / dt,  # [mm s-1]
        'throughfall': throughfall / dt,  # [mm s-1]
        }

    states = {
        'volumetric_water': new_volumetric_water,  # [m3 m-3]
        'water_potential': new_matrix_potential,  # [m]
        'water_content': new_gravimetric_water,  # [g g-1]
        'water_storage': new_water_content,  # [kg m-2] or [mm]
        'temperature': new_temperature,  # [degC]
        'hydraulic_conductivity': hydraulic_conductivity,  # [m s-1]
        'thermal_conductivity': thermal_conductivity,  # [W m-1 K-1]
        }

    return fluxes, states

def water_exchange(dt,
                   water_storage,
                   properties,
                   forcing):
    """
    Args:
        dt (float): timestep [s]
        water_content (float): [kg m-2]
        properties (dict): forestfloor object properties
        states (dict):
            'air_temperature': [degC]
            'wind_speed': [m s-1]
            'precipitation': [mm s-1]
            'h2o': [mol mol-1]
            'pond_storage': [mm]
            'soil_water_potential': [m]
            'soil_hydraulic_conductivity'
            'soil_depth'

    Returns:
        water_content (float): [kg m-2]
    """

    if dt == 0.0:
        dt = dt + EPS

    # [kg m-2] or [mm]
    water = water_storage
    max_water = properties['max_water_content'] * properties['dry_mass']
    min_water = properties['min_water_content'] * properties['dry_mass']

    water = min(max_water, water)
    water = max(min_water, water)

    max_recharge = max(max_water - water, 0.0)
    max_recharge = min(max_water, max_recharge)

    # [kg m-2 s-1] or [mm s-1]
    max_recharge_rate = max_recharge / dt

    # [kg m-2 s-1] or [mm s-1]
    max_evaporation_rate = (water_storage - (min_water + EPS)) / dt

    if np.isinf(max_evaporation_rate) or max_evaporation_rate < 0.0:
        max_evaporation_rate = 0.0

    # [kg m-2 s-1] or [mm s-1]
    max_condensation_rate = -((max_water - EPS) - water_storage) / dt

    if np.isinf(max_condensation_rate) or max_condensation_rate > 0.0:
        max_condensation_rate = 0.0

    # boundary layer conductances for H2O, heat and CO2
    # [mol m-2 s-1]
    conductance_to_air = moss_atm_conductance(
        forcing['wind_speed'],
        properties['roughness_height'],
        dT=0.0
    )

    # water vapor conductance from moss to air
    # Relative conductance is from Williams and Flanagan (1996), Oecologia.
    # Crosses 1 at appr. when water content/dry mass is 8.8
    relative_conductance = min(
        1.0,
        (0.1285 * water_storage / properties['dry_mass'] - 0.1285)
    )

    # [mol m-2 s-1]
    conductance_to_air_h2o = (
        conductance_to_air['h2o'] * relative_conductance
    )

    # [kg m-2 s-1]
    evaporation_rate = (
        conductance_to_air_h2o
        * (611.0 / forcing['air_pressure']
           * np.exp(17.502 * forcing['air_temperature'] / (forcing['air_temperature'] + 240.0))
           - forcing['h2o'])
        * MOLAR_MASS_H2O
    )

    evaporation_rate = min(evaporation_rate, max_evaporation_rate)
    evaporation_rate = max(evaporation_rate, max_condensation_rate)

    max_recharge_rate = max(max_recharge_rate + evaporation_rate, 0.0)
    max_recharge = max_recharge_rate * dt

    # take into account that maximum water capacity is not fully achieved

    # [mm] or [kg m-2]
    interception = (
        max_recharge
        * (1.0 - np.exp(-(1.0 / max_water)
                        * forcing['precipitation'] * dt))
    )

    # [kg m-2 s-1] or [mm s-1]
    interception_rate = interception / dt

    # [kg m-2 s-1] or [mm s-1]
    max_recharge = max(
        max_recharge - interception,
        0.0
    )

    pond_recharge = (
        max_recharge
        * (1.0 - np.exp(-(1.0 / max_water)
                        * forcing['soil_pond_storage'] * dt))
    )

    # [kg m-2 s-1] or [mm s-1]
    pond_recharge_rate = pond_recharge / dt

    # [kg m-2 s-1] or [mm s-1]
    max_recharge_rate = max(
        max_recharge - pond_recharge,
        0.0
    ) / dt

    # [g g-1]
    water_content = (water / properties['dry_mass'])

    # [m m-3]
    volumetric_water = (water_content / WATER_DENSITY
                        * properties['bulk_density'])

    # [m]
    water_potential = convert_hydraulic_parameters(
        volumetric_water,
        properties['water_retention'],
        'volumetric_water'
    )

    # [m s-1]
    hydraulic_conductivity = hydraulic_conduction(
        water_potential,
        properties['water_retention']
    )

    # [kg m-2 s-1] or [mm s-1]
    capillary_rise = capillarity(
        dt=dt,
        properties=properties,
        hydraulic_conductivity=hydraulic_conductivity,
        water_potential=water_potential,
        water_content=water_content,
        soil_hydraulic_conductivity=forcing['soil_hydraulic_conductivity'],
        soil_water_potential=forcing['soil_water_potential'],
        soil_depth=forcing['soil_depth']
    )

    # [kg m-2 s-1] or [mm s-1]
    capillary_rise = min(capillary_rise, max_recharge_rate)

    # [kg m-2 s-1] or [mm s-1]
    max_recharge_rate = max(max_recharge_rate - capillary_rise, 0.0)

    # calculate mass balance of water

    # [kg m-2 s-1] or [mm s-1]
    dy_water = (
        interception_rate
        + pond_recharge_rate
        + capillary_rise
        - evaporation_rate
    )

    # [kg m-2] or [mm]
    new_water_storage = dy_water + water_storage

    # [g g-1]
    new_water_content = new_water_storage / properties['dry_mass']

    # [m3 m-3]
    new_volumetric_water = (
        new_water_content
        / WATER_DENSITY
        * properties['bulk_density']
    )

    new_water_potential = convert_hydraulic_parameters(
        new_volumetric_water,
        properties['water_retention'],
        'volumetric_water'
    )

    hydraulic_conductivity = hydraulic_conduction(
        new_water_potential,
        properties['water_retention']
    )

    thermal_conductivity = thermal_conduction(new_volumetric_water)

    fluxes = {
        'evaporation': evaporation_rate,  # [mm s-1]
        'capillar_rise': capillary_rise,  # [mm s-1]
        'pond_recharge': pond_recharge_rate,  # [mm s-1]
        'throughfall': forcing['precipitation'] - interception_rate  # [mm s-1]
    }

    states = {
        'volumetric_water': new_volumetric_water,  # [m3 m-3]
        'water_potential': new_water_potential,  # [m]
        'water_content': new_water_content,  # [g g-1]
        'water_storage': new_water_storage,  # [kg m-2] or [mm]
        'hydraulic_conductivity': hydraulic_conductivity,  # [m s-1]
        'thermal_conductivity': thermal_conductivity,  # [W m-1 K-1]
    }

    return fluxes, states


def capillarity(dt,
                properties,
                hydraulic_conductivity,
                water_potential,
                water_content,
                soil_hydraulic_conductivity,
                soil_water_potential,
                soil_depth):
    r""" Estimates liquid water flux from soil to moss
    :math:`q = -k(\\frac{\partial h}{\partial z}+1)` [mm s\ :sup:`-1`] to moss
    as capillary rise from underlying soil. Capillary rise is further limited
    by bryophytes available water storage capacity.

    For water flow in soil, multiply capillary_rise by bryophyte
    ground_coverage and add that to the root_sink.

    Args:
        dt: [s]
        properties: dictionary containing characteristics of BryoType object
            * 'height'
            * 'max_water_content'
            * 'dry_mass'
        hydraulic_conductivity: [m s\ :sup:`-1`\ ]
        water_potential: [m]
        water_content: [g g\ :sup:`-1`\ ]
        soil_hydraulic_conductivity: [m s\ :sup:`-1`\ ]
            from 1st calculation node
        soil_water_potential: [m]
            from 1st calculation node
        soil_depth: [m]

    Returns:
        float:
            capillary rise in [mm s\ :sup:`-1`\ ]
    """

    # [m/s] geometric average of hydrological conductivities of
    # soil and bryophyte
    conductivities = hydraulic_conductivity * soil_hydraulic_conductivity
    k = np.power(conductivities, 0.5)

    # [mm/s] or [kg/(m2 s)]
    capillary_rise = (1.0e3
                      * max(0.0,
                            - k * ((water_potential - soil_water_potential)
                                   / (properties['height'] - soil_depth) + 1.0)))

    capillary_rise = min(
        capillary_rise,
        (properties['max_water_content'] - water_content)
        * properties['dry_mass'] / dt)

    return capillary_rise


def saturation_vapor_pressure(temperature):
    r""" Calculates saturation vapor pressure over water surface.

    Args:
        temperature: [\ :math:`^{\circ}`\ C]

    Returns:
        float: saturation vapor pressure in [Pa]
    """

    # [Pa]
    return 611.0 * np.exp((17.502 * temperature) / (temperature + 240.97))


def vapor_conductance_porous_media(volumetric_water,
                                   porosity,
                                   depth,
                                   temperature,
                                   ambient_pressure=101300.0):
    r""" Estimates molecular conductance of CO\ :sub:`2` and H\ :sub:`2`\ O in
    porous media.

    Following asumptions are made:
        1. all gas exchange is due to molecular diffusion
        2. relative diffusvity in soil is to that of air (D/D\ :sub:`o`\ ) as
           described by Millington and Quirk (1961)

    Args:
        volumetric_water: [m\ :sup:`3` m\ :sup:`-3`\ ]
        porosity: [-]
        depth: [m]
            Transport distance
        ambient_pressure: [Pa]

    Returns:
        dictionary:
            molecular conductances for
                * 'CO2': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'H2O': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
    """

    # [K]
    temperature = temperature + DEG_TO_KELVIN

    # [m3/m3], air filled porosity
    afp = np.maximum(0.0, porosity - volumetric_water)
    # afp = max(0, porosity - volumetric_water)

    # [mol m-3], air molar density
    cair = ambient_pressure / (GAS_CONSTANT * temperature)

    # D/Do, diffusivity in porous media relative to that in free air,
    # Millington and Quirk (1961)
    relative_diffusivity = (np.power((temperature / DEG_TO_KELVIN), 1.75)
                            * np.power(afp, 10.0/3.0) / porosity**2)

    # [mol/(m2 s)], through moss depth: mol m-1 s-1 x m-1
    conductance_h2o = cair * MOLECULAR_DIFFUSIVITY_H2O * relative_diffusivity / depth
    conductance_co2 = cair * MOLECULAR_DIFFUSIVITY_CO2 * relative_diffusivity / depth

    return {
        'co2': conductance_co2,
        'h2o': conductance_h2o
        }


def thermal_conduction(volumetric_water, method='donnel'):
    r""" Estimates heat conductivity (km) of bryophyte layer.

    By default organic matter heat conductivity is calculated by using equation
    by O'Donnel et al. (2009, Soil Sci. 174).

    Args:
        volumetric_water: [m\ :sup:`3` m\ :sup:`-3`\ ]
        flag (optional):
            optional methods are:
                * 'campbell' (Campbell et al. (1985))
                * 'constant' (0.25)
                * 'lauren' (Lauren (1999, table 4))

    Returns:
        float: heat conductivity in [W m\ :sup:`-1` K\ :sup:`-1`\ ]
    """

    method = method.lower()

    heat_conductivity = None

    if method == 'donnel':  # O'Donnell
        heat_conductivity = np.minimum(
            0.6,
            3.0e-2 + 5.0e-1 * volumetric_water)

        heat_conductivity = np.maximum(4.0e-2, heat_conductivity)

    elif method == 'campbell':
        heat_conductivity = (
            0.4 + 0.5
            * volumetric_water(0.4 - 0.06)
            * np.exp(-(1. * volumetric_water)**4))

    elif method == 'constant':
        heat_conductivity = 0.25

    elif method == 'lauren':
        heat_conductivity = -0.004 + 0.609 * volumetric_water

    # [W/(m K)]
    return heat_conductivity


def soil_boundary_layer_conductance(u, z, zo, Ta, dT, P=101300.):
    """
    Computes soil surface boundary layer conductance (mol m-2 s-1)
    assuming velocity profile logarithmic between z and z0.
    INPUT: u - mean velocity (m/s)
           z - height of mean velocity u (m)
           zo - soil surface roughness length for momentum (m)
           Ta - ambient temperature (degC)
           dT - soil surface-air temperature difference (degC)
           P - pressure(Pa)
    OUTPUT: boundary-layer conductances (mol m-2 s-1)
        gb_h - heat (mol m-2 s-1)
        gb_c- CO2 (mol m-2 s-1)
        gb_v - H2O (mol m-2 s-1)
    Based on Daamen & Simmons (1996). Note: gb decreases both in
    unstable and stable conditions compared to near-neutral;
    nonfeasible?
    Samuli Launiainen, 18.3.2014
    to python by Kersti
    """

    u = np.maximum(u, EPS)

    rho_air = 44.6*(P / 101300.0)*(273.15 / (Ta + 273.13))  # molar density of air [mol/m3]

    delta = 5.0 * GRAVITY * z * dT / ((Ta + 273.15) * u**2)
    if delta > 0:
        d = -0.75
    else:
        d = -2
    rb = (np.log(z/zo))**2 / (0.4**2*u)*(1 + delta)**d

    gb_h = rho_air * 1 / rb
    gb_v = MOLECULAR_DIFFUSIVITY_H2O / THERMAL_DIFFUSIVITY_AIR * gb_h
    gb_c = MOLECULAR_DIFFUSIVITY_CO2 / THERMAL_DIFFUSIVITY_AIR * gb_h

    return gb_h, gb_c, gb_v


def moss_atm_conductance_old(wind_speed, roughness_height):
    r""" Estimates boundary layer conductance of bryophyte canopy.

    Simple version without free convection

    Wind speed should represent vertical wind speed at ca. 20 cm above moss
    canopy (parametrsization derived from wind-tunnel studies). Roughness
    lengths scale is equal to the 'characteristic vertical height of
    individual bryophyte shoots' (typically order of 3-10 mm).

    Estimated CO\ :sub:`2` and heat conductances are related to water vapor
    by the ratio of molecular diffusivities
    (Campbell and Norman, 1998, eq. 7.29-7.33).

    References:
        Rice et al., 2001.
            Significance of variation in bryophyte canopy structure.
            Amer. J. Bot. 88:1568-1576.
        Rice, 2006.
            Towards an integrated undestanding of Bryophyte performance:
            the dimensions of space and time.
            Lindbergia 31:41-53.

    Args:
        wind_speed: [m s\ :sup:`-1`\ ]
        roughness_height: [m]

    Returns:
        dictionary:
            boundary layer conductances for
                * 'co2': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'h2o': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'heat': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
    NOTE:
        Add here approximation for free convection;
        in Campbell & Norman sect 7:

        g_free['heat'] = 0.05 * [(Tsur - Tair) / d]**0.25,
        where d [m] is characteristic dimension.
        If assume moss stems as cylinders, d ~stem height.
        g_free['h2o'] = 1.09 x g_free['heat'],
        g_free['co2'] = 0.75 x g_free['co2']

        The pain in the ass is the characteristic dimension which is hard
        to define... for range of d and Tsur - Tair we get values that
        may indicate uncertainty range... in any case these are comparable
        to conductances due to forced convection (see fig)
    """

    schmidt_number_h2o = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_H2O
    reynolds_number = wind_speed * roughness_height / AIR_VISCOSITY

    # Rice et al. (2001) eq. 1
    conductance_h2o = (
        AIR_DENSITY
        * np.power(10, -3.18)
        * np.power(reynolds_number, 1.61)
        * MOLECULAR_DIFFUSIVITY_H2O
        / roughness_height
        * np.power(schmidt_number_h2o, 1.0/3.0))

    # [mol m-2 s-1], differ from H2O by ratio of molecular diffusivities
    conductance_co2 = (MOLECULAR_DIFFUSIVITY_CO2
                       / MOLECULAR_DIFFUSIVITY_H2O
                       * conductance_h2o)

    conductance_heat = (THERMAL_DIFFUSIVITY_AIR
                        / MOLECULAR_DIFFUSIVITY_H2O
                        * conductance_h2o)

    return {
        'co2': conductance_co2,
        'h2o': conductance_h2o,
        'heat': conductance_heat
        }


def moss_atm_conductance(wind_speed, roughness_height, dT=0.0, atten_factor=0.25):
    r""" Estimates boundary layer conductance of bryophyte canopy for paralell
    forced and free convection.

    Wind speed should represent vertical wind speed at ca. 20 cm above moss
    canopy (parametrsization derived from wind-tunnel studies). Roughness
    lengths scale is equal to the 'characteristic vertical height of
    individual bryophyte shoots' (typically order of 3-10 mm).

    Estimated CO\ :sub:`2` and heat conductances are related to water vapor
    by the ratio of molecular diffusivities
    (Campbell and Norman, 1998, eq. 7.29-7.33).

    References:
        Rice et al., 2001.
            Significance of variation in bryophyte canopy structure.
            Amer. J. Bot. 88:1568-1576.
        Rice, 2006.
            Towards an integrated undestanding of Bryophyte performance:
            the dimensions of space and time.
            Lindbergia 31:41-53.
        Schuepp, 1980.
            Observations on the use of analytical and numerical models for the
            description of transfer to porous surface vegetation such as
            lichen.
            Boundary-Layer Meteorol. 29: 59-73.
        Kondo & Ishida, 1997

    Args:
        wind_speed: [m s\ :sup:`-1`\ ]
        roughness_height: [m]
        deltaT: [degC], moss-air temperature difference
        atten_factor: [-] dimensionless attenuation factor for continuous moss carpets

    Returns:
        dictionary:
            boundary layer conductances for
                * 'co2': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'h2o': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'heat': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]

    """

    Sc_v = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_H2O
    Sc_c = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_CO2
    Pr = AIR_VISCOSITY / THERMAL_DIFFUSIVITY_AIR
    Re = wind_speed * roughness_height / AIR_VISCOSITY


    # Rice et al. (2001) eq. 1: ShSc**0.33 = CRe**n, where C=6.6e-4 and n=1.61.
    # however, exponent n for individual species is <1.53 so use median values of model
    # fitted to individual species.
    C = 0.0067
    n = 1.27

    Sh_v = atten_factor * C*Re**n * Sc_v**0.33 # Sherwood numbner for H2O

    conductance_h2o = Sh_v * MOLECULAR_DIFFUSIVITY_H2O / roughness_height # ms-1

    # free convection as parallell pathway, based on Condo and Ishida, 1997.
    b = 2.2e-3 #ms-1K-1 b=1.1e-3 for smooth, 3.3e-3 for rough surface
    dT = np.maximum(dT, 0.0)
    gfree = Sc_v / Pr * b * dT**0.33  # mol m-2 s-1

    # [mol m-2 s-1]
    conductance_h2o = (conductance_h2o + gfree) * AIR_DENSITY


    # [mol m-2 s-1], differ from H2O by ratio of molecular diffusivities
    conductance_co2 = Sc_c / Sc_v * conductance_h2o

    conductance_heat = Pr / Sc_v * conductance_h2o

    return {
        'co2': conductance_co2,
        'h2o': conductance_h2o,
        'heat': conductance_heat
        }


def evaporation_through(properties,
                        volumetric_water,
                        moss_temperature,
                        air_temperature,
                        atm_partial_pressure_h2o,
                        wind_speed,
                        ambient_pressure,
                        soil_temperature,
                        soil_hydraulic_head,
                        soil_hydraulic_conductivity,
                        soil_depth):
    r""" Estimates soil evaporation rate through bryophyte layer.

    Evaporation in bryophyte layer is limited either by atmospheric demand
    and transport or soil supply of water.

    Water vapor flow from soil to air must overcome two resistances
    in series: 1. molecular diffusion through porous living moss, 2. molecular
    diffusion from moss canopy to 1st caluclation node in the atomosphere. The
    2nd resistance is assumed to be equal to conductance from wet moss canopy.

    Method does not use state variables of BryoType instance due to iteration.

    Args:
        properties (dict): characteristics of BryoType object
        volumetric_water: [m\ :sup:`3` m\ :sup:`-3`\ ]
        air_temperature: [\ :math:`^{\circ}`\ C]
        atm_partial_pressure_h2o: [Pa]
        wind_speed: [m s\ :sup:`-1`\ ]
        ambient_pressure: [Pa]
        soil_temperature: [\ :math:`^{\circ}`\ C]
            from 1st calculation node
        soil_hydraulic_head: [m]
            from 1st calculation node
        soil_hydraulic_conductivity: [m s\ :sup:`-1`\ ]
            from 1st calculation node
        soil_depth: [m]
            as negative value

    Returns:
        float:
            evaporation in [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
    """

    # [mol/(m2 s)]
    # "from soil surface to moss canopy height"
    # change air_temperature to bryo_temperature
    moss_conductance = vapor_conductance_porous_media(
        volumetric_water,
        properties['porosity'],
        properties['height'],
        moss_temperature,
        ambient_pressure)

    # [mol/(m2 s)]
    # "from moss canopy to the atmosphere"
    temp_difference = moss_temperature - air_temperature
    atm_conductance = moss_atm_conductance(wind_speed,
                                           properties['roughness_height'],
                                           dT=temp_difference)

    # [mol/(m2 s)], two resistors in series
    conductance_h2o = (
        moss_conductance['h2o'] * atm_conductance['h2o']
        / (moss_conductance['h2o'] + atm_conductance['h2o']))

    # Assuming soil is saturated (rh = 1),calculate the maximum evaporation rate

    # [mol/(m2 s)]
    # atmospheric evaporative demand
    evaporative_demand = (conductance_h2o
                          * (saturation_vapor_pressure(soil_temperature) - atm_partial_pressure_h2o)
                          / ambient_pressure)

    # [-, fraction]
    relative_humidity = min(1.0, atm_partial_pressure_h2o
                         / saturation_vapor_pressure(air_temperature))

    # [m], in equilibrium with atmospheric relative humidity
    atm_hydraulic_head = (
        GAS_CONSTANT
        * (DEG_TO_KELVIN + air_temperature)
        * np.log(relative_humidity)
        / (MOLAR_MASS_H2O * GRAVITY))

    # [mol/(m2 s)]: 1e3 is water density kgm-3: 1e-3 kgm-3 / (kg/mol) x m/s = mol m-2 s-1
    evaporative_supply = max(0.0,
        1e3 / MOLAR_MASS_H2O * soil_hydraulic_conductivity
        * ((atm_hydraulic_head - soil_hydraulic_head) / soil_depth - 1.0))

    # [mol/(m2 s)]
    # return min(evaporative_demand, evaporative_supply)
    return {'evaporative_demand': evaporative_demand,
            'evaporative_supply': evaporative_supply,
            'soil_evaporation': min(evaporative_demand, evaporative_supply)}


def effective_saturation(value, water_retention, input_unit):
    r""" Calculates effective saturation from volumetric water
    content (theta, [m3 m-3]) or from water potential (psi, [m])

    Args:
        value (float):
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is VOLUMETRIC_WATER
            * [m] if input is WATER_POTENTIAL
        water_retention (dict):
            if input_unit == 'water_potential'
                'alpha'
                'n'
            if input_unit == VOLUMETRIC_WATER
                'theta_r'
                'theta_s'
        input_unit (str): VOLUMETRIC_WATER or WATER_POTENTIAL
    Returns:
        float: [-]
    """

    options = set(['water_potential', 'volumetric_water'])

    if isinstance(input_unit, str):
        input_unit = input_unit.lower()
    else:
        raise ValueError("Input unit has to be string")

    if input_unit not in options:
        raise ValueError("Input unit options are: ".format(*options))


    if input_unit == 'water_potential':
        n = water_retention['n']
        m = 1.0 - np.divide(1.0, n)
        psi = 100.0 * np.minimum(value, 0.0)

        eff_sat = (1.0 + (water_retention['alpha'] * abs(psi)) ** n) ** -m

    elif input_unit == 'volumetric_water':

        theta_r = water_retention['theta_r']
        theta_s = water_retention['theta_s']

        theta = np.minimum(value, theta_s)
        theta = np.maximum(theta, theta_r)

        eff_sat = np.divide(((theta - theta_r) + EPS), (theta_s - theta_r) + EPS)

    return eff_sat


def convert_hydraulic_parameters(value, water_retention, input_unit):
    r""" Converts between hydraulic parameters volumetric water content and
        water potential.

    Note:
        In case of heterogenous poresystem, linear interpolation is used
        to convert effective saturation to water potential. Water retention
        curve with 500 points with logarithmic spacing in which water potential
        is ranging from -1e-6 to -1e2. In case of effective saturation is out
        of interpolation bounds, water potential is set to 0.0 m in lower bound
        or -1e2 m in upper bound.

    For conversion, VanGenuchten-model for water retention is used with
    following parameters:
            - saturated water content (theta_s, [cm\ :sup:`3` cm\ :sup:`-3`\ ])
            - residual water content (theta_r, [cm\ :sup:`3` cm\ :sup:`-3`\ ])
            - air entry suction (alpha, [cm\ :sup:`-1`\ ])
            - pore size distribution (n, [-])

    Volumetric water content is restricted to be between residual water
    content and saturated water content.

    Water potential is restricted to be between -10000 m and 0.0 m.

    Args:
        value:
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is VOLUMETRIC_WATER
            * [m] if input is WATER_POTENTIAL
        water_retention (dict): water retension parameters of BryoType object
        input_unit:
                * VOLUMETRIC_WATER
                * WATER_POTENTIAL
    Returns:
        float:
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input was WATER_POTENTIAL
            * [m] if input was VOLUMETRIC_WATER
    """

    if isinstance(input_unit, str):
        input_unit.lower()

    theta_s = water_retention['theta_s']  # sat. water content [m3m-3]
    theta_r = water_retention['theta_r']  # residual water content [m3m-3]

    if input_unit == 'water_potential':
        eff_sat = effective_saturation(value,
                                       water_retention,
                                       input_unit)

        return eff_sat * (theta_s - theta_r) + theta_r

    elif input_unit == 'volumetric_water':

        alpha = water_retention['alpha']  # micropore air entry suction [cm-1]
        n = water_retention['n']  # micropore shape parameter [-]

        eff_sat = effective_saturation(value, water_retention, input_unit)
        inv_effs = 1.0 / eff_sat
        m = np.divide(n, n - 1.0)

        # [m]
        psi = 1e-2 * 1.0 / -alpha * (inv_effs ** m - 1.0) ** (1.0 / n)

        if isinstance(psi, list):
            psi = np.array(list)

        if isinstance(psi, np.ndarray):
            psi[psi <= -1000.0] = -1000.0
            #psi = np.where(psi <= -100.0, -100.0, psi)
        else:
            if psi <= -1000.0:
                return -1000.0

        return psi


def hydraulic_conduction(water_potential, water_retention, method=None):
    r""" Estimates hydraulic conductivity.

    Hydraulic conductivity is based on Voortman et al. (2014) and
    Voortman et al. (2015) for xerophilious mosses.

    The implementation is valid for unimodal and multimodal poresystems.
    Multimodal Mualem-VanGenuchten is used as introduced by
    Priesack and Durner (2006).

    09.08.2018 AJK:
    At the moment hydraulic conduction used only for micropore
    system (capillary pore system) and it is only needed to calculate capillary
    rise. Macropore system does not retain water and it is assumed that
    macropore system does not restrict hydraulic conduction.

    References:
        Voortman et al. (2014) Hydrological Processes, 28, 6251-6264
        Priesack and Durner (2006) Vadose Zone Journal, 5, 121-124

    Lauren:
    Estimation is based on Lauren (1999, table 4) and parametrization is valid
    from -4 to -80 kPa. Assumption is that in the dry states
    (lower than -80 kPa) there are no hydraulic conductivity.

    Args:
        water_potential: [m]
        water_retention (list):
            parameters of water retention curve
                0. saturated water content (theta_s) [m\ :sup:`3` m :sup:`-3`\ ]
                1. residual water content (theta_r) [m\ :sup:`3` m :sup:`-3`\ ]
                2. air entry suction (alpha) [cm\ :sup:`-1`]
                3. pore size distribution (n) [-]
                4. saturated conductivity (K_s) [m s\ :sup:`-1`]
                5. pore connectivity (l) [-]

    Returns:
        float: hydraulic conductivity in [m s\ :sup:`-1`\ ]
    """

    if isinstance(method, str):
        method.lower()

    saturated_conductivity = water_retention['saturated_conductivity']

    if method == 'lauren':
        # kPa (standard gravity 10 m/s2)
        saturated_conductivity = water_retention['saturated_conductivity']
        water_potential = -10.0 * water_potential

        # Assuming that in dry states there is no hydraulic conductivity
        if water_potential > 80.0:
            return 0.0

        if water_potential < 4.0:
            water_potential = 4.0

        # relative conductivity respect to 4 kPa
        conductivity = (
            np.power(10.0, 1.0/(-0.62 + 0.26 * np.log10(water_potential)))
            / np.power(10.0, 1.0/(-0.62 + 0.26 * np.log10(4)))
            )
        # [m/s]
        return conductivity * saturated_conductivity

    else:
        # Possibility to add more poresystems

        psi = 100.0 * np.minimum(water_potential, 0.0)

        micropores = {'alpha': water_retention['alpha'],
                      'n': water_retention['n']}

        poresystems = [micropores]

        coefficients = []
        denominators = []
        nominators = []

        for idx in range(len(poresystems)):
            m = 1.0 - np.divide(1.0, poresystems[idx]['n'])
            alpha = poresystems[idx]['alpha']

            sat_eff = effective_saturation(psi,
                                           poresystems[idx],
                                           'water_potential')

            coefficients.append(sat_eff)

            denom = 1.0 - np.power(sat_eff, 1.0 / m)
            denom = 1.0 - np.power(denom, m)
            denominators.append(alpha * denom)

            nominators.append(alpha)

        pore_connectivity = water_retention['pore_connectivity']

        coefficient = np.power(np.sum(coefficients, axis=0), pore_connectivity)

        denominator = np.power(np.sum(denominators, axis=0), 2.0)

        nominator = np.power(np.sum(nominators, axis=0), 2.0)

        return saturated_conductivity * coefficient * (denominator / nominator)


def emitted_longwave_radiation(temperature, properties=None):
    r""" Estimates emitted longwave radiation

    Args:
        temperature (float): [W m\ :sup:`-2`]
        properties (dict/float): properties dictionary or emissivity
    Returns:
        (float): [W m\ :sup:`-2`]
    """
    if isinstance(properties, dict):
        emissivity = properties['optical_properties']['emissivity']

    elif isinstance(properties, float):
        emissivity = properties

    else:
        emissivity = 0.98

    emitted_longwave_radiation = (
        emissivity * STEFAN_BOLTZMANN
        * np.power((temperature + DEG_TO_KELVIN), 4.0))

    return emitted_longwave_radiation


def bryophyte_shortwave_albedo(water_content, properties=None):
    r""" Estimates hydration status scaled albedo for PAR and NIR regions.

    The effect of water content of spectral properties are based on studies by
    Vogelmann and Moss (1993) and Fernandes (1999)on Sphagnum cuspidatum and
    Pleurozium schreberi, respectively.

    The albedo is scaled specific reflectance for PAR (400-700 nm) or
    NIR (750-1400) regions. The scaling coefficient is common for both
    PAR and NIR and it is based on relationship between normalized
    reflectaces and hydration status. The species specific albedo is
    assumed to represent a reflectance when bryophyte is in full hydration.

    If bryophyte's properties are not given, estimation is based on generic
    fit of water content against reflectances separately for PAR and NIR.
    Fits are based on studies by Vogelmann and Moss (1993), and
    Fernandes (1999) on Sphagnum cuspidatum and Pleurozium schreberi,
    respectively.

    References:
        Vogelmann and Moss (1993)
            Remote Sensing of Environment 45:273-279.
        Fernandes (1999)
            PhD thesis entitled: 'Scale influences of surface
            parametrization on modelled boreal carbon and water budgets'

    Args:
        water_content (float): [g g\ :sup:`-2`]
        max_water_content (float): [g g\ :sup:`-2`]

    Returns:
            list: reflectances
                0. par albedo
                1. nir albedo
    """

    if properties is None:

        def reflectance(x, a, b):
            r""" Describes relationship between water content and reflectance

            Args:
                x (float): water_content
                a (float): fitting parameter
                b (float): fitting parameter

            Returns:
                percentage (float)
            """
            return np.abs(a * np.power(x, b))

        # Original reflectance in percents
        if isinstance(water_content, float):

            if water_content < 0.3:
                nir = 68.41/100.0

            else:
                nir = reflectance(water_content, 47.66, -0.3002) / 100.0

            par = reflectance(water_content, 8.84, -0.1303) / 100.0

        else:

            nir = np.empty(np.shape(water_content))
            nir = reflectance(water_content, 47.66, -0.3002) / 100.0
            nir[water_content < 0.3] = 68.41 / 100.0

            par = reflectance(water_content, 8.84, -0.1303) / 100.0

        return {'PAR': par, 'NIR': nir}

    else:
        albedo_nir = properties['optical_properties']['albedo_NIR']
        albedo_par = properties['optical_properties']['albedo_PAR']
        normalized_water_content = water_content / properties['max_water_content']

        scaling_coefficient = 0.2 / (1 - 0.9 * np.power(1.8703, -normalized_water_content))

        albedo_par = scaling_coefficient * albedo_par
        albedo_nir = scaling_coefficient * albedo_nir

        return {'PAR': albedo_par, 'NIR': albedo_nir}
