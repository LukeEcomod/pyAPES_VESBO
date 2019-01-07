#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 08:44:28 2018


NOTE:
    migrated to python3
    - only absolute imports

@author: ajkieloaho
"""

import logging
import numpy as np

from canopy.constants import EPS
from canopy.constants import STEFAN_BOLTZMANN, LATENT_HEAT, DEG_TO_KELVIN
from canopy.constants import MOLAR_MASS_H2O, WATER_DENSITY, GRAVITY
from canopy.constants import SPECIFIC_HEAT_AIR, SPECIFIC_HEAT_H2O, GAS_CONSTANT
from canopy.micromet import e_sat
from .heat_and_water import soil_boundary_layer_conductance

logger = logging.getLogger(__name__)


class Baresoil():
    """ Represents bare soil cover-soil-atmosphere interactions

    """

    def __init__(self, properties, initial_conditions=None):
        self.properties = properties
        self.coverage = properties['ground_coverage']

        if initial_conditions is not None:
            self.temperature = initial_conditions['temperature']
        else:
            self.temperature = 10.

        self.old_temperature = self.temperature

    def update(self):
        """ Update states to new states after iteration
        """
        self.old_temperature = self.temperature

    def restore(self):
        """ Restores new states back to states before iteration.
        """

        self.temperature = self.old_temperature

    def run(self, dt, forcing, parameters, controls):
        """ Calculates one timestep and updates states of Baresoil instance
        """

        states, fluxes = heat_balance(
            forcing=forcing,
            parameters=parameters,
            controls=controls,
            properties=self.properties,
            temperature=self.temperature)

        # update state variables
        self.temperature = states['temperature']

        return fluxes, states


def heat_balance(forcing, parameters, controls, properties, temperature):
    r""" Solves bare soil surface temperature

    Uses linearized energy balance equation from soil conditions from
    previous timestep

    Args:
        forcing (dict):
            'wind_speed': [m s-1]
            'air_temperature': [degC]
            'h2o': [mol mol\ :sup:`-1`\ ]
            'air_pressure': [Pa]
            'forestfloor_temperature': [degC]
            'soil_tempereture': [degC]
            'soil_water_potential': [m]
            'par': [W m-2]
            'nir': [W m-2]
            'lw_dn': [W m-2]
            'lw_up': [W m-2]
        parameters (dict):
            'soil_hydraulic_conductivity': [m s-1]
            'soil_thermal_conductivity': []
            'height': [m] height to the first canopy calculation node
            'soil_depth': [m] depth to the first soil calculation node
        controls (dict):
            'energy_balance': boolean
        properties (dict):
            'optical_properties':
                'emissivity'
                'albedo_PAR'
                'albedo_NIR'
    """

    U = forcing['wind_speed']
    T = forcing['air_temperature']
    P = forcing['air_pressure']
    T_ave = forcing['forestfloor_temperature']
    T_soil = forcing['soil_temperature']
    h_soil = forcing['soil_water_potential']

    z_soil = parameters['soil_depth']
    Kt = parameters['soil_thermal_conductivity']
    Kh = parameters['soil_hydraulic_conductivity']

    soil_emi = properties['optical_properties']['emissivity']
    # radiative conductance [mol m-2 s-1]
    gr = 4.0 * soil_emi * STEFAN_BOLTZMANN * T_ave**3 / SPECIFIC_HEAT_AIR

    if controls['energy_balance']:  # energy balance switch

        albedo_par = properties['optical_properties']['albedo_PAR']
        albedo_nir = properties['optical_properties']['albedo_NIR']

        # absorbed shortwave radiation
        SW_gr = (1 - albedo_par) * forcing['par'] + (1 - albedo_nir) * forcing['nir']

        # net longwave radiation
        LWn = forcing['lw_dn'] - forcing['lw_up']

        # initial guess for surface temperature
        surface_temperature = temperature
    else:

        SW_gr, LWn = 0.0, 0.0
        # set temperature to average of soil and air
        surface_temperature = forcing['air_temperature']
        # geometric mean of air_temperature and soil_temperature
#        surface_temperature = (
#            np.power(forcing['air_temperature'] * forcing['soil_temperature'], 0.5)
#        )

    dz_soil = - z_soil

# change this either to baresoil temperature or baresoil old_temperature

    # boundary layer conductances for forcing['h2o'] and heat [mol m-2 s-1]
    gb_h, _, gb_v = soil_boundary_layer_conductance(
        u=forcing['wind_speed'],
        z=parameters['height'],
        zo=properties['roughness_length'],
        Ta=forcing['air_temperature'],
        dT=0.0,
        P=forcing['air_pressure']
    )  # OK to assume dt = 0.0?

    # Maximum LE
    # atm pressure head in equilibrium with atm. relative humidity
    es_a, _ = e_sat(T)
    RH = min(1.0, forcing['h2o'] * P / es_a)  # air relative humidity above ground [-]
    h_atm = GAS_CONSTANT * (DEG_TO_KELVIN + T) * np.log(RH)/(MOLAR_MASS_H2O * GRAVITY)  # [m]
    # maximum latent heat flux constrained by h_atm
    LEmax = max(0.0,
                -LATENT_HEAT * Kh * (h_atm - h_soil - z_soil) / dz_soil * WATER_DENSITY / MOLAR_MASS_H2O)  # [W/m2]

    # LE demand
    # vapor pressure deficit between leaf and air, and slope of vapor pressure curve at T
    es, s = e_sat(surface_temperature)
    Dsurf = es / P - forcing['h2o']  # [mol/mol] - allows condensation
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
        Told = surface_temperature
        if controls['energy_balance']:
            # solve leaf temperature [degC]
            surface_temperature = (
                (SW_gr + LWn + SPECIFIC_HEAT_AIR*gr*T_ave
                 + SPECIFIC_HEAT_AIR*gb_h*T - LE + LATENT_HEAT*s*gb_v*Told
                 + Kt / dz_soil * T_soil)
                / (SPECIFIC_HEAT_AIR*(gr + gb_h) + LATENT_HEAT*s*gb_v + Kt / dz_soil)
            )

            err = abs(surface_temperature - Told)
            es, s = e_sat(surface_temperature)
            Dsurf = es / P - forcing['h2o']  # [mol/mol] - allows condensation
            s = s / P  # [mol/mol/degC]
            LE = LATENT_HEAT * gb_v * Dsurf

            if LE > LEmax:
                LE = LEmax
                s = 0.0
            if iterNo == itermax:
                logger.debug('Maximum number of iterations reached: T_baresoil = %.2f, err = %.2f',
                             surface_temperature, err)
        else:
            err = 0.0

    if (abs(surface_temperature - temperature) > 20 or np.isnan(surface_temperature)):  # into iteration loop? chech photo or interception
        logger.debug('Unrealistic baresoil temperature %.2f set to previous value %.2f: %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f',
                     surface_temperature, temperature,
                     U, T, forcing['h2o'], P, T_ave, T_soil, h_soil, SW_gr, LWn)
        surface_temperature = temperature
        es, s = e_sat(surface_temperature)
        Dsurf = es / P - forcing['h2o']  # [mol/mol] - allows condensation
        LE = LATENT_HEAT * gb_v * Dsurf
        if LE > LEmax:
            LE = LEmax

    """ --- energy and water fluxes --- """
    # sensible heat flux [W m-2]
    Hw = SPECIFIC_HEAT_AIR * gb_h * (surface_temperature - T)
    # non-isothermal radiative flux [W m-2]
    Frw = SPECIFIC_HEAT_AIR * gr *(surface_temperature - T_ave)
    # ground heat flux [W m-2]
    Gw = Kt / dz_soil * (surface_temperature - T_soil)
    # evaporation rate [mol m-2 s-1]
    Ep = LE / LATENT_HEAT  #gb_v * Dsurf

    # energy closure
    closure = SW_gr + LWn - Hw - LE - Gw

    fluxes = {
        'latent_heat': LE,
        'energy_closure': closure,
        'evaporation': Ep,
        'radiative_flux': Frw,
        'sensible_heat': Hw,
        'ground_heat': Gw
    }

    states = {
            'temperature': surface_temperature
            }

    return states, fluxes

# EOF
