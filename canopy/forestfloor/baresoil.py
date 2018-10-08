#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 08:44:28 2018

@author: ajkieloaho
"""

import numpy as np

from canopy.constants import EPS
from canopy.constants import STEFAN_BOLTZMANN, LATENT_HEAT, DEG_TO_KELVIN
from canopy.constants import MOLAR_MASS_H2O, WATER_DENSITY, GRAVITY
from canopy.constants import SPECIFIC_HEAT_AIR, SPECIFIC_HEAT_H2O, GAS_CONSTANT
from canopy.micromet import e_sat

from heat_and_water import soil_boundary_layer_conductance


class BaresoilModel(object):
    """ Represents bare soil cover-soil-atmosphere interactions

    """

    def __init__(self, properties, initial_temperature=None):
        self.properties = properties

        self.coverage = properties['ground_coverage']
        self.poros = properties['porosity']

        self.albedo = {
                'PAR': properties['optical_properties']['albedo_PAR'],
                'NIR': properties['optical_properties']['albedo_NIR']
                }

        self.emissivity = properties['optical_properties']['emissivity']
        self.soil_zr = properties['roughness_length']

        if initial_temperature is not None:
            self.temperature = initial_temperature
        else:
            self.temperature = 10.

        self.old_temperature = self.temperature

    def update(self):
        """ Update states to new states after iteration
        """
        self.old_temperature = self.temperature

    def run(self, dt, forcing):
        """ Calculates one timestep and updates states of BaresoilModel instance
        """

        states, fluxes = heat_balance(forcing, self.properties)

        return fluxes, states

def heat_balance(forcing, properties):
    """ Solves bare soil surface temperature

    Uses linearized energy balance equation from soil conditions from
    previous timestep

    Args:
        forcing (dict)
        properties (dict)
    """

    z_can = forcing['height']
    U = forcing['wind_speed']
    T = forcing['air_temperature']
    H2O = forcing['h2o']
    P = forcing['air_pressure']
    T_ave = forcing['forestfloor_temperature']

    zr = properties['roughness_length']
    T_soil = forcing['soil_temperature']
    h_soil = forcing['soil_water_potential']
    z_soil = forcing['depth']
    Kh = forcing['soil_hydraulic_conductivity']
    Kt = forcing['soil_thermal_conductivity']

    Par_gr = forcing['par']
    Nir_gr = forcing['nir']

    LWn = forcing['lw_net']
    Ebal = forcing['Ebal']  # energy balance switch

    albedo_par = properties['optical_properties']['albedo_PAR']
    albedo_nir = properties['optical_properties']['albedo_NIR']
    soil_emi = properties['optical_properties']['emissivity']

    dz_soil = - z_soil
    # initial guess for surface temperature
    T_surf = T
    # boundary layer conductances for H2O and heat [mol m-2 s-1]
    gb_h, _, gb_v = soil_boundary_layer_conductance(u=U, z=z_can, zo=zr, Ta=T, dT=0.0, P=P)  # OK to assume dt = 0.0?
    # radiative conductance [mol m-2 s-1]
    gr = 4.0 * soil_emi * STEFAN_BOLTZMANN * T_ave**3 / SPECIFIC_HEAT_AIR

    # absorbed shortwave radiation
    SW_gr = (1 - albedo_par) * Par_gr + (1 - albedo_nir) * Nir_gr

    # Maximum LE
    # atm pressure head in equilibrium with atm. relative humidity
    es_a, _ = e_sat(T)
    RH = H2O * P / es_a  # air relative humidity above ground [-]
    h_atm = GAS_CONSTANT * (DEG_TO_KELVIN + T) * np.log(RH)/(MOLAR_MASS_H2O * GRAVITY)  # [m]
    # maximum latent heat flux constrained by h_atm
    LEmax = -LATENT_HEAT * Kh * (h_atm - h_soil - z_soil) / dz_soil * WATER_DENSITY / MOLAR_MASS_H2O  # [W/m2]

    # LE demand
    # vapor pressure deficit between leaf and air, and slope of vapor pressure curve at T
    es, s = e_sat(T_surf)
    Dsurf = es / P - H2O  # [mol/mol] - allows condensation
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
        Told = T_surf
        if Ebal:
            # solve leaf temperature [degC]
            T_surf = (SW_gr + LWn + SPECIFIC_HEAT_AIR*gr*T_ave + SPECIFIC_HEAT_AIR*gb_h*T - LE + LATENT_HEAT*s*gb_v*Told
                      + Kt / dz_soil * T_soil) / (SPECIFIC_HEAT_AIR*(gr + gb_h) + LATENT_HEAT*s*gb_v + Kt / dz_soil)
            err = np.nanmax(abs(T_surf - Told))
#            print ('iterNo', iterNo, 'err', err, 'T_surf', T_surf)
            es, s = e_sat(T_surf)
            Dsurf = es / P - H2O  # [mol/mol] - allows condensation
            s = s / P  # [mol/mol/degC]
            LE = LATENT_HEAT * gb_v * Dsurf
            if LE > LEmax:
                LE = LEmax
                s = 0.0
            if iterNo == itermax:
                print 'Maximum number of iterations reached in surface energy module'
                print('err', err, 'T_surf', T_surf)
        else:
            err = 0.0

    """ --- energy and water fluxes --- """
    # sensible heat flux [W m-2]
    Hw = SPECIFIC_HEAT_AIR * gb_h * (T_surf - T)
    # non-isothermal radiative flux [W m-2]
    Frw = SPECIFIC_HEAT_AIR * gr *(T_surf - T_ave)
    # ground heat flux [W m-2]
    Gw = Kt / dz_soil * (T_surf - T_soil)
    # evaporation rate [mol m-2 s-1]
    Ep = gb_v * Dsurf

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
            'temperature': T_surf
            }

    return states, fluxes


# EOF