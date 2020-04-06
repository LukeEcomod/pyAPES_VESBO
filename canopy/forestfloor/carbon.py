#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:23:28 2018

Module contains functions to calculate carbon exchange in bryophytes.

Note:
    migrated to python3
    - no changes made

@author: ajkieloaho
"""

import numpy as np
from canopy.constants import EPS
from canopy.forestfloor.bryophoto import photo_farquhar, relative_capacity, conductance

class BryophyteFarquhar(object):
    r"""Bryophyte photosynthesis and respiration using community-level Farquhar-model
        Gives fluxes per unit ground area.
    """
    def __init__(self, para, carbon_pool=0):
        self.photopara = para
        self.carbon_pool = carbon_pool

    def co2_exchange(self, Qp, Ca, T, w, wstar=None):
        """
        computes net CO2 exchange of moss community
        Args:
            Qp - incident PAR [umol m-2 s-1]
            Ca - ambient CO2 [ppm]
            T - moss temperature (degC)
            w - moss water content (g/g)
            wstar - delayed water content (g/g) for desiccation recovery
        Returns:
            An - net CO2 exchange An = -A + Rd, <0 is uptake [umol m-2 s-1]
            A - photosynthesis rate [umol m-2 s-1]
            Rd - dark respiration rate [umol m-2 s-1]
            Cc - internal CO2 (ppm)
            g - total conductance for CO2 [umol m-2 s-1]
        Note:
            units of An, A, Rd and g depend on whether Vcmax etc. are given
            per ground area, per dry mass or something else...
        """
        p = self.photopara.copy()
        cap, rcap = relative_capacity(p, w, wstar)
        p['Vcmax'] *= cap
        p['Jmax'] *= cap
        p['alpha'] *= cap
        p['Rd'] *= rcap

        # conductance (mol m-2 s-1)
        g = conductance(p, w)

        # solve Anet and Cc iteratively until Cc converges
        err = 10^9
        Cc = 0.8*Ca

        while err > 1e-3:
            Cco = Cc
            An, Rd, _, _ = photo_farquhar(p, Qp, Cc, T)

            Cc = Ca - An / g  # new Cc
            Cc = 0.5*(Cco + Cc)
            err = np.nanmax(abs(Cc - Cco))

        return {'net_co2': -An,
                'photosynthesis': An + Rd,
                'respiration': Rd,
                'internal_co2': Cc,
                'conductance_co2': g
                }

class BryophyteCarbon(object):
    r"""Bryophyte community photosynthesis and respiration using simple
    light-response with empirical temperature and moisture functions.

    Gives fluxes per unit ground area.
    """
    def __init__(self, para, carbon_pool=0):
        self.amax = para['photosynthesis']['amax']
        self.b = para['photosynthesis']['b']
        self.moisture_coeff = para['photosynthesis']['moisture_coeff']
        self.temperature_coeff = para['photosynthesis']['temperature_coeff']

        self.r10 = para['respiration']['r10']
        self.q10 = para['respiration']['q10']

        self.max_water_content = para['max_water_content'] # g g-1
        self.carbon_pool = carbon_pool

    def carbon_exchange(self,
                        temperature,
                        water_content,
                        incident_par):
        r""" Estimates photosynthesis and respiration rates of bryophyte layer.

        Photosynthesis is restricted by both tissue water content
        (dry conditions) and excess water film on leaves (diffusion limitation)
        as in Williams and Flanagan (1996). Water content
        coefficients are 3rd order polynomial fitted to the data represented by
        Williams and Flanagan (1996) and used to calculate effect of water
        content on photosynthesis.

        Empirical modifier of photosynthesis due to water content assumes that
        both light-limited and Rubisco-limited assimilation of carbon are
        affected similarly. This seems to apply for Pleurozium and Sphagnum
        when normalized water content is used as scaling. Assumes that
        there is always 5 percents left in photosynthetic capacity.Empirical
        modifier of photosynthesis due to temperature is based on
        late growing season presented in Fig. 2 in Frolking et al. (1996).

        References:
            Frolking et al. (1996)
                Global Change Biology 2:343-366
            Williams and Flanagan (1996)
                Oecologia 108:38-46

        Args:
            water_content (float): [g g-1]
            temperature (float): [degC]
            incident_par (float): [W m\ :sup:`-2`]

        Returns:
            dictionary:
                * 'photosynthesis_rate':
                  [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
                * 'respiration_rate':
                  [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
                 * 'co_flux': net co2 exchange, <0 uptake
                  [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
        """
        # check inputs
        # [umol/(m2 s)]
        incident_par = np.maximum(EPS, 4.56 * incident_par)

        normalized_water_content = water_content / self.max_water_content


        # hyperbolic light response at community level [umolm-2(ground) s-1]
        light_response = self.amax * incident_par / (self.b + incident_par )

        # moisture and temperature responses [-]
        water_modifier = (self.moisture_coeff[3]
                          + self.moisture_coeff[2] * normalized_water_content
                          + self.moisture_coeff[1] * normalized_water_content ** 2.0
                          + self.moisture_coeff[0] * normalized_water_content ** 3.0)

        water_modifier = np.maximum(0.05, water_modifier)

        temperature_modifier = (self.temperature_coeff[3]
                                + self.temperature_coeff[2] * temperature
                                + self.temperature_coeff[1] * temperature ** 2.0
                                + self.temperature_coeff[0] * temperature ** 3.0)

        temperature_modifier = np.maximum(0.01, temperature_modifier)

        temperature_modifier = np.minimum(1.0, temperature_modifier)

        # [umol m-2 (leaf) s-1]
        photosynthetic_rate = light_response * water_modifier * temperature_modifier

        # --- respiration rate [umol m-2 (ground) s-1]

        """ replace with smooth function and move as parameter """
        if water_content < 7.0:
            water_modifier_respiration = (
                -0.45 + 0.4 * water_content
                - 0.0273 * water_content ** 2)
        else:
            water_modifier_respiration = (
                -0.04 * water_content + 1.38)

        # effect of water content is in the range 0.01 to 1.0
        water_modifier_respiration = np.maximum(0.01, np.minimum(1.0, water_modifier_respiration))

        # r = r10 * Q10^((T-10) / 10) [umol m-2 s-1]
        respiration_rate = (
            self.r10 * self.q10**((temperature - 10.0) / 10.0) * water_modifier_respiration
            )

        # umol m-2 (ground) s-1
        return {
            'photosynthesis': photosynthetic_rate,
            'respiration': respiration_rate,
            'net_co2': -photosynthetic_rate + respiration_rate
            }

class OrganicRespiration():
    """
    Litter layer respiration. Values per unit ground area
    """
    def __init__(self, para, carbon_pool=0.0):
        self.r10 = para['r10']
        self.q10 = para['q10']
        #self.moisture_coeff = para['respiration']['moisture_coeff']
        self.carbon_pool = 0.0

    def respiration(self, temperature, volumetric_water):
        """
        respiration rate [umol m-2 (ground) s-1]
        """
        r = self.r10 * np.power(self.q10, (temperature - 10.0) / 10.0)

        # add moisture response
        fW = 1.0
        r = r * fW
        return {'photosynthesis': 0.0,
                'respiration': r,
                'net_co2': r
               }


class SoilRespiration():
    """
    Soil respiration
    """
    def __init__(self, para, weights=1):
        # base rate [umol m-2 s-1] and temperature sensitivity [-]
        self.r10 = para['r10']
        self.q10 = para['q10']

        # moisture response of Skopp et al. 1990
        self.moisture_coeff = para['moisture_coeff']

        if weights is not None:
            # soil respiration computed in layers and weighted
            self.weights = weights
            self.Nlayers = len(weights)

    def respiration(self, soil_temperature, volumetric_water, volumetric_air):
        """ Soil respiration beneath forestfloor

        Heterotrophic and autotrophic respiration rate (CO2-flux) based on
        Pumpanen et al. (2003) Soil.Sci.Soc.Am

        Restricts respiration by soil moisuture as in
        Skopp et al. (1990), Soil.Sci.Soc.Am

        Args:

            Ts - soil temperature [degC]
            Wliq - soil vol. moisture content [m3 m-3]
        Returns:
            soil respiration rate [umol m-2 s-1]

        """
        # Skopp limitparam [a,b,d,g] for two soil types
        # sp = {'Yolo':[3.83, 4.43, 1.25, 0.854],
        #       'Valentine': [1.65,6.15,0.385,1.03]}

        # unrestricted respiration rate
        x = self.r10 * np.power(self.q10, (soil_temperature - 10.0) / 10.0)

        # moisture response (substrate diffusion, oxygen limitation)
        f = np.minimum(self.moisture_coeff[0] * volumetric_water**self.moisture_coeff[2],
                       self.moisture_coeff[1] * volumetric_air**self.moisture_coeff[3])
        f = np.minimum(f, 1.0)

        respiration = x * f

        if hasattr(self, 'weights'):
            respiration = sum(self.weights * respiration[0:self.Nlayers])

        return respiration
