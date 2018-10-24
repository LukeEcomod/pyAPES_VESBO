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


def soil_respiration(properties, Ts, Wliq, Wair):
    """ Soil respiration beneath forestfloor

    Heterotrophic and autotrophic respiration rate (CO2-flux) based on
    Pumpanen et al. (2003) Soil.Sci.Soc.Am

    Restricts respiration by soil moisuture as in
    Skopp et al. (1990), Soil.Sci.Soc.Am

    Args:
        properties (dict):
            'R10'
            'Q10'
            'poros'
            'limitpara'
        Ts - soil temperature [degC]
        Wliq - soil vol. moisture content [m3 m-3]
    Returns:
        rsoil - soil respiration rate [umol m-2 s-1]
        fm - relative modifier (Skopp et al.)
    """
    # Skopp limitparam [a,b,d,g] for two soil types
    # sp = {'Yolo':[3.83, 4.43, 1.25, 0.854],
    #       'Valentine': [1.65,6.15,0.385,1.03]}

    limitpara = properties['limitpara']
    r10 = properties['R10']
    q10 = properties['Q10']

    # unrestricted respiration rate
    base_respiration = r10 * np.power(q10, (Ts - 10.0) / 10.0)

    # moisture response (substrate diffusion, oxygen limitation)
    modifier = np.minimum(limitpara[0] * Wliq**limitpara[2],
                          limitpara[1] * Wair**limitpara[3])  # ]0...1]
    modifier = np.minimum(modifier, 1.0)

    respiration = base_respiration * modifier

    return respiration, modifier


def carbon_exchange(properties,
                    water_content,
                    temperature,
                    incident_par,
                    water_content_coeff=None,
                    temperature_coeff=None):
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

    .. plot::
        import matplotlib.pyplot as plt
        import numpy as np
        water_content_coeff = [6.4355, -14.0605, 9.1867, -0.8720]
        x1 = np.linspace(0, 1, 100)
        y1 = [max(0.05, (water_content_coeff[3]
            + water_content_coeff[2] * i
            + water_content_coeff[1] * i ** 2
            + water_content_coeff[0] * i ** 3)) for i in x1]
        temperature_coeff = [-4.3e-5, -8.3e-4, 0.08, 0.1]
        x2 = np.linspace(0, 30, 100)
        y2  = [max(0.01, (temperature_coeff[3]
            + temperature_coeff[2] * i
            + temperature_coeff[1] * i ** 2
            + temperature_coeff[0] * i ** 3)) for i in x2]
        fig, ax1 = plt.subplots()
        ax2 = ax1.twiny()
        fig.subplots_adjust(bottom=0.2)
        ax2.set_frame_on(True)
        ax2.patch.set_visible(False)
        ax2.xaxis.set_ticks_position('bottom')
        ax2.xaxis.set_label_position('bottom')
        ax2.spines['bottom'].set_position(('outward', 40))
        ax1.plot(x1, y1, 'b-', label='relative water content')
        ax2.plot(x2, y2, 'r-', label='temperature')
        ax1.set_xlabel('Relative water content')
        ax2.set_xlabel('Temperature')
        ax1.set_ylabel('Relative effect on photosynthesis')
        hand1, label1 = ax1.get_legend_handles_labels()
        hand2, label2 = ax2.get_legend_handles_labels()
        plt.legend(hand1+hand2,label1+label2, frameon=False)
        plt.show()

    References:
        Frolking et al. (1996)
            Global Change Biology 2:343-366
        Williams and Flanagan (1996)
            Oecologia 108:38-46

    Args:
        incident_par (float): [W m\ :sup:`-2`]
            total photosynthetically active radiation (PAR) at ground level
        water_content_coeff (optional): list of fitting parameters
        temperature_coeff (optional): list of fitting parameters

    Returns:
        dictionary:
            * 'photosynthesis_rate':
              [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
            * 'respiration_rate':
              [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
    """
    # check inputs
    # [umol/(m2 s)]
    incident_par = np.maximum(EPS, 4.56 * incident_par)

    # [-, fraction]
    normalized_water_content = (
        water_content / properties['max_water_content'])

    # photosynthetic light response parameters
    amax = properties['photosynthesis']['Amax']
    b = properties['photosynthesis']['b']

    # respiration parameters
    q10 = properties['respiration']['Q10']
    r10 = properties['respiration']['Rd10']

    if water_content_coeff is None:
        water_content_coeff = [6.4355, -14.0605, 9.1867, -0.8720]

    if temperature_coeff is None:
        temperature_coeff = [-4.3e-5, -8.3e-4, 0.08, 0.1]

    # SL 12.9.2018: LET'S USE COMMUNITY-LEVEL LIGHT-RESPONSE
    # hyperbolic light response [umolm-2(leaf) s-1]
    light_response = amax * incident_par / (b + incident_par )

    # moisture and temperature responses [-]
    water_modifier = (water_content_coeff[3]
                      + water_content_coeff[2] * normalized_water_content
                      + water_content_coeff[1] * normalized_water_content ** 2.0
                      + water_content_coeff[0] * normalized_water_content ** 3.0)

    water_modifier = np.maximum(0.05, water_modifier)

    temperature_modifier = (temperature_coeff[3]
                            + temperature_coeff[2] * temperature
                            + temperature_coeff[1] * temperature ** 2.0
                            + temperature_coeff[0] * temperature ** 3.0)

    temperature_modifier = np.maximum(0.01, temperature_modifier)

    temperature_modifier = np.minimum(1.0, temperature_modifier)

    # [umol m-2 (leaf) s-1]
    community_photosynthetic_rate = light_response * water_modifier * temperature_modifier

    # upscale to field scale [umol/(m2 ground s)]
    photosynthetic_rate = community_photosynthetic_rate * properties['ground_coverage']

    # respiration rate [umol m-2 (ground) s-1]
    if water_content < 7.0:
        water_modifier_respiration = (
            -0.45 + 0.4 * water_content
            - 0.0273 * water_content ** 2)
    else:
        water_modifier_respiration = (
            -0.04 * water_content + 1.38)

    # effect of water content is in the range 0.01 to 1.0
    water_modifier_respiration = np.maximum(0.01, np.minimum(1.0, water_modifier_respiration))

    # r = r10 * Q10^((T-10) / 10) [umol/(m2 s)]
    respiration_rate = (
        r10 * q10**((temperature - 10.0) / 10.0) * water_modifier_respiration
        )

    # [umol/(m2 ground s)]
    respiration_rate = (
        respiration_rate * properties['ground_coverage'])

    return {
        "photosynthesis_rate": photosynthetic_rate,
        "respiration_rate": respiration_rate
        }


def carbon_exchange_old(properties,
                        water_content,
                        temperature,
                        incident_par,
                        attenuation_coefficient=0.8,
                        water_content_coeff=None,
                        temperature_coeff=None):
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

    .. plot::

        import matplotlib.pyplot as plt
        import numpy as np

        water_content_coeff = [6.4355, -14.0605, 9.1867, -0.8720]
        x1 = np.linspace(0, 1, 100)
        y1 = [max(0.05, (water_content_coeff[3]
            + water_content_coeff[2] * i
            + water_content_coeff[1] * i ** 2
            + water_content_coeff[0] * i ** 3)) for i in x1]

        temperature_coeff = [-4.3e-5, -8.3e-4, 0.08, 0.1]
        x2 = np.linspace(0, 30, 100)
        y2  = [max(0.01, (temperature_coeff[3]
            + temperature_coeff[2] * i
            + temperature_coeff[1] * i ** 2
            + temperature_coeff[0] * i ** 3)) for i in x2]

        fig, ax1 = plt.subplots()
        ax2 = ax1.twiny()

        fig.subplots_adjust(bottom=0.2)

        ax2.set_frame_on(True)
        ax2.patch.set_visible(False)
        ax2.xaxis.set_ticks_position('bottom')
        ax2.xaxis.set_label_position('bottom')
        ax2.spines['bottom'].set_position(('outward', 40))

        ax1.plot(x1, y1, 'b-', label='relative water content')
        ax2.plot(x2, y2, 'r-', label='temperature')

        ax1.set_xlabel('Relative water content')
        ax2.set_xlabel('Temperature')
        ax1.set_ylabel('Relative effect on photosynthesis')

        hand1, label1 = ax1.get_legend_handles_labels()
        hand2, label2 = ax2.get_legend_handles_labels()

        plt.legend(hand1+hand2,label1+label2, frameon=False)
        plt.show()

    References:
        Frolking et al. (1996)
            Global Change Biology 2:343-366
        Williams and Flanagan (1996)
            Oecologia 108:38-46

    Args:
        incident_par (float): [W m\ :sup:`-2`]
            total photosynthetically active radiation (PAR) at ground level
        attenuation_coefficient: [-]
            assumed diffuse radiation
        water_content_coeff (optional): list of fitting parameters
        temperature_coeff (optional): list of fitting parameters

    Returns:
        dictionary:
            * 'photosynthesis_rate':
              [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
            * 'respiration_rate':
              [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
    """
    # check inputs
    # [umol/(m2 s)]
    incident_par = max(EPS, 4.56 * incident_par)

    # [-, fraction]
    normalized_water_content = (
        water_content / properties['max_water_content'])

    # photosynthetic light response parameters
    amax = properties['photosynthesis'][0]
    b = amax / (2.0 * properties['photosynthesis'][1])

    # respiration parameters
    r10 = properties['respiration']['Rd10']
    q10 = properties['respiration']['Q10']

    if water_content_coeff is None:
        water_content_coeff = [6.4355, -14.0605, 9.1867, -0.8720]

    if temperature_coeff is None:
        temperature_coeff = [-4.3e-5, -8.3e-4, 0.08, 0.1]

    # divide canopy into 10 equall layers from top and compute par at each layer
    # effective_lai --> attenuation of light within the moss
    effective_lai = (properties['leaf_area_index'] / properties['ground_coverage'])
    cumulative_lai = np.linspace(effective_lai / 10.0, effective_lai, 10)

    # par at each layer
    par = incident_par * np.exp(-attenuation_coefficient * cumulative_lai)

    # hyperbolic light response [umolm-2(leaf) s-1]
    light_response = amax * par / (b * (par + EPS))

    # moisture and temperature responses [-]
    water_modifier = (water_content_coeff[3]
                      + water_content_coeff[2] * normalized_water_content
                      + water_content_coeff[1] * normalized_water_content ** 2.0
                      + water_content_coeff[0] * normalized_water_content ** 3.0)

    water_modifier = max(0.05, water_modifier)

    temperature_modifier = (temperature_coeff[3]
                            + temperature_coeff[2] * temperature
                            + temperature_coeff[1] * temperature ** 2.0
                            + temperature_coeff[0] * temperature ** 3.0)

    temperature_modifier = max(0.01, temperature_modifier)

    temperature_modifier = min(1.0, temperature_modifier)

    # [umol m-2 (leaf) s-1]
    leaf_photosynthetic_rate = light_response * water_modifier * temperature_modifier

    # upscale to field scale [umol/(m2 ground s)]
    photosynthetic_rate = ((effective_lai / 10.0) * sum(leaf_photosynthetic_rate)
                           * properties['ground_coverage'])

    # respiration rate [umol m-2 (ground) s-1]
    if water_content < 7.0:
        water_modifier_respiration = (
            -0.45 + 0.4 * water_content
            - 0.0273 * water_content ** 2)
    else:
        water_modifier_respiration = (
            -0.04 * water_content + 1.38)

    # effect of water content is in the range 0.01 to 1.0
    water_modifier_respiration = max(0.01, min(1.0, water_modifier_respiration))

    # r = r10 * Q10^((T-10) / 10) [umol/(m2 s)]
    respiration_rate = (
        r10 * q10**((temperature - 10.0) / 10.0) * water_modifier_respiration
        )

    # [umol/(m2 ground s)]
    respiration_rate = (
        respiration_rate * properties['ground_coverage'])

    return {
        "photosynthesis_rate": photosynthetic_rate,
        "respiration_rate": respiration_rate
        }