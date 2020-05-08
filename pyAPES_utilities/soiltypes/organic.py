# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:23:48 2018

@author: Kersti Haahti


Parameters for Lettosuo peat soil
"""

plot=False

# depth of layer bottom [m], soil surface at 0.0
zh = [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1., -1.5, -2.0]
N = len(zh)

# TEST
porosity = [0.943, 0.882, 0.882, 0.882, 0.882, 0.882, 0.882, 0.882, 0.882, 0.882, 0.882, 0.882]
residual_water_content = [0.002, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104]
pf_alpha = [0.202, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044]
pf_n = [1.349, 1.349, 1.349, 1.349, 1.349, 1.349, 1.349, 1.349, 1.349, 1.349, 1.349, 1.349]

Kvsat = [4.97E-05, 3.21E-05, 2.07E-05, 1.34E-05, 8.63E-06, 5.57E-06, 3.60E-06, 2.32E-06, 1.50E-06, 9.68E-07, 2.61E-07, 1.16E-07]
Khmult = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat

Khsat = [Kvsat[i] * Khmult[i] for i in range(N)]

soil_properties = {'pF': {  # vanGenuchten water retention parameters
                         'ThetaS': porosity,
                         'ThetaR': residual_water_content,
                         'alpha': pf_alpha,  # [cm-1]
                         'n': pf_n
                         },
                  'saturated_conductivity_vertical': Kvsat,  # saturated vertical hydraulic conductivity [m s-1]
                  'saturated_conductivity_horizontal': Khsat,  # saturated horizontal hydraulic conductivity [m s-1]
                  'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
                  'solid_composition': {  # fractions of solid volume [-]
                               'organic': [1.0 for i in range(N)],
                               'sand': [0.0 for i in range(N)],
                               'silt': [0.0 for i in range(N)],
                               'clay': [0.0 for i in range(N)]
                               },
                  'freezing_curve': [0.5 for i in range(N)],  # freezing curve parameter
                  'bedrock': {
                              'solid_heat_capacity': 2.16e6,  # [J m-3 (solid) K-1]
                              'thermal_conductivity': 3.0  # thermal conductivity of non-porous bedrock [W m-1 K-1]
                              }
                  }
