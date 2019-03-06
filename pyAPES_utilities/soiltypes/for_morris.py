# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 09:13:46 2019

@author: L1656
"""

""" organic soil profile """
N = 12
soil_organic = {
        'zh': [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1., -1.5, -2.0],
        'properties': {
              'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.915, 0.905, 0.893, 0.893, 0.854, 0.854, 0.854, 0.854, 0.854, 0.854, 0.854, 0.854],
                    'ThetaR': [0.058, 0.074, 0.083, 0.083, 0.063, 0.063, 0.063, 0.063, 0.063, 0.063, 0.063, 0.063],
                    'alpha': [0.075, 0.051, 0.038, 0.038, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022],
                    'n': [1.335, 1.344, 1.341, 1.341, 1.283, 1.283, 1.283, 1.283, 1.283, 1.283, 1.283, 1.283]
                    },
#              'saturated_conductivity_vertical': [2.0E-04, 5.0E-05, 2.0E-05, 5.75E-06, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07],  # saturated vertical hydraulic conductivity [m s-1]
#              'saturated_conductivity_horizontal': [10 * 2.0E-04, 10 * 5.0E-05, 5 * 2.0E-05, 5 * 5.75E-06, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07],  # saturated horizontal hydraulic conductivity [m s-1]
              'saturated_conductivity_vertical': [2.07E-05, 1.09E-05, 5.75E-06, 5.75E-06, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07],
              'saturated_conductivity_horizontal': [10 * 2.07E-05, 10 * 1.09E-05, 5 * 5.75E-06, 5 * 5.75E-06, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07, 8.43E-07],
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
        }