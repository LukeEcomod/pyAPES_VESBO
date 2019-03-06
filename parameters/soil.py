# -*- coding: utf-8 -*-
"""
Soil submodel parameters

"""

""" grid and soil properties """

grid = {#thickness of computational layer [m]
        'dz': [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
               0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
               0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2],
        # bottom depth of layers with different characteristics [m]
        'zh': [-0.05, -0.11, -0.35, -10.0]
        }

soil_properties = {'pF': {  # vanGenuchten water retention parameters
                        'ThetaS': [0.915, 0.50, 0.50, 0.50],
                        'ThetaR': [0.058, 0.08, 0.08, 0.08],
                        'alpha': [0.075, 0.06, 0.06, 0.06],
                        'n': [1.335, 1.35, 1.35, 1.35]
                        },
                  'saturated_conductivity_vertical': [2.42E-05, 1.23E-05, 1.23E-05, 1.23E-05],  # saturated vertical hydraulic conductivity [m s-1]
                  'saturated_conductivity_horizontal': [2.42E-05, 1.23E-05, 1.23E-05, 1.23E-05],  # saturated horizontal hydraulic conductivity [m s-1]
                  'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
                  'solid_composition': {
                         'organic': [0.1611, 0.0714, 0.1091, 0.028],
                         'sand': [0.4743, 0.525, 0.5037, 0.5495],
                         'silt': [0.3429, 0.3796, 0.3641, 0.3973],
                         'clay': [0.0217, 0.0241, 0.0231, 0.0252]
                         },
                  'freezing_curve': [0.2, 0.5, 0.5, 0.5],  # freezing curve parameter
                  'bedrock': {
                              'solid_heat_capacity': 2.16e6,  # [J m-3 (solid) K-1]
                              'thermal_conductivity': 3.0  # thermal conductivity of non-porous bedrock [W m-1 K-1]
                              }
                  }

""" water model specs """
water_model = {'solve': True,
               'type': 'Richards',  #'Equilibrium', #'Richards',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
               'pond_storage_max': 0.05,  #  maximum pond depth [m]
               'initial_condition': {  # (dict) initial conditions
                       'ground_water_level': -2.0,  # groundwater depth [m]
                       'pond_storage': 0.  # initial pond depth at surface [m]
                       },
               'lower_boundary': {  # lower boundary condition (type, value, depth)
                       'type': 'head_oneway',
                       'value': -2.0,
#                       'type': 'impermeable',
#                       'value': None,
#                       'depth': -2.0
                       },
               'drainage_equation': {  # drainage equation and drainage parameters
                       'type': None,  #
#                       'type': 'Hooghoudt',  #
#                       'depth': 1.0,  # drain depth [m]
#                       'spacing': 45.0,  # drain spacing [m]
#                       'width': 1.0,  # drain width [m]
                       }
                }

""" heat model specs """
heat_model = {'solve': True,
              'initial_condition': {
                      'temperature': 4.0,  # initial soil temperature [degC]
                      },
              'lower_boundary': {  # lower boundary condition (type, value)
                      'type': 'temperature',
                      'value': 4.0
                      },
              }

spara = {'grid': grid,
         'soil_properties': soil_properties,
         'water_model': water_model,
         'heat_model': heat_model}