# -*- coding: utf-8 -*-
"""
Soilprofile
"""

def get_spara(soiltype):
    """ grid """

    grid = {'dz': [0.1],
            'zh': [-0.1]
            }

    soil_properties = {'pF': {  # vanGenuchten water retention parameters
                             'ThetaS': [0.89],
                             'ThetaR': [1e-12],
                             'alpha': [0.96],
                             'n': [1.24]
                             },
                      'saturated_conductivity_vertical': [2.42e-05],  # saturated vertical hydraulic conductivity [m s-1]
                      'saturated_conductivity_horizontal': [2.42e-05],  # saturated horizontal hydraulic conductivity [m s-1]
                      'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
                      'solid_composition': {  # fractions of solid volume [-]
                                   'organic': [0.2076],
                                   'sand': [0.2621],
                                   'silt': [0.5152],
                                   'clay': [0.0151]
                                   },
                      'freezing_curve': [0.2],  # freezing curve parameter
                      'bedrock': {
                                  'solid_heat_capacity': 2.16e6,  # [J m-3 (solid) K-1]
                                  'thermal_conductivity': 3.0  # thermal conductivity of non-porous bedrock [W m-1 K-1]
                                  }
                      }

    """ water model specs """
    water_model = {'solve': False,
                   'type': 'Richards',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
                   'pond_storage_max': 0.05,  #  maximum pond depth [m]
                   'initial_condition': {  # (dict) initial conditions
#                                         'ground_water_level': -1.0,  # groundwater depth [m]
                                         'pond_storage': 0.  # initial pond depth at surface [m]
                                         },
                   'lower_boundary': {  # lower boundary condition (type, value, depth)
                                      'type': 'impermeable',
                                      'value': None,
                                      'depth': -2.0
                                      },
                   'drainage_equation': {  # drainage equation and drainage parameters
                                         'type': 'Hooghoudt',
                                         'depth': 1.0,  # drain depth [m]
                                         'spacing': 45.0,  # drain spacing [m]
                                         'width': 1.0,  # drain width [m]
                                         }
                    }
    
    """ heat model specs """
    heat_model = {'solve': False,
                  'initial_condition': {
                                       'temperature': 4.0,  # initial soil temperature [degC]
                                       },
                  'lower_boundary': {  # lower boundary condition (type, value)
                                    'type': 'temperature',
                                    'value': 4.0
                                    },
            }

    spara = {'grid': grid, 'soil_properties': soil_properties, 'water_model': water_model, 'heat_model': heat_model}

    return spara


#from .soiltypes.organic import soil_properties as organic_soil_properties
#from .soiltypes.mineral import soil_properties as mineral_soil_properties
#
#def get_spara(soiltype):
#    """ grid """
#    # thickness of layers with different characteristics [m]
#    thickness = [0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0]
#    N = len(thickness)
#    # number of nodes in each layer
#    nodes = [10, 5, 5, 2, 2, 5, 5]
#    # computational layer depths [m]
#    dz = sum([[thickness[i] / nodes[i] for k in range(nodes[i])] for i in range(N)], [])
#    
#    if soiltype.upper() == 'ORGANIC':
#        soil_properties, zh = organic_soil_properties()
#    else:
#        soil_properties, zh = mineral_soil_properties(soiltype.upper())
#    
#    grid = {'dz': dz,
#            'zh': zh
#            }
#    
#    """ water model specs """
#    water_model = {'solve': True,
#                   'type': 'Richards',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
#                   'pond_storage_max': 0.05,  #  maximum pond depth [m]
#                   'initial_condition': {  # (dict) initial conditions
#                                         'ground_water_level': -0.2,  # groundwater depth [m]
#                                         'pond_storage': 0.  # initial pond depth at surface [m]
#                                         },
#                   'lower_boundary': {  # lower boundary condition (type, value, depth)
#                                      'type': 'impermeable',
#                                      'value': None,
#                                      'depth': -1.8
#                                      },
#                   'drainage_equation': {  # drainage equation and drainage parameters
#                                         'type': 'Hooghoudt',
#                                         'depth': 1.0,  # drain depth [m]
#                                         'spacing': 45.0,  # drain spacing [m]
#                                         'width': 1.0,  # drain width [m]
#                                         }
#                    }
#    
#    """ heat model specs """
#    heat_model = {'solve': True,
#                  'initial_condition': {
#                                       'temperature': 4.0,  # initial soil temperature [degC]
#                                       },
#                  'lower_boundary': {  # lower boundary condition (type, value)
#                                    'type': 'temperature',
#                                    'value': 4.0
#                                    },
#            }
#
#    spara = {'grid': grid, 'soil_properties': soil_properties, 'water_model': water_model, 'heat_model': heat_model}
#
#    return spara
