# -*- coding: utf-8 -*-
"""
Soilprofile
"""
from soiltypes.organic import soil_properties as organic_soil_properties
from soiltypes.mineral import soil_properties as mineral_soil_properties

def get_spara(soiltype):
    """ grid """
    # thickness of layers with different characteristics [m]
    thickness = [0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0]
    N = len(thickness)
    # number of nodes in each layer
    nodes = [10, 5, 5, 2, 2, 5, 5]
    # computational layer depths [m]
    dz = sum([[thickness[i] / nodes[i] for k in range(nodes[i])] for i in range(N)], [])
    
    if soiltype.upper() == 'ORGANIC':
        soil_properties, zh = organic_soil_properties()
    else:
        soil_properties, zh = mineral_soil_properties(soiltype.upper())
    
    grid = {'dz': dz,
            'zh': zh
            }
    
    """ water model specs """
    water_model = {'solve': True,
                   'type': 'Richards',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
                   'pond_storage_max': 0.05,  #  maximum pond depth [m]
                   'initial_condition': {  # (dict) initial conditions
                                         'ground_water_level': -0.2,  # groundwater depth [m]
                                         'pond_storage': 0.  # initial pond depth at surface [m]
                                         },
                   'lower_boundary': {  # lower boundary condition (type, value, depth)
                                      'type': 'impermeable',
                                      'value': None,
                                      'depth': -1.8
                                      },
                   'drainage_equation': {  # drainage equation and drainage parameters
                                         'type': 'Hooghoudt',
                                         'depth': 1.0,  # drain depth [m]
                                         'spacing': 45.0,  # drain spacing [m]
                                         'width': 1.0,  # drain width [m]
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

    spara = {'grid': grid, 'soil_properties': soil_properties, 'water_model': water_model, 'heat_model': heat_model}

    return spara