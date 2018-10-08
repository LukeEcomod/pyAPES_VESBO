# -*- coding: utf-8 -*-
"""
Soilprofile
"""
import numpy as np
from utilities import fit_pF, peat_hydrol_properties

""" grid """
# thickness of layers with different characteristics [m]
thickness = [0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0]
N = len(thickness)
# depth of layer bottom [m], soil surface at 0.0
zh = -np.cumsum(thickness)
# number of nodes in each layer
nodes = [10, 5, 5, 2, 2, 5, 5]
# computational layer depths [m] --- USE AS INPUT??
dz = sum([[thickness[i] / nodes[i] for k in range(nodes[i])] for i in range(N)], [])

grid = {'dz': dz,
        'zh':zh
        }

""" soil properties """
# pf based on bulk density
bd = [0.0811, 0.140, 0.174, 0.171, 0.119, 0.119, 0.119]
pf_para, _ = peat_hydrol_properties(bd)#, fig=True
porosity = [pf_para[k][0] for k in range(N)]
residual_water_content = [pf_para[k][1] for k in range(N)]
pf_alpha = [pf_para[k][2] for k in range(N)]
pf_n = [pf_para[k][3] for k in range(N)]
# Hydraulic conductivity [m s-1]
Kvsat = [1.7e-4, 2e-5, 5e-5, 5e-6, 3e-6, 1e-6, 1e-7]  # vertical
Khmult = [10.0, 10.0, 5.0, 5.0, 5.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat
Khsat = [Kvsat[i] * Khmult[i] for i in range(N)]

soil_properties = {'pF': {  # vanGenuchten water retention parameters
                         'ThetaS': np.array(porosity),
                         'ThetaR': residual_water_content,
                         'alpha': pf_alpha,
                         'n': pf_n
                         },
                  'saturated_conductivity_vertical': Kvsat,  # saturated vertical hydraulic conductivity [m s-1]
                  'saturated_conductivity_horizontal': Khsat,  # saturated horizontal hydraulic conductivity [m s-1]
                  'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
                  'fractions': {  # fractions of solid volume [-]
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

""" water model specs """
water_model = {'solve': True,
               'type': 'Richards',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
               'pond_storage_max': 0.01,  #  maximum pond depth [m]
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
                                      'spacing': 40.0,  # drain spacing [m]
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