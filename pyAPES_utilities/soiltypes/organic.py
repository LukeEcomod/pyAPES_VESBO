# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:23:48 2018

@author: L1656
"""
from pyAPES_utilities.parameter_utilities import fit_pF, peat_hydrol_properties

plot=False

# depth of layer bottom [m], soil surface at 0.0
zh = [-0.1, -0.2, -0.3, -0.4, -0.5, -1.0, -2.0]
N = len(zh)
# pf based on bulk density
bd = [0.0811, 0.140, 0.174, 0.171, 0.119, 0.119, 0.119]
pf_para, _ = peat_hydrol_properties(bd, fig=plot, labels=['layer ' + str(i) for i in range(N)])
porosity = [pf_para[k][0] for k in range(N)]
residual_water_content = [pf_para[k][1] for k in range(N)]
pf_alpha = [pf_para[k][2] for k in range(N)]
pf_n = [pf_para[k][3] for k in range(N)]
# Hydraulic conductivity [m s-1]
Kvsat = [1.7e-4, 2e-5, 5e-5, 5e-6, 3e-6, 1e-6, 1e-7]  # vertical
Khmult = [10.0, 10.0, 5.0, 5.0, 5.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat
Khsat = [Kvsat[i] * Khmult[i] for i in range(N)]

soil_properties = {'pF': {  # vanGenuchten water retention parameters
                         'ThetaS': porosity,
                         'ThetaR': residual_water_content,
                         'alpha': pf_alpha,
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
