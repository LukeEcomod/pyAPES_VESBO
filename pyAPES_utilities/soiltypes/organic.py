# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:23:48 2018

@author: L1656


Parameters for Lettosuo peat soil
"""
from pyAPES_utilities.parameter_utilities import fit_pF, peat_hydrol_properties

plot=False

# depth of layer bottom [m], soil surface at 0.0
zh = [-0.1, -0.2, -0.3, -0.4, -0.5, -1.0, -2.0]
N = len(zh)
# pf based on bulk density
bd = [0.1047, 0.1454, 0.1591,0.1300, 0.1119, 0.1591, 0.1591]
vp = [1, 5.75, 4.5, 4.5, 4.75, 8, 8]

pf_para, Ksat = peat_hydrol_properties(vp,  var='H', fig=plot, labels=['layer ' + str(i) for i in range(N)], ptype='C')
#pf_para, Ksat = peat_hydrol_properties(bd,  var='bd', fig=plot, labels=['layer ' + str(i) for i in range(N)], ptype='C')

## raw humus from Laihos measurements
## heads [kPa]
#head = [0.01, 0.3, 0.981, 4.905, 9.81, 33.0, 98.1]
## volumetric water content [%]
#watcont = [[94.69, 49.42, 29.61, 21.56, 20.05, 17.83, 16.54]] # 0-10 cm, näyte 1
##watcont = [[91.98, 66.70, 57.49, 39.95, 34.41, 29.83, 28.39]] # 0-10 cm, näyte 2
#pf_para[0] = fit_pF(head, watcont, fig=plot, percentage=True, kPa=True)[0]

## Laiho Lettosuo
#head = [0.01, 0.3, 0.981, 4.905, 9.81, 33.0, 98.1]
## volumetric water content
#watcont = [[94.69, 49.42, 29.61, 21.56, 20.05, 17.83, 16.54], # 0-10 cm
#           [91.41, 66.26, 56.98, 45.58, 41.44, 39.32, 37.89], # 10-20 cm
#           [89.12, -999, 72.83, 63.97, 54.40, 50.15, 48.80], # 20-30 cm
#           [89.46, -999, 82.46, 76.79, 66.93, 63.61, 62.53], # 30-40 cm
#           [92.22, -999, 87.06, 78.02, 74.76, 72.77, 71.70]]  # 40-50 cm
#pf_para=fit_pF(head, watcont, fig=False,percentage=True, kPa=True)

porosity = [pf_para[k][0] for k in range(N)]
residual_water_content = [pf_para[k][1] for k in range(N)]
pf_alpha = [pf_para[k][2] for k in range(N)]
pf_n = [pf_para[k][3] for k in range(N)]
# function of humification and depth (Päivänen 1973)
#Kvsat = [5.0e-5, 2.3e-5, 2.3e-5, 1.3e-5, 5.9e-6, 1.7e-6, 1.7e-6]
# TEST
Kvsat = Ksat
#Kvsat[-2:] = 1.0e-7

Khmult = [30.0, 20.0, 10.0, 5.0, 1.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat

# TEST
Kvsat = [4.97E-05, 3.21E-05, 2.07E-05, 1.34E-05, 8.63E-06, 2.32E-06, 2.61E-07]
Khmult = [30.0, 10.0, 5.0, 1.0, 1.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat

Khsat = [Kvsat[i] * Khmult[i] for i in range(N)]

#Khmult1 = [20.0, 10.0, 5.0, 5.0, 1.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat
#Khmult2 = [20.0, 10.0, 5.0, 5.0, 5.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat
#Khmult3 = [30.0, 20.0, 5.0, 1.0, 1.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat
#Khmult4 = [30.0, 20.0, 5.0, 5.0, 1.0, 1.0, 1.0]  # horizontal Khsat = Khmult * Kvsat

soil_properties = {'pF': {  # vanGenuchten water retention parameters
                         'ThetaS': porosity,
                         'ThetaR': residual_water_content,
                         'alpha': pf_alpha,
                         'n': pf_n
                         },
                  'saturated_conductivity_vertical': Kvsat,  # saturated vertical hydraulic conductivity [m s-1]
                  'saturated_conductivity_horizontal': Khsat,
#                                                        [Kvsat[i] * Khmult1[i] for i in range(N)]),
#                                                        [Kvsat[i] * Khmult2[i] for i in range(N)],
#                                                        [Kvsat[i] * Khmult3[i] for i in range(N)],
#                                                        [Kvsat[i] * Khmult4[i] for i in range(N)]),# saturated horizontal hydraulic conductivity [m s-1]
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
