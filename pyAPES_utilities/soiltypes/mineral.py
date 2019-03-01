# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:44:19 2018

@author: L1656

Soil parameters for Ct, Vt, Mt, and OMt forests in different soil horizons (O, A, B1, B2, and C)

theta, and h in horizons A-C form Jauhiainen (2004), and in horizon O from Jauhiainen (2004) for Vt and OMt,
and Lauren and Heiskanen (2001) for Ct and Mt. Porosity is checked from Westman (1991) for O-horizon

Khsat are for pine forest (Launianen, 2015)

vOrg, vSand, vSilt, vClay from Westman (1990)
"""

from pyAPES_utilities.parameter_utilities import fit_pF
from matplotlib import pyplot as plt

def soil_properties(sitetype, plot=False):
    horizons = ['O', 'A', 'B1', 'B2', 'C']
    if sitetype.upper() == 'OMT':
        zh = [-0.046, -0.18, -0.279, -0.378, -10.0]
        Khsat = [2.42e-05, 2.08e-06, 3.06e-06, 3.06e-06, 4.17e-06]
        watcont = [[0.78, 0.634, 0.279, -999, 0.237, 0.205, 0.190, -999, 0.104, 0.067],  # Jauhiainen (2004), porosity from Westman (1990) 
                   [-999, 0.497, 0.439, -999, 0.342, 0.273, 0.217, -999, 0.112, 0.065],  # Jauhiainen (2004)
                   [-999, 0.495, 0.456, -999, 0.366, 0.286, 0.232, -999, 0.128, 0.082],  # Jauhiainen (2004)
                   [-999, 0.458, 0.413, -999, 0.322, 0.255, 0.206, -999, 0.094, 0.065],  # Jauhiainen (2004)
                   [-999, 0.406, 0.357, -999, 0.296, 0.254, 0.216, -999, 0.101, 0.065]]  # Jauhiainen (2004)
        vOrg = [0.2699, 0.0990, 0.0946, 0.0946, 0.0698]
        vSand = [0.3262, 0.4026, 0.4045, 0.4045, 0.4156]
        vSilt = [0.3684, 0.4547, 0.4569, 0.4569, 0.4694]
        vClay = [0.0355, 0.0438, 0.0440, 0.0440, 0.0452]
    elif sitetype.upper() == 'MT':
        zh = [-0.042, -0.095, -0.2025, -0.31, -10.0]
        Khsat = [2.42e-05, 2.08e-06, 3.06e-06, 3.06e-06, 4.17e-06]
        watcont = [[0.89, -999, 0.5111, 0.3997,-999,-999, 0.2954, 0.2331, 0.1696, -999],  # Lauren and Heiskanen (2001), porority from Westman (1990)
                   [-999, 0.514, 0.442, -999, 0.360, 0.298, 0.234, -999, 0.125, 0.073],  # Jauhiainen (2004)
                   [-999, 0.526, 0.477, -999, 0.397, 0.323, 0.246, -999, 0.141, 0.089],  # Jauhiainen (2004)
                   [-999, 0.462, 0.413, -999, 0.345, 0.286, 0.222, -999, 0.112, 0.070],  # Jauhiainen (2004)
                   [-999, 0.413, 0.362, -999, 0.312, 0.274, 0.217, -999, 0.107, 0.055]]  # Jauhiainen (2004)
        vOrg = [0.2076, 0.0982, 0.1288, 0.1288, 0.0654]
        vSand = [0.2621, 0.2983, 0.2882, 0.2882, 0.3091]
        vSilt = [0.5152, 0.5863, 0.5664, 0.5664, 0.6076]
        vClay = [0.0151, 0.0172, 0.0166, 0.0166, 0.0178]
    elif sitetype.upper() == 'VT':
        zh = [-0.052, -0.095, -0.272, -0.439, -10.0]
        Khsat = [2.42e-05, 2.08e-06, 3.06e-06, 3.06e-06, 4.17e-06]
        watcont = [[0.91, 0.578, 0.344, -999, 0.250, 0.201, 0.185, -999, 0.143, 0.123],  # Jauhiainen (2004), porosity from Westman (1990)
                   [-999, 0.431, 0.391, -999, 0.293, 0.224, 0.167, -999, 0.085, 0.051],  # Jauhiainen (2004)
                   [-999, 0.514, 0.475, -999, 0.364, 0.283, 0.215, -999, 0.130, 0.071],  # Jauhiainen (2004)
                   [-999, 0.431, 0.391, -999, 0.293, 0.224, 0.167, -999, 0.085, 0.051],  # Jauhiainen (2004)
                   [-999, 0.404, 0.334, -999, 0.248, 0.170, 0.117, -999, 0.043, 0.027]]  # Jauhiainen (2004)
        vOrg = [0.1611, 0.0714, 0.1091, 0.1091, 0.0280]
        vSand = [0.4743, 0.5250, 0.5037, 0.5037, 0.5495]
        vSilt = [0.3429, 0.3796, 0.3641, 0.3641, 0.3973]
        vClay = [0.0217, 0.0241, 0.0231, 0.0231, 0.0252]
    elif sitetype.upper() == 'CT':
        zh = [-0.041, -0.07, -0.213, -0.356, -10.0]
        Khsat = [2.42e-05, 2.08e-06, 3.06e-06, 3.06e-06, 4.17e-06]
        watcont = [[0.94, -999, 0.4186, 0.3241, -999, -999, 0.2476, 0.2000, 0.1590, -999],  # Lauren and Heiskanen (2001), porosity from Westman (1990)
                   [-999, 0.545, 0.474, -999, 0.332, 0.266, 0.210,-999,  0.117, 0.066],  # Jauhiainen (2004)
                   [-999, 0.498, 0.436, -999, 0.313, 0.248, 0.203, -999, 0.123, 0.072],  # Jauhiainen (2004)
                   [-999, 0.412, 0.328, -999, 0.215, 0.166, 0.131, -999, 0.075, 0.040],  # Jauhiainen (2004)
                   [-999, 0.357, 0.246, -999, 0.123, 0.068, 0.048, -999, 0.026, 0.017]]  # Jauhiainen (2004)
        vOrg = [0.1063, 0.0930, 0.1063, 0.1063, 0.0356]
        vSand = [0.3147, 0.3194	, 0.3147, 0.3147, 0.3396]
        vSilt = [0.5745, 0.5831, 0.5745, 0.5745, 0.6200]
        vClay = [0.0045, 0.0046, 0.0045, 0.0045, 0.0049]
    else:
        raise ValueError("Unknown sitetype %s" % sitetype)

    N = len(zh)

    # wrc characteristics
    head = [0.001, 1., 10., 30.,32., 63., 100., 300., 1000., 16000.]
    pf_para = fit_pF(head, watcont, fig=plot, labels=horizons)
    porosity = [pf_para[k][0] for k in range(N)]
    residual_water_content = [pf_para[k][1] for k in range(N)]
    pf_alpha = [pf_para[k][2] for k in range(N)]
    pf_n = [pf_para[k][3] for k in range(N)]

    soil_properties = {'pF': {  # vanGenuchten water retention parameters
                             'ThetaS': porosity,
                             'ThetaR': residual_water_content,
                             'alpha': pf_alpha,
                             'n': pf_n
                             },
                      'saturated_conductivity_vertical': Khsat,  # saturated vertical hydraulic conductivity [m s-1]
                      'saturated_conductivity_horizontal': Khsat,  # saturated horizontal hydraulic conductivity [m s-1]
                      'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
                      'solid_composition': {  # volumetric fractions of solid volume [-]
                                   'organic': vOrg,
                                   'sand': vSand,
                                   'silt': vSilt,
                                   'clay': vClay
                                   },
                      'freezing_curve': [0.2, 0.5, 0.5, 0.5, 0.5],  # freezing curve parameter
                      'bedrock': {
                                  'solid_heat_capacity': 2.16e6,  # [J m-3 (solid) K-1]
                                  'thermal_conductivity': 3.0  # thermal conductivity of non-porous bedrock [W m-1 K-1]
                                  }
                      }

    return soil_properties, zh