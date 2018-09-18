# -*- coding: utf-8 -*-
"""
SOIL DISCRETIZATION AND PARAMETERS
"""
import numpy as np
from parameter_utils import fit_pF, peat_hydrol_properties

# thickness of layers with different characteristics
thickness = [0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 1.0]
# number of nodes in each layer
nodes = [10, 5, 5, 2, 2, 5, 5]

# ---- Layer characteristics ----

# Water retention
## PÄIVÄNEN SEDGE SATTUMANVARAISESTI
#watcont = [[94.3, 68, 47.9, 35.5, 22, 18.4, 16.7, 12.6, 7.9, 6.3, -999],
#           [91.7, 82.9, 61.8, 35.9, 31.7, 25.2, 23.4, 19.2, 17.3, 14.4, -999],
#           [90.6, 86.2, 56.4, 36.4, 33.4, 29.8, 26.8, 23.5, 20.2, 16.4, -999],
#           [89.7, 85, 74.5, 53.4, 36, 29, 24.7, 22.1, 17.6, 14.7, -999],
#           [87.3, 85.4, 77.8, 64, 41.4, 28.8, 23.4, 22.7, 21.9, 17, -999],
#           [89.3, 86.5, 80.7, 52.5, 45.6, 35.4, 32, 25.1, 20.6, 18.4, -999],
#           [91, 89.9, 84.7, 60, -999, 33.8, 27.2, 29.3, -999, 17.2, 12.9]]
#head = [0.0001, 1, 3.2, 10, 20, 60, 100, 200, 500, 1000, 1500]
#pF_para = fit_pF(head, watcont, fig=False)
## ------------ LAIHO 2015 ------------
#watcont = [[94.69, 49.42, 29.61, 21.56, 20.05, 17.83, 16.54],
#           [91.41, 66.26, 56.98, 45.58, 41.44, 39.32, 37.89],
#           [89.12, -999, 72.83, 63.97, 54.40, 50.15, 48.80],
#           [89.46, -999, 82.46, 76.79, 66.93, 63.61, 62.53],
#           [92.22, -999, 87.06, 78.02, 74.76, 72.77, 71.70],
#           [92.22, -999, 87.06, 78.02, 74.76, 72.77, 71.70],
#           [92.22, -999, 87.06, 78.02, 74.76, 72.77, 71.70]]
#head = [0.0001, 0.3, 0.981, 4.905, 9.81, 33.0, 98.1]
## -------------------------------------
## Fit water retention parameters
#pF_para = fit_pF(head, watcont, fig=False)

# pF based on bulk density
bd = [0.0811, 0.140, 0.174, 0.171, 0.119, 0.119, 0.119]
pF_para, _ = peat_hydrol_properties(bd)#, fig=True)

# Hydraulic conductivity [m s-1]
Kvsat = [1.7e-4, 2e-5, 5e-5, 5e-6, 3e-6, 1e-6, 1e-7]  # vertical
Khmult = [10.0, 10.0, 5.0, 5.0, 5.0, 1.0, 1.0]  # horizontal

# Save to arrays
z = np.array([])
ThetaS = np.array([])
ThetaR = np.array([])
alpha = np.array([])
n = np.array([])

Ksat = np.array(Kvsat)
Khsat = np.array(Kvsat)*np.array(Khmult)

zh = np.array([0.0])

for k in range(0, len(thickness)):
    z = np.append(z, -np.arange(zh[k] + thickness[k]/(2 * nodes[k]), zh[k] + thickness[k], thickness[k] / nodes[k]))
    ThetaS = np.append(ThetaS, pF_para[k][0])
    ThetaR = np.append(ThetaR, pF_para[k][1])
    alpha = np.append(alpha, pF_para[k][2])
    n = np.append(n,  pF_para[k][3])
    zh = np.append(zh, zh[k] + thickness[k])

zh = -zh[1:]

# Soil model parameters
spara = {
        'z': z,
        'zh':zh,
        'pF': {
                'ThetaS': ThetaS, 
                'ThetaR': ThetaR, 
                'alpha': alpha,     # 1/cm
                'n': n
                },                  # (dict): vanGenuchten water retention parameters; scalars or arrays of len(z)
        'Ksat': Ksat,               # (float/array): saturated vertical hydraulic conductivity [ms-1]
        'Khsat': Khsat,              # (float/array): saturated horizontal hydraulic conductivity [ms-1] - used in drainage equation
        'Csv': 250000.0 * np.ones(len(zh)),                # (float/array): dry soil vol. heat capacity [J m-3 (total volume) K-1]
        'vOrg': 1.0 * np.ones(len(zh)),               # (float/array): organic matter fraction of solid volume [-]
        'vSand': 0.0 * np.ones(len(zh)),              # (float/array): sand fraction of solid volume [-]
        'vSilt': 0.0 * np.ones(len(zh)),              # (float(array): silt fraction of solid volume [-]
        'vClay': 0.0 * np.ones(len(zh)),              # (float(array): clay fraction of solid volume [-]
        'fp': 0.5 * np.ones(len(zh)),                 # (float/array): freezing curve parameter
        'max_pond': 0.01,           # (float) maximum pond depth [m]
        'ini_cond': {               # (dict): inputs are floats or arrays of len(z)
                'gwl': -0.2,        # (float) [m] or (float/array) Wtot', vol. water content [-] or 'h', matrix water potential [m]
                'T': 4.0,          # soil temperature [degC]
                'pond': 0.0
                },                  # [m] initial pond depth at surface
        'lbc_heat': {
                'type': 'temperature',
                'value': 4.0
                },           # lower boundary condition type for heat
        'lbc_water': {
                'type': 'impermeable',
                'value': None,
                'depth': -2.0
                },                  # lower boundary condition type for water (type, value, depth)
        'homogenous': False,         # (boolean): True assumes vertically homogenous profile and float inputs
        'solve_heat': True,        # (boolean): True solves heatflow
        'solve_water': True,        # (boolean): True solves waterflow
        'solve_water_type':'Equilibrium', #'Richards', # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
        'Bedrock': {
                'Cv': 2160000.0,
                'Lambda': 3.0
                },
        'drainage_equation': {
                'type': 'Hooghoudt',
                'depth': 1.0,
                'spacing': 40.0,
                'width': 1.0,
                }                   # (dict): Drainage equation and drainage parameters
        }

"""
spara = {
        'z': z, # np.arange(-0.005, -2.0, -0.01),
        'pF': {
                'ThetaS': 0.88, 
                'ThetaR': 0.093, 
                'alpha': 0.029,     # 1/cm
                'n': 1.34
                },                  # (dict): vanGenuchten water retention parameters; scalars or arrays of len(z)
        'Ksat': 1e-5,               # (float/array): saturated vertical hydraulic conductivity [ms-1]
        'Khsat': 1e-5,              # (float/array): saturated horizontal hydraulic conductivity [ms-1] - used in drainage equation
        'Csv': -1.0,                # (float/array): dry soil vol. heat capacity [J m-3 (total volume) K-1]
        'vOrg': -1.0,               # (float/array): organic matter fraction of solid volume [-]
        'vSand': -1.0,              # (float/array): sand fraction of solid volume [-]
        'vSilt': -1.0,              # (float(array): silt fraction of solid volume [-]
        'vClay': -1.0,              # (float(array): clay fraction of solid volume [-]
        'fp': -1.0,                 # (float/array): freezing curve parameter
        'max_pond': 0.01,           # (float) maximum pond depth [m]
        'ini_cond': {               # (dict): inputs are floats or arrays of len(z)
                'gwl': -0.2,        # (float) [m] or (float/array) Wtot', vol. water content [-] or 'h', matrix water potential [m]
                'T': -1.0,          # soil temperature [degC]
                'pond': 0.0
                },                  # [m] initial pond depth at surface
        'lbc_heat': -1.0,           # lower boundary condition type for heat
        'lbc_water': {
                'type': 'impermeable',
                'value': None,
                'depth': -2.0,
                },                  # lower boundary condition type for water (type, value, depth)
        'homogenous': True,         # (boolean): True assumes vertically homogenous profile and float inputs
        'solve_heat': False,        # (boolean): True solves heatflow
        'solve_water': True,        # (boolean): True solves waterflow
        'solve_water_type':'Equilibrium',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
        'Bedrock': {
                'Cv': 2160000.0,
                'Lambda': 3.0
                },
        'drainage_equation': {
                'type': 'Hooghoudt',
                'depth': 1.0,
                'spacing': 40.0,
                'width': 1.0,
                }                   # (dict): Drainage equation and drainage parameters
        }
"""