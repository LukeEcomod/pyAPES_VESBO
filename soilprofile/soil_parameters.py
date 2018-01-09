# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 14:05:24 2017

@author: L1656
"""
import numpy as np

para = {
        'z': np.arange(-0.01, -2.0, -0.02),
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
        'solve_water': True,        # (boolean): True solves waterflow
        'solve_water_type': '',  # 'Equilibrium',  # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
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


        