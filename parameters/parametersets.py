#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

# Modifications to some parameters

def get_parameters(scenario):
    # spefify as one values (same for all simulations) or tuple of length 'count'
    if scenario.upper() == 'TEST':
        parameters = {
            'count': 2,
            'scenario': 'test',
            'general':{
                'start_time' : "2005-05-01",
                'end_time' : "2005-09-30",
                'forc_filename' : "Hyytiala/FIHy_forcing_2005-2010.dat"
            },
            'canopy': {
                'planttypes': {
                    'decid': {
                        'LAImax': (1.2, 3.0),
                    },
                },
            },
        }
        return parameters
    elif scenario.upper() == 'BYPASS_SOIL':
        parameters = {
            'count': 1,
            'scenario': 'bypass_soil',
            'general':{
                'start_time' : "2005-06-01",
                'end_time' : "2005-06-05",
                'forc_filename' : "Hyytiala/FIHy_forcing_2005-2010.dat"
            },
            'soil':{  # one layer of 10 cm
                'grid': {
                    'dz': [0.1],
                    'zh': [-0.1]
                    },
                'soil_properties': {
                    'pF': {
                        'ThetaS': [0.80],
                        'ThetaR': [0.01],
                        'alpha': [0.70],
                        'n': [1.25]
                        },
                    'saturated_conductivity_vertical': [2.42E-05],
                    'saturated_conductivity_horizontal': [2.42E-05],
                    'solid_composition': {
                         'organic': [0.1611],
                         'sand': [0.4743],
                         'silt': [0.3429],
                         'clay': [0.0217]
                         },
                    'freezing_curve': [0.2]
                    },
                'water_model':{
                    'solve': False
                    },
                'heat_model':{
                    'solve': False
                    },
                }
        }
        return parameters

    else:
        raise ValueError("Unknown parameterset!")

# EOF