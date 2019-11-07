#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

ranges = {}

def get_parameters(scenario):
    # spefify as one values (same for all simulations) or tuple of length 'count'
    if scenario.upper() == 'CASE_1':
        parameters = {
                'count': 1,
                'general':{
                        'start_time' : "2005-05-01",
                        'end_time' : "2005-11-01"
                        },
                'canopy': {
                        'ctr': { # controls
                            'WMA': True,  # well-mixed assumption
                            'Ebal': False,  # no energy balance
                                },
                        },
                'soil': {
                    'water_model': {
                            'solve': False,
#                            'type': 'Equilibrium',
                            },
                    'heat_model':{
                            'solve': False,
                            },
                        },
                    }
        return parameters
    else:
        raise ValueError("Unknown parameterset!")

def iterate_parameters(parameters, default, count):
    """ Going through recursively senstivity nested parameter dictionary.
    Args:
        paramters (dict): nested dictionary
    """

    for key, value in parameters.items():
        if key == 'count':
            continue
        elif isinstance(value, dict):
            if key in default:
                default[key] = iterate_parameters(value, default[key], count)
            else:
                print(key + ' not in parameters')

        else:
            if isinstance(value, tuple):
                default[key] = value[count]
            else:
                default[key] = value

    return default
