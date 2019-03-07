#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

from pyAPES_utilities.soiltypes.for_morris import soil_organic

ranges = {}

parameters = {
        'count': 1,
        'canopy': {
            'forestfloor': {
                        'bryophytes': {
                                'hylocomium': {
                                        'ground_coverage': 1.0,
                                        },
                                'sphagnum': {
                                        'ground_coverage': 0.0,
                                        },
                                'pleurozium': {
                                        'ground_coverage': 0.0,
                                        }
                                },
                        'litter': {
                                'ground_coverage': 0.0
                                },
                        'baresoil': {
                                'ground_coverage': 0.0
                                }
                        },
            'planttypes': {
                        'pine': {
                                'LAImax': 5.0
                                },
                        'spruce': {
                                'LAImax': 0.0
                                },
                        'decidious': {
                                'LAImax': 0.0
                                },
                        'shrubs': {
                                'LAImax': 0.5
                                },
                        },
                },
        'soil': {
                'grid': {
                        'zh': soil_organic['zh']
                        },
                'soil_properties': soil_organic['properties'],
                'water_model': {
                        'initial_condition':{
                                'ground_water_level': -0.2
                                },
                        'lower_boundary': {  # lower boundary condition (type, value, depth)
                               'type': 'impermeable',
                               'value': None,
                               'depth': -2.0
                               },
                       'drainage_equation': {  # drainage equation and drainage parameters
                               'type': 'Hooghoudt',  #
                               'depth': 0.2,  # drain depth [m]
                               'spacing': 200.0,  # drain spacing [m]
                               'width': 1.0,  # drain width [m]
                               }
                        }
                }
        }

parameters2 = {
        'count': 1,
        'canopy': {
            'forestfloor': {
                        'bryophytes': {
                                'hylocomium': {
                                        'ground_coverage': 1.0,
                                        },
                                'sphagnum': {
                                        'ground_coverage': 0.0,
                                        },
                                'pleurozium': {
                                        'ground_coverage': 0.0,
                                        }
                                },
                        'litter': {
                                'ground_coverage': 0.0
                                },
                        'baresoil': {
                                'ground_coverage': 0.0
                                }
                        },
            'planttypes': {
                        'pine': {
                                'LAImax': 5.0
                                },
                        'spruce': {
                                'LAImax': 0.0
                                },
                        'decidious': {
                                'LAImax': 0.0
                                },
                        'shrubs': {
                                'LAImax': 0.5
                                },
                        },
                },
        }

def iterate_parameters(parameters, default, count):
    """ Going through recursively senstivity nested parameter dictionary.

    Args:
        paramters (dict): nested dictionary
    """

    for key, value in parameters.items():
        if key == 'count':
            continue
        elif isinstance(value, dict):
            default[key] = iterate_parameters(value, default[key], count)

        else:
            if isinstance(value, tuple):
                default[key] = value[count]
            else:
                default[key] = value

    return default
