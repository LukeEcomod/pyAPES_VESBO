#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

ranges = {}


# parameter sets of simulations
parameters = {
        'count': 4,
        'canopy': {
#                'forestfloor': {
#                        'bryophytes': {
#                                'hylocomium': {
#                                        'ground_coverage': (1.0, 1.0),
#                                        'height': (0.05, 0.05)
#                                        },
#                                'pleurozium': {
#                                        'ground_coverage': (0.0, 0.0),
#                                        'height': (0.04, 0.04)
#                                        }
#                                        },
#                        'baresoil': {
#                                'ground_coverage': (0.0, 0.0)
#                                }
#                                },
                'planttypes': {
                        'pine': {
                                'LAImax': ([2.1], [2.1], [1.5], [1.5])
                                },
                        'shrubs': {
                                'LAImax': ([0.6], [0.4], [0.2], [0.0])
                                }
                        }
        }
    }


def iterate_parameters(parameters, default, count):
    """ Going through recursively senstivity nested parameter dictionary.

    Args:
        paramters (dict): nested dictionary
    """

    for key, value in parameters.iteritems():
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