#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

ranges = {}

def get_parameters(name=None):

    if name is None:
        # parameter sets of simulations
        parameters = {
                'count': 1,
                'canopy': {
                        'forestfloor': {
                                'bryophytes': {
                                        'hylocomium': {
                                                'ground_coverage': (1.0)
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.0)
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (0.0)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.0)
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.0)
                                        }
                                },
                        'planttypes': {
                                'pine': {
                                        'LAImax': ([0.3 * 4.1])
                                        },
                                'spruce': {
                                        'LAImax': ([0.65 * 4.1])
                                        },
                                'decidious': {
                                        'LAImax': ([0.05 * 4.1])
                                        },
                                'shrubs': {
                                        'LAImax': ([0.5])
                                        }
                                }
                }
            }

    else:
        raise ValueError('Unknown sensitivity parametrization %s' % name.upper())

    return parameters


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
