#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

ranges = {}

parameters = {
        'count': 1,
        'canopy': {
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
