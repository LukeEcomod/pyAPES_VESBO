#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

ranges = {}


# parameter sets of simulations
parameters = {
        'canopy': {
                'forestfloor': {
                        'bryophytes': {
                                'hylocomium': {
                                        'ground_coverage': [0.5, 0.4],
                                        'height': [0.05, 0.05]
                                        },
                                'pleurozium': {
                                        'ground_coverage': [0.4, 0.5],
                                        'height': [0.04, 0.04]
                                        }
                                        },
                        'bareground': {
                                'ground_coverage': [0.1, 0.1]
                                }
                                },
                'planttypes': {
                        'pine': {
                                'LAImax': [[2.1], [2.1]]
                                },
                        'shrubs': {
                                'LAImax': [[0.7], [0.7]]
                                }
                        }
        }
    }


def iterate_parameters(parameters, param_space):
    """ Going through recursively senstivity nested parameter dictionary.

    Args:
        paramters (dict): nested dictionary
    """

    for key, value in parameters.iteritems():
        if isinstance(value, dict):
            iterate_parameters(value, param_space)
        else:
            