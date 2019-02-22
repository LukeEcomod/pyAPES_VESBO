#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

from pyAPES_utilities.parameter_utilities import lad_profiles
from pyAPES_utilities.soiltypes.organic import soil_properties, zh
from parameters.canopy import grid

ranges = {}

# normed leaf area density profiles
fdir = 'pyAPES_utilities/runkolukusarjat/'
quantiles = [1.0]  # quantiles used in creating tree lad profiles (groups of same planttype)
hs = 0.5  # height of understory shrubs [m]
control = lad_profiles(grid, fdir + 'letto2014.txt', quantiles, hs, plot=False)
partial = lad_profiles(grid, fdir + 'letto2016_partial.txt', quantiles, hs, plot=False)
clearcut = lad_profiles(grid, fdir + 'letto2016_clearcut.txt', quantiles, hs, plot=False)

# spefify as one values (same for all simulations) or tuple of length 'count'
lettosuo_parameters = {
        'count': 3,
        'canopy': {
                'forestfloor': {
                        'bryophytes': {
                                'hylocomium': {
                                        'ground_coverage': 0.0,
                                        },
                                'sphagnum': {
                                        'ground_coverage': 0.0,
                                        },
                                'pleurozium': {
                                        'ground_coverage': (0.5, 0.5, 0.1),
                                        }
                                },
                        'litter': {
                                'ground_coverage': (0.5, 0.5, 0.7)
                                },
                        'baresoil': {
                                'ground_coverage': (0.0, 0.0, 0.2)
                                }
                        },
                'planttypes': {
                        'pine': {
                                'LAImax': (control['lai']['pine'], partial['lai']['pine'], clearcut['lai']['pine']),
                                'lad': (control['lad']['pine'], partial['lad']['pine'], clearcut['lad']['pine'])
                                    },
                        'spruce': {
                                'LAImax': (control['lai']['spruce'], partial['lai']['spruce'], clearcut['lai']['spruce']),
                                'lad': (control['lad']['spruce'], partial['lad']['spruce'], clearcut['lad']['spruce'])
                                    },
                        'decidious': {
                                'LAImax': (control['lai']['decid'], partial['lai']['decid'], clearcut['lai']['decid']),
                                'lad': (control['lad']['decid'], partial['lad']['decid'], clearcut['lad']['decid'])
                                    },
                        'shrubs': {
                                'LAImax': ([0.8], [0.8], [0.5]),
                                'lad': (control['lad']['shrubs'], partial['lad']['shrubs'], clearcut['lad']['shrubs'])
                                    },
                        },
                },
        'soil': {
                'grid': {
                        'zh': zh
                        },
                'soil_properties': soil_properties
                }
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
