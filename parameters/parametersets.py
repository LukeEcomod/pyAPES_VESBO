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
                'start_time' : "2005-06-01",
                'end_time' : "2005-06-10",
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

    elif scenario.upper() == 'KRYCKLAN':
        from pyAPES_utilities.parameter_utilities import single_lad_profiles
        from parameters.SmearII import grid
        # normed leaf area density profiles
        fdir = 'pyAPES_utilities/runkolukusarjat/'
        hs = 0.5  # height of understory shrubs [m]
        stand = single_lad_profiles(grid, fdir + 'Krycklan_c2.txt', hs, plot=False, biomass_function='marklund')

        parameters = {
                'count': 1,
                'scenario': 'krycklan',
                'general':{
                    'start_time' : "2005-06-01",
                    'end_time' : "2005-06-10",
                    'forc_filename' : "Hyytiala/FIHy_forcing_2005-2010.dat"
                },
                'canopy': {
                        'loc': {
                                'lat': 64.26,
                                'lon': 19.77
                                },
                        'radiation':{
                                'Par_alb': 0.1,
                                'Nir_alb': 0.43,
                                },
                        'interception': {
                                'wmax': 0.35
                                },
                        'planttypes': {
                                'pine': {
                                        'LAImax': 0.31 * 4.8,
                                        'lad': stand['lad']['pine']
                                        },
                                'spruce': {
                                        'LAImax': 0.64 * 4.8,
                                        'lad': stand['lad']['spruce']
                                        },
                                'decid': {
                                        'LAImax': 0.05 * 4.8,
                                        'lad': stand['lad']['decid']
                                        },
                                'shrubs': {
                                        'LAImax': 0.6,
                                        'lad': stand['lad']['shrubs']
                                        }
                                }
                }
            }
    else:
        raise ValueError("Unknown parameterset!")

    return parameters

# EOF