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
                    'start_time' : "2019-07-01",
                    'end_time' : "2019-07-30",
                    'forc_filename' : "Svartberget/Svartberget_forcing_2019.csv" # "Svartberget/Svartberget_forcing_2014_2016.csv"
                },
                'canopy': {
                        'loc': {
                                'lat': 64.26,
                                'lon': 19.77
                                },
                        'radiation':{
                                'Par_alb': 0.1,
                                'Nir_alb': 0.39,
                                },
                        'interception': {
                                'wmax': 0.35,
                                'wmaxsnow': 0.7
                                },
                        'planttypes': {
                                'pine': {
                                        'LAImax': 0.31 * 4.8,
                                        'lad': stand['lad']['pine'],
                                        'phenop': {
                                            'smax': 15.0  # Kolari 2014
                                            },
                                        'photop': {
                                            'Vcmax': 50.,
                                            'Jmax': 98.,
                                            'Rd': 1.2,
                                            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                                                'Vcmax': [72., 200., 649.],
                                                'Jmax': [50., 200., 646.],
                                                },
                                            'g1': 2.5,
                                            'g0': 4.0e-3, #(-4.0e-3,1.0e-5,4.0e-3),
                                            },
                                        },
                                'spruce': {
                                        'LAImax': 0.64 * 4.8,
                                        'lad': stand['lad']['spruce'],
                                        'phenop': {
                                            'smax': 15.0  # Kolari 2014
                                            },
                                        'photop': {
                                            'Vcmax': 60.,
                                            'Jmax': 113.,
                                            'Rd': 1.4,
                                            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                                                'Vcmax': [72., 200., 649.],
                                                'Jmax': [50., 200., 646.],
                                                },
                                            'g1': 2.5,
                                            'g0': 4.0e-3, #(-4.0e-3,1.0e-5,4.0e-3),
                                            },
                                        },
                                'decid': {
                                        'LAImax': 0.05 * 4.8,
                                        'lad': stand['lad']['decid'],
                                        'phenop': {
                                            'smax': 15.0  # Kolari 2014
                                            },
                                        'photop': {
                                            'Vcmax': 45.,
                                            'Jmax': 89.,
                                            'Rd': 1.0,
                                            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                                                'Vcmax': [72., 200., 649.],
                                                'Jmax': [50., 200., 646.],
                                                },
                                            'g1': 4.5,
                                            'g0': 1.0e-2, #(-4.0e-3,1.0e-5,1.0e-2),
                                            },
                                        },
                                'shrubs': {
                                        'LAImax': 0.6,
                                        'lad': stand['lad']['shrubs'],
                                        'phenop': {
                                            'smax': 15.0  # Kolari 2014
                                            },
                                        'photop': {
                                            'Vcmax': 40.,
                                            'Jmax': 79.,
                                            'Rd': 0.9,
                                            'kn': 0.0,
                                            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                                                'Vcmax': [72., 200., 649.],
                                                'Jmax': [50., 200., 646.],
                                                },
                                            'g1': 4.5,
                                            'g0': 1.0e-2, #(-4.0e-3,1.0e-5,1.0e-2),
                                            },
                                        }
                                },
                        'forestfloor': {
                                'snowpack': {
                                    'initial_conditions': {
                                        'snow_water_equivalent':20.0
                                        }
                                    }
                                }
                }
            }
    else:
        raise ValueError("Unknown parameterset!")

    return parameters

# EOF