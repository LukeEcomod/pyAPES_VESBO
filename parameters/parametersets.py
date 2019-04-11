#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

from pyAPES_utilities.parameter_utilities import single_lad_profiles
from pyAPES_utilities.soiltypes.organic import soil_properties, zh
from parameters.canopy import grid

ranges = {}

# normed leaf area density profiles
fdir = 'pyAPES_utilities/runkolukusarjat/'
hs = 0.5  # height of understory shrubs [m]
control = single_lad_profiles(grid, fdir + 'letto2014.txt', hs, plot=False, biomass_function='marklund_mod')
partial = single_lad_profiles(grid, fdir + 'letto2016_partial.txt', hs, plot=False, biomass_function='marklund_mod')
clearcut = single_lad_profiles(grid, fdir + 'letto2016_clearcut.txt', hs, plot=False, biomass_function='marklund_mod')

# spefify as one values (same for all simulations) or tuple of length 'count'
lettosuo_parameters = {
        'count': 3,
        'canopy': {
                'loc': {
                        'lat': 60.63,
                        'lon': 23.95
                        },
                'micromet': {
                        'dPdx': 0.0
                        },
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
                                'ground_coverage': (0.5, 0.5, 0.9)
                                },
                        'baresoil': {
                                'ground_coverage': (0.0, 0.0, 0.0)
                                }
                        },
                'planttypes': {
                        'pine': {
                                'LAImax': (control['lai']['pine'], partial['lai']['pine'], clearcut['lai']['pine']),
                                'lad': (control['lad']['pine'], partial['lad']['pine'], clearcut['lad']['pine']),
                                'rootp': {
                                        'root_depth': 0.2
                                        }
                                    },
                        'spruce': {
                                'LAImax': (control['lai']['spruce'], partial['lai']['spruce'], clearcut['lai']['spruce']),
                                'lad': (control['lad']['spruce'], partial['lad']['spruce'], clearcut['lad']['spruce']),
                                'rootp': {
                                        'root_depth': 0.2
                                        }
                                    },
                        'decidious': {
                                'LAImax': (control['lai']['decid'], partial['lai']['decid'], clearcut['lai']['decid']),
                                'lad': (control['lad']['decid'], partial['lad']['decid'], clearcut['lad']['decid']),
                                'rootp': {
                                        'root_depth': 0.2
                                        }
                                    },
                        'shrubs': {
                                'LAImax': (0.8, 0.8, 0.5),
                                'lad': (control['lad']['shrubs'], partial['lad']['shrubs'], clearcut['lad']['shrubs']),
                                'rootp': {
                                        'root_depth': 0.2
                                        }
                                    },
                        },
                },
        'soil': {
                'grid': {
                        'zh': zh
                        },
                'soil_properties': soil_properties,
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
                               'depth': 0.8,  # drain depth [m]
                               'spacing': 45.0,  # drain spacing [m]
                               'width': 1.0,  # drain width [m]
                               }
                        }
                }
        }

lettosuo_parameters_clc = {
        'count': 2,
        'canopy': {
                'loc': {
                        'lat': 60.63,
                        'lon': 23.95
                        },
                'micromet': {
                        'dPdx': 0.0,
                        },
                'forestfloor': {
                        'bryophytes': {
                                'hylocomium': {
                                        'ground_coverage': 0.0,
                                        },
                                'sphagnum': {
                                        'ground_coverage': 0.0,
                                        },
                                'pleurozium': {
                                        'ground_coverage': 0.1,
                                        }
                                },
                        'litter': {
                                'ground_coverage': 0.9
                                },
                        'baresoil': {
                                'ground_coverage': 0.0
                                }
                        },
                'planttypes': {
                        'pine': {
                                'LAImax': clearcut['lai']['pine'],
                                'lad': clearcut['lad']['pine'],
                                'rootp': {
                                        'root_depth': 0.2
                                        }
                                    },
                        'spruce': {
                                'LAImax': clearcut['lai']['spruce'],
                                'lad': clearcut['lad']['spruce'],
                                'rootp': {
                                        'root_depth': 0.2
                                        }
                                    },
                        'decidious': {
                                'LAImax': clearcut['lai']['decid'],
                                'lad': clearcut['lad']['decid'],
                                'rootp': {
                                        'root_depth': 0.2
                                        }
                                    },
                        'shrubs': {
                                'LAImax': 0.5,
                                'lad': clearcut['lad']['shrubs'],
                                'rootp': {
                                        'root_depth': 0.2
                                        }
                                    },
                        },
                },
        'soil': {
                'grid': {
                        'zh': zh
                        },
                'soil_properties': soil_properties,
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
                               'depth': 0.8,  # drain depth [m]
                               'spacing': 45.0,  # drain spacing [m]
                               'width': 1.0,  # drain width [m]
                               }
                        }
                }
        }

alpha=0.3

lettosuo_parameters_ctrl = {
        'count': 1,
        'canopy': {
#                'ctr': { # speed up!
#                        'WMA': True,  # well-mixed assumption
#                        'Ebal': False,  # no energy balance
#                        },
                'loc': {
                        'lat': 60.63,
                        'lon': 23.95
                        },
                'micromet': {
                        'dPdx': 0.0
                        },
                'radiation':{
                        'Par_alb': 0.1,
                        'Nir_alb': 0.43
                        },
                'interception':{
                        'wmaxsnow': 0.75e-3, #1.6e-03, 0.75e-3, 0.75e-3),  # maximum interception storage capacity for snow [m per unit of LAI]
                        'Tmin': 0.0,  # temperature below which all is snow [degC]
                        'Tmax': 2.0, #2.0, 1.0, 2.0),  # temperature above which all is water [degC]
                        },
                'forestfloor': {
                        'bryophytes': {
                                'hylocomium': {
                                        'ground_coverage': 0.0,
                                        },
                                'sphagnum': {
                                        'ground_coverage': 0.0,
                                        },
                                'pleurozium': {
                                        'ground_coverage': 0.5,
                                        }
                                },
                        'litter': {
                                'ground_coverage': 0.5
                                },
                        'baresoil': {
                                'ground_coverage': 0.0
                                }
                        },
                'planttypes': {
                        'pine': {
                                'LAImax': control['lai']['pine'],
                                'lad': control['lad']['pine'],
                                'rootp': {
                                        'root_depth': 0.2
                                        },
                                'photop': {
                                    'Vcmax': 55.,
                                    'Jmax': 108.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                    'Rd': 1.3,  # 0.023*Vcmax
                                    'tresp': {
                                        'Vcmax': [72., 200., 649.],  # (Kattge and Knorr, 2007)
                                        'Jmax': [50., 200., 646.],  # (Kattge and Knorr, 2007)
                                        'Rd': [33.0]
                                        },
                                    'alpha': alpha,
                                    'theta': 0.7,
                                    'g1': 2.3,
                                    'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                    },
                                    },
                        'spruce': {
                                'LAImax': control['lai']['spruce'],
                                'lad': control['lad']['spruce'],
                                'rootp': {
                                        'root_depth': 0.2
                                        },
                                'photop': {
                                    'Vcmax': 60.,
                                    'Jmax': 118.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                    'Rd': 1.4,  # 0.023*Vcmax
                                    'tresp': {
                                        'Vcmax': [72., 200., 649.],  # (Kattge and Knorr, 2007)
                                        'Jmax': [50., 200., 646.],  # (Kattge and Knorr, 2007)
                                        'Rd': [33.0]
                                        },
                                    'alpha': alpha,
                                    'theta': 0.7,
                                    'g1': 2.3,
                                    'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                    },
                                    },
                        'decidious': {
                                'LAImax': control['lai']['decid'],
                                'lad': control['lad']['decid'],
                                'rootp': {
                                        'root_depth': 0.2
                                        },
                                'photop': {
                                    'Vcmax': 45.,
                                    'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                    'Rd': 1.0,  # 0.023*Vcmax
                                    'tresp': {
                                        'Vcmax': [72., 200., 649.],  # (Kattge and Knorr, 2007)
                                        'Jmax': [50., 200., 646.],  # (Kattge and Knorr, 2007)
                                        'Rd': [33.0]
                                        },
                                    'alpha': alpha,
                                    'theta': 0.7,
                                    'g1': 4.0,
                                    'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                    },
                                    },
                        'shrubs': {
                                'LAImax': 0.8,
                                'lad': control['lad']['shrubs'],
                                'rootp': {
                                        'root_depth': 0.2
                                        },
                                'photop': {
                                    'Vcmax': 45.,
                                    'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                    'Rd': 1.0,  # 0.023*Vcmax
                                    'tresp': {
                                        'Vcmax': [72., 200., 649.],  # (Kattge and Knorr, 2007)
                                        'Jmax': [50., 200., 646.],  # (Kattge and Knorr, 2007)
                                        'Rd': [33.0]
                                        },
                                    'alpha': alpha,
                                    'theta': 0.7,
                                    'g1': 4.0,
                                    'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                    },
                                    },
                        },
                },
        'soil': {
                'grid': {
                        'zh': zh
                        },
                'soil_properties': soil_properties,
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
                               'depth': 0.8,  # drain depth [m]
                               'spacing': 45.0,  # drain spacing [m]
                               'width': 1.0,  # drain width [m]
                               }
                        }
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
