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
control = single_lad_profiles(grid, fdir + 'letto2014.txt', hs, plot=False, biomass_function='aleksis_combination')  #, biomass_function='marklund_mod')
partial = single_lad_profiles(grid, fdir + 'letto2016_partial.txt', hs, plot=False, biomass_function='aleksis_combination')  #, biomass_function='marklund_mod')
clearcut = single_lad_profiles(grid, fdir + 'letto2016_clearcut.txt', hs, plot=False, biomass_function='aleksis_combination')  #, biomass_function='marklund_mod')

alpha=0.3
root_depth=0.2

def get_parameters(scenario):
    # spefify as one values (same for all simulations) or tuple of length 'count'
    if scenario.upper() == 'ALL':
        lettosuo_parameters = {
            'count': 3,
            'general':{
                        'start_time' : "2009-10-01",
                        'end_time' : "2019-01-01"
                        },
            'canopy': {
#                    'ctr': { # speed up!
#                            'WMA': True,  # well-mixed assumption
#                            'Ebal': False,  # no energy balance
#                            },
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
                            'wmax': 0.35e-03,
                            'wmaxsnow': 1.4e-03,
                            'c_rain': 1.05,
                            'c_snow': 1.3
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
                                            'ground_coverage': (0.4, 0.4, 0.0),
                                            }
                                    },
                            'litter': {
                                    'ground_coverage': (0.6, 0.6, 1.0)
                                    },
                            'baresoil': {
                                    'ground_coverage': (0.0, 0.0, 0.0)
                                    },
                            'snowpack': {
                                    'kmelt': 0.75*2.31e-8
                                    }
                            },
                    'planttypes': {
                            'pine': {
                                    'LAImax': (control['lai']['pine'], partial['lai']['pine'], clearcut['lai']['pine']),
                                    'lad': (control['lad']['pine'], partial['lad']['pine'], clearcut['lad']['pine']),
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 55.,
                                        'Jmax': 108.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.3,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'spruce': {
                                    'LAImax': (control['lai']['spruce'], partial['lai']['spruce'], clearcut['lai']['spruce']),
                                    'lad': (control['lad']['spruce'], partial['lad']['spruce'], clearcut['lad']['spruce']),
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 60.,
                                        'Jmax': 118.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.4,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'decid': {
                                    'LAImax': (control['lai']['decid'], partial['lai']['decid'], clearcut['lai']['decid']),
                                    'lad': (control['lad']['decid'], partial['lad']['decid'], clearcut['lad']['decid']),
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 4.0,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'shrubs': {
                                    'LAImax': (0.8, 0.8, 0.5),
                                    'lad': (control['lad']['shrubs'], partial['lad']['shrubs'], clearcut['lad']['shrubs']),
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
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
    #                        'type': 'Equilibrium',
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
        return lettosuo_parameters
    elif scenario.upper() == 'ALL2':
        lettosuo_parameters = {
            'count': 3,
            'general':{
                        'start_time' : "2009-10-01",
                        'end_time' : "2019-01-01"
                        },
            'canopy': {
#                    'ctr': { # speed up!
#                            'WMA': True,  # well-mixed assumption
#                            'Ebal': False,  # no energy balance
#                            },
                    'loc': {
                            'lat': 60.63,
                            'lon': 23.95
                            },
                    'micromet': {
                            'dPdx': 0.0
                            },
                    'radiation':{
                            'Par_alb': 0.12, # 0.11, # 0.1,
                            'Nir_alb': 0.55, # 0.46, # 0.43
                            },
                    'interception':{
                            'wmax': 0.35e-03,
                            'wmaxsnow': 1.4e-03,
                            'c_rain': 1.05,
                            'c_snow': 1.3
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
                                            'ground_coverage': (0.4, 0.4, 0.0),
                                            }
                                    },
                            'litter': {
                                    'ground_coverage': (0.6, 0.6, 1.0)
                                    },
                            'baresoil': {
                                    'ground_coverage': (0.0, 0.0, 0.0)
                                    },
                            'snowpack': {
                                    'kmelt': 0.75*2.31e-8
                                    }
                            },
                    'planttypes': {
                            'pine': {
                                    'LAImax': (control['lai']['pine'], partial['lai']['pine'], clearcut['lai']['pine']),
                                    'lad': (control['lad']['pine'], partial['lad']['pine'], clearcut['lad']['pine']),
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 55.,
                                        'Jmax': 108.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.3,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 4.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'spruce': {
                                    'LAImax': (control['lai']['spruce'], partial['lai']['spruce'], clearcut['lai']['spruce']),
                                    'lad': (control['lad']['spruce'], partial['lad']['spruce'], clearcut['lad']['spruce']),
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 60.,
                                        'Jmax': 118.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.4,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 4.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'decid': {
                                    'LAImax': (control['lai']['decid'], partial['lai']['decid'], clearcut['lai']['decid']),
                                    'lad': (control['lad']['decid'], partial['lad']['decid'], clearcut['lad']['decid']),
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 4.5,
                                        'g0': 1.0e-2,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'shrubs': {
                                    'LAImax': (0.8, 0.8, 0.5),
                                    'lad': (control['lad']['shrubs'], partial['lad']['shrubs'], clearcut['lad']['shrubs']),
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 40.,
                                        'Jmax': 79.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 0.9,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 4.5,
                                        'g0': 1.0e-2,  # this needs to be small, otherwise tr during dry conditions too high..
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
    #                        'type': 'Equilibrium',
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
        return lettosuo_parameters
    elif scenario.upper() == 'CLEARCUT':
        lettosuo_parameters = {
            'count': 4,
            'general':{
                        'start_time' : "2015-10-01",
                        'end_time' : "2019-01-01"
                        },
            'canopy': {
                    'ctr': { # speed up!
                            'WMA': True,  # well-mixed assumption
                            'Ebal': False,  # no energy balance
                            },
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
                            'wmax': 0.35e-03,
                            'wmaxsnow': 1.4e-03,
                            'c_rain': 1.05,
                            'c_snow': 1.3
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
                                            'ground_coverage': 0.0,
                                            }
                                    },
                            'litter': {
                                    'ground_coverage': 1.0,
                                    },
                            'baresoil': {
                                    'ground_coverage': 0.0,
                                    },
                            'snowpack': {
                                    'kmelt': 0.75*2.31e-8
                                    }
                            },
                    'planttypes': {
                            'pine': {
                                    'LAImax': clearcut['lai']['pine'],
                                    'lad': clearcut['lad']['pine'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 55.,
                                        'Jmax': 108.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.3,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'spruce': {
                                    'LAImax': clearcut['lai']['spruce'],
                                    'lad': clearcut['lad']['spruce'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 60.,
                                        'Jmax': 118.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.4,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'decid': {
                                    'LAImax': clearcut['lai']['decid'],
                                    'lad': clearcut['lad']['decid'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 4.0,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'shrubs': {
                                    'LAImax': (0.2, 0.4, 0.8, 1.0),
                                    'lad': clearcut['lad']['shrubs'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
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
    #                        'type': 'Equilibrium',
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
                                   },
                            },
                    }
            }
        return lettosuo_parameters
    elif scenario.upper() == 'PARTIAL':
        lettosuo_parameters = {
            'count': 3,
            'general':{
                        'start_time' : "2015-10-01",
                        'end_time' : "2019-01-01"
                        },
            'canopy': {
                    'ctr': { # speed up!
                            'WMA': True,  # well-mixed assumption
                            'Ebal': False,  # no energy balance
                            },
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
                            'wmax': 0.35e-03,
                            'wmaxsnow': 1.4e-03,
                            'c_rain': 1.05,
                            'c_snow': 1.3
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
                                            'ground_coverage': 0.4,
                                            }
                                    },
                            'litter': {
                                    'ground_coverage': 0.6,
                                    },
                            'baresoil': {
                                    'ground_coverage': 0.0,
                                    },
                            'snowpack': {
                                    'kmelt': 0.75*2.31e-8
                                    }
                            },
                    'planttypes': {
                            'pine': {
                                    'LAImax': partial['lai']['pine'],
                                    'lad': partial['lad']['pine'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 55.,
                                        'Jmax': 108.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.3,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'spruce': {
                                    'LAImax': partial['lai']['spruce'],
                                    'lad': partial['lad']['spruce'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 60.,
                                        'Jmax': 118.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.4,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'decid': {
                                    'LAImax': partial['lai']['decid'],
                                    'lad': partial['lad']['decid'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 4.0,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'shrubs': {
                                    'LAImax': (0.4, 0.8, 1.2),
                                    'lad': partial['lad']['shrubs'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
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
    #                        'type': 'Equilibrium',
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
                                   },
                            },
                    }
            }
        return lettosuo_parameters
    if scenario.upper() == 'CONTROL':
        lettosuo_parameters = {
            'count': 3,
            'general':{
                        'start_time' : "2009-10-01",
                        'end_time' : "2019-01-01"
                        },
            'canopy': {
#                    'ctr': { # speed up!
#                            'WMA': True,  # well-mixed assumption
#                            'Ebal': False,  # no energy balance
#                            },
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
                            'wmax': 0.35e-03,
                            'wmaxsnow': 1.4e-03,
                            'c_rain': 1.05,
                            'c_snow': 1.3
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
                                            'ground_coverage': 0.4,
                                            }
                                    },
                            'litter': {
                                    'ground_coverage': 0.6,
                                    },
                            'baresoil': {
                                    'ground_coverage': 0.0,
                                    },
                            'snowpack': {
                                    'kmelt': 0.75*2.31e-8
                                    }
                            },
                    'planttypes': {
                            'pine': {
                                    'LAImax': control['lai']['pine'],
                                    'lad': control['lad']['pine'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 55.,
                                        'Jmax': 108.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.3,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'spruce': {
                                    'LAImax': control['lai']['spruce'],
                                    'lad': control['lad']['spruce'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 60.,
                                        'Jmax': 118.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.4,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 2.5,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'decid': {
                                    'LAImax': control['lai']['decid'],
                                    'lad': control['lad']['decid'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
                                        'alpha': alpha,
                                        'theta': 0.7,
                                        'g1': 4.0,
                                        'g0': 1.0e-3,  # this needs to be small, otherwise tr during dry conditions too high..
                                        },
                                    },
                            'shrubs': {
                                    'LAImax': (0.5, 0.8, 1.1),
                                    'lad': control['lad']['shrubs'],
                                    'rootp': {
                                            'root_depth': root_depth,
                                            },
                                    'photop': {
                                        'Vcmax': 45.,
                                        'Jmax': 89.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
                                        'Rd': 1.0,  # 0.023*Vcmax
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
    #                        'type': 'Equilibrium',
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
        return lettosuo_parameters
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
