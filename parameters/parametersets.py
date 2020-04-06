#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

# Previously used parameter set based on drained organic soil in Lettosuo
#from pyAPES_utilities.soiltypes.organic import soil_properties, zh
# Revised parameter set based on properties on sphagnum peat used in SpaFHy_peat (Kersti Haahti)
from pyAPES_utilities.soiltypes.sphagnum_peat import soil_properties, zh
from tools.utilities import lad_constant
from parameters.canopy import z

def get_parameters(scenario):
    # spefify as one values (same for all simulations) or tuple of length 'count'
    if scenario.upper() == 'DEGERO':
        parameters = {
            'count': 1,
            'scenario': 'Degero',
            'general':{
                'start_time' : "2014-06-01",
                'end_time' : "2014-06-05",
                'forc_filename' : "Degero_forcing_2014_2016.csv"
            },
            'canopy': {
                # 'ctr': { # controls
                #     'WMA': True,  # well-mixed assumption
                #     'Ebal': False,  # no energy balance
                #   },
                'planttypes': {
                    'plant1': {
                        'LAImax': 0.6,
                        'lad': lad_constant(z, 1.0, 0.0)
                    },
                    'plant2': {
                        'LAImax': 0.0
                    },
                },
                'forestfloor': {
                    'bottom_layer_types': {
                        'litter': {
                            'coverage': 0.0
                        },
                        'forest_moss': {
                            'coverage': 0.0
                        },
                        'sphagnum': {
                            'coverage': 1.0,
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
                        # 'type': 'Equilibrium',
                        'initial_condition':{
                            'ground_water_level': -0.05
                        },
                    'lower_boundary': {  # lower boundary condition (type, value, depth)
                    'type': 'impermeable',
                        'value': None,
                        'depth': -2.0
                    },
                    'drainage_equation': {  # drainage equation and drainage parameters
                        'type': 'Hooghoudt',  #
                        'depth': 0.0,  # drain depth [m]
                        'spacing': 100.0,  # drain spacing [m]
                        'width': 1.0,  # drain width [m]
                    }
                }
            }
        }

        return parameters

    else:
        raise ValueError("Unknown parameterset!")

# EOF