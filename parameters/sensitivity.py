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
                        'interception': {
                                'wmax': (0.35e-03)
                                },
                        'forestfloor': {
                                'bryophytes': {
                                        'hylocomium': {
                                                'ground_coverage': (0.0)
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.0)
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (1.0)
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
                                        'LAImax': ([0.31 * 5.15])
                                        },
                                'spruce': {
                                        'LAImax': ([0.64 * 5.15])
                                        },
                                'decidious': {
                                        'LAImax': ([0.05 * 5.15])
                                        },
                                'shrubs': {
                                        'LAImax': ([0.6])
                                        }
                                }
                }
            }
    elif name=='Krycklan_sensitivity':
        N = 27
        LAI_tot = ([5.15], [5.15], [5.15], [5.15], [5.15], [5.15], [5.15], [5.15], [5.15], [3.4], [3.4], [3.4], [3.4], [3.4], [3.4], [3.4], [3.4], [3.4], [6.9], [6.9], [6.9], [6.9], [6.9], [6.9], [6.9], [6.9], [6.9])
        # parameter sets of simulations
        parameters = {
                'count': 27,
                'canopy': {
                        'interception': {
                                'wmax': (0.00035, 0.0002, 0.0005, 0.00035, 0.0002, 0.0005, 0.00035, 0.0002, 0.0005, 0.00035, 0.0002, 0.0005, 0.00035, 0.0002, 0.0005, 0.00035, 0.0002, 0.0005, 0.00035, 0.0002, 0.0005, 0.00035, 0.0002, 0.0005, 0.00035, 0.0002, 0.0005)
                                },
                        'forestfloor': {
                                'bryophytes': {
                                        'hylocomium': {
                                                'ground_coverage': (0.0,) * N
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.0,) * N
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (1.0,) * N
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.0,) * N
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.0,) * N
                                        }
                                },
                        'planttypes': {
                                'pine': {
                                        'LAImax': tuple([[0.31*LAI_tot[i][0]] for i in range(N)])
                                        },
                                'spruce': {
                                        'LAImax': tuple([[0.64*LAI_tot[i][0]] for i in range(N)])
                                        },
                                'decidious': {
                                        'LAImax': tuple([[0.05*LAI_tot[i][0]] for i in range(N)])
                                        },
                                'shrubs': {
                                        'LAImax': ([0.6], [0.6], [0.6], [0.4], [0.4], [0.4], [0.8], [0.8], [0.8], [0.6], [0.6], [0.6], [0.4], [0.4], [0.4], [0.8], [0.8], [0.8], [0.6], [0.6], [0.6], [0.4], [0.4], [0.4], [0.8], [0.8], [0.8])
                                        }
                                }
                }
            }
    elif name=='Krycklan_sensitivity2':
        N = 3
        LAI_tot = ([5.15], [3.4], [6.9])
        # parameter sets of simulations
        parameters = {
                'count': 3,
                'canopy': {
                        'interception': {
                                'wmax': (0.00035, 0.0002, 0.0005)
                                },
                        'forestfloor': {
                                'bryophytes': {
                                        'hylocomium': {
                                                'ground_coverage': (0.0,) * N
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.0,) * N
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (1.0,) * N
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.0,) * N
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.0,) * N
                                        }
                                },
                        'planttypes': {
                                'pine': {
                                        'LAImax': tuple([[0.31*LAI_tot[i][0]] for i in range(N)])
                                        },
                                'spruce': {
                                        'LAImax': tuple([[0.64*LAI_tot[i][0]] for i in range(N)])
                                        },
                                'decidious': {
                                        'LAImax': tuple([[0.05*LAI_tot[i][0]] for i in range(N)])
                                        },
                                'shrubs': {
                                        'LAImax': ([0.6], [0.4], [0.8])
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
