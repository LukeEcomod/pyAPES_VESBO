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
                                                'ground_coverage': (0.3),
#                                                'height': (0.06)
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.3),
#                                                'height': (0.045)
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (0.0),
#                                                'height': (0.095)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.2)
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.2)
                                        }
                                },
                }
            }

    elif name.upper() == 'CLEARCUT':
        parameters = {
                'count': 1,
                'canopy': {
                        'forestfloor': {
                                'bryophytes': {
                                        'hylocomium': {
                                                'ground_coverage': (0.0, 0.0)
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.0, 0.0)
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (0.07, 0.07)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.71, 0.71)
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.22, 0.22)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.10], [0.49])
                                        }
                                }
                }
            }

    elif name.upper() == 'PARTIAL':
        parameters = {
                'count': 2,
                'canopy': {
                        'forestfloor': {
                                'bryophytes': {
                                        'hylocomium': {
                                                'ground_coverage': (0.0, 0.0)
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.0, 0.0)
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (0.49, 0.49)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.51, 0.51)
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.0, 0.0)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.66], [0.79])
                                        }
                                }
                }
            }

    elif name.upper() == 'CONTROL':
        parameters = {
                'count': 1,
                'canopy': {
                        'forestfloor': {
                                'bryophytes': {
                                        'hylocomium': {
                                                'ground_coverage': (0.0),
                                                'height': (0.06)
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.0),
                                                'height': (0.045)
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (0.49),
                                                'height': (0.095)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.51)
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.0)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.66])
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
