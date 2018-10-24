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
                'count': 2,
                'canopy': {
                        'forestfloor': {
                                'bryophytes': {
                                        'hylocomium': {
                                                'ground_coverage': (1.0, 1.0),
                                                'height': (0.05, 0.05)
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (0.0, 0.0),
                                                'height': (0.04, 0.04)
                                                }
                                                },
                                'baresoil': {
                                        'ground_coverage': (0.0, 0.0)
                                        }
                                        },
                        'planttypes': {
                                'pine': {
                                        'LAImax': ([2.1], [2.1])
                                        },
                                'shrubs': {
                                        'LAImax': ([0.7], [0.0])
                                        }
                                }
                }
            }

    elif name.upper() == 'CLEARCUT':
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
                                                'ground_coverage': (0.0),
                                                'height': (0.095)
                                                }
                                        },
                                'baresoil': {
                                        'ground_coverage': (1.0)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.7])
                                        }
                                }
                }
            }

    elif name.upper() == 'PARTIAL':
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
                                                'ground_coverage': (0.3),
                                                'height': (0.095)
                                                }
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.7)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.7])
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
                                                'ground_coverage': (0.6),
                                                'height': (0.095)
                                                }
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.4)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.7])
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

    for key, value in parameters.iteritems():
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