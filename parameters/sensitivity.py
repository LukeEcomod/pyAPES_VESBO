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
                                                'ground_coverage': (0.0),
#                                                'height': (0.06)
                                                },
                                        'sphagnum': {
                                                'ground_coverage': (0.0),
#                                                'height': (0.045)
                                                },
                                        'pleurozium': {
                                                'ground_coverage': (0.0),
#                                                'height': (0.095)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.0)
                                        },
                                'baresoil': {
                                        'ground_coverage': (1.0)
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
                                                'ground_coverage': (0.07),
                                                'height': (0.095)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.71)
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.22)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.49])
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
                                                'ground_coverage': (0.5),
                                                'height': (0.095)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.5)
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.0)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.8])
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
                                                'ground_coverage': (0.4),
                                                'height': (0.095)
                                                }
                                        },
                                'litter': {
                                        'ground_coverage': (0.6)
                                        },
                                'baresoil': {
                                        'ground_coverage': (0.0)
                                        }
                                },
                        'planttypes': {
                                'shrubs': {
                                        'LAImax': ([0.65])
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