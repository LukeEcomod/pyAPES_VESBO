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
                                },
#                        'planttypes': {
#                                'pine': {  # Medlyn et al. 2002
#                                        'photop': {
#                                                'Vcmax': (67.33),
#                                                'Jmax': (70.77),
#                                                'tresp': {
#                                                        'Vcmax': ([69.8, 200.0, 659.9]),
#                                                        'Jmax': ([100.3, 147.9, 511.0]),
#                                                        }
#                                                }
#                                        }
#                                'spruce': {  # Tarvainen et al. 2013
#                                        'photop': {
#                                                'Vcmax': (69.7),
#                                                'Jmax': (130.2),
#                                                'tresp': {
#                                                        'Vcmax': ([53.2, 200.0, 640.0]),
#                                                        'Jmax': ([38.4, 200.0, 655.5]),
#                                                        }
#                                                }
#                                        }
#                                'decidious': {  # Medlyn et al. 2002
#                                        'photop': {
#                                                'Vcmax': (101.9),
#                                                'Jmax': (111.89),
#                                                'tresp': {
#                                                        'Vcmax': ([63.8, 200.0, 655.0]),
#                                                        'Jmax': ([108.5, 156.8, 543.2]),
#                                                        }
#                                                }
#                                        }
#                                }
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
#    if name is None:
#        # parameter sets of simulations
#        parameters = {
#                'count': 2,
#                'canopy': {
#                        'forestfloor': {
#                                'bryophytes': {
#                                        'hylocomium': {
#                                                'ground_coverage': (1.0, 1.0),
#                                                'height': (0.05, 0.05)
#                                                },
#                                        'pleurozium': {
#                                                'ground_coverage': (0.0, 0.0),
#                                                'height': (0.04, 0.04)
#                                                }
#                                                },
#                                'baresoil': {
#                                        'ground_coverage': (0.0, 0.0)
#                                        }
#                                        },
#                        'planttypes': {
#                                'pine': {
#                                        'LAImax': ([2.1], [2.1])
#                                        },
#                                'shrubs': {
#                                        'LAImax': ([0.7], [0.0])
#                                        }
#                                }
#                }
#            }
#
#    elif name.upper() == 'CLEARCUT':
#        parameters = {
#                'count': 1,
#                'canopy': {
#                        'forestfloor': {
#                                'bryophytes': {
#                                        'hylocomium': {
#                                                'ground_coverage': (0.0),
#                                                'height': (0.06)
#                                                },
#                                        'sphagnum': {
#                                                'ground_coverage': (0.0),
#                                                'height': (0.045)
#                                                },
#                                        'pleurozium': {
#                                                'ground_coverage': (0.0),
#                                                'height': (0.095)
#                                                }
#                                        },
#                                'baresoil': {
#                                        'ground_coverage': (1.0)
#                                        }
#                                },
#                        'planttypes': {
#                                'shrubs': {
#                                        'LAImax': ([0.7])
#                                        }
#                                }
#                }
#            }
#
#    elif name.upper() == 'PARTIAL':
#        parameters = {
#                'count': 1,
#                'canopy': {
#                        'forestfloor': {
#                                'bryophytes': {
#                                        'hylocomium': {
#                                                'ground_coverage': (0.0),
#                                                'height': (0.06)
#                                                },
#                                        'sphagnum': {
#                                                'ground_coverage': (0.0),
#                                                'height': (0.045)
#                                                },
#                                        'pleurozium': {
#                                                'ground_coverage': (0.3),
#                                                'height': (0.095)
#                                                }
#                                        },
#                                'baresoil': {
#                                        'ground_coverage': (0.7)
#                                        }
#                                },
#                        'planttypes': {
#                                'shrubs': {
#                                        'LAImax': ([0.7])
#                                        }
#                                }
#                }
#            }
#
#    elif name.upper() == 'CONTROL':
#        parameters = {
#                'count': 1,
#                'canopy': {
#                        'forestfloor': {
#                                'bryophytes': {
#                                        'hylocomium': {
#                                                'ground_coverage': (0.0),
#                                                'height': (0.06)
#                                                },
#                                        'sphagnum': {
#                                                'ground_coverage': (0.0),
#                                                'height': (0.045)
#                                                },
#                                        'pleurozium': {
#                                                'ground_coverage': (0.6),
#                                                'height': (0.095)
#                                                }
#                                        },
#                                'baresoil': {
#                                        'ground_coverage': (0.4)
#                                        }
#                                },
#                        'planttypes': {
#                                'shrubs': {
#                                        'LAImax': ([0.7])
#                                        }
#                                }
#                }
#            }

#    else:
#        raise ValueError('Unknown sensitivity parametrization %s' % name.upper())
#
#    return parameters

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
