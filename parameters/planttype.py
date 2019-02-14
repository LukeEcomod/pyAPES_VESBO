# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 10:23:58 2018

@author: Kersti Haahti
"""
from parameters.utilities import lad_profiles

def get_planttypes(dbhfile, grid):
    """
    Returns: 
        planttypes (list):
            i. planttype_i (dict):
                'name' (str): name of planttype
                'LAImax' (list): leaf area index of planttype groups
                'lad' (list of arrays): normalized leaf area density profiles of planttype groups
                'phenop' (dict): parameters for seasonal cycle of phenology
                    'Xo': initial delayed temperature [degC]
                    'fmin': minimum photocapacity [-]
                    'Tbase': base temperature [degC]
                    'tau': time constant [days]
                    'smax': threshold for full acclimation [degC]
                'laip' (dict): parameters forleaf-area seasonal dynamics
                    'lai_min': minimum LAI, fraction of annual maximum [-]
                    'lai_ini': initial LAI fraction, if None lai_ini = Lai_min * LAImax
                    'DDsum0': degreedays at initial time [days]
                    'Tbase': base temperature [degC]
                    'ddo': degreedays at bud burst [days]
                    'ddur': duration of recovery period [days]
                    'sdl':  daylength for senescence start [h]
                    'sdur': duration of decreasing period [days]
                'photop' (dict): leaf gas-exchange parameters
                    'Vcmax': maximum carboxylation velocity [umolm-2s-1]
                    'Jmax': maximum rate of electron transport [umolm-2s-1]
                    'Rd': dark respiration rate [umolm-2s-1]
                    'alpha': quantum yield parameter [mol/mol]
                    'theta': co-limitation parameter of Farquhar-model
                    'La': stomatal parameter (Lambda, m, ...) depending on model
                    'm':
                    'g0': residual conductance for CO2 [molm-2s-1]
                    'kn': used to scale photosynthetic capacity (vertical N gradient)
                    'beta':  co-limitation parameter of Farquhar-model
                    'drp':
                    'tresp' (dict): temperature sensitivity parameters
                        'Vcmax': [Ha, Hd, dS]; activation energy [kJmol-1], deactivation energy [kJmol-1],  entropy factor [J mol-1]
                        'Jmax': [Ha, Hd, dS];
                        'Rd': [Ha]; activation energy [kJmol-1)]
                'leafp' (dict): leaf properties
                    'lt': leaf lengthscale [m]
                    'par_alb': leaf Par albedo [-]
                    'nir_alb': leaf Nir albedo [-]
                    'emi': leaf emissivity [-]
                'rootp' (dict): root zone properties
                    'root_depth': root depth [m]
                    'beta': shape parameter for root distribution model
                    'RAI_LAI_multiplier': multiplier for total fine root area index (RAI = 2*LAImax)
                    'fine_radius': fine root radius [m]
                    'radial_K': maximum bulk root membrane conductance in radial direction [s-1]
    """

    # normed leaf area density profiles
    quantiles = [1.0]  # quantiles used in creating tree lad profiles (groups of same planttype)

    hs = 0.5  # height of understory shrubs [m]
    lai_shrubs = [0.7]

    lad_pine, lad_spruce, lad_decid, lad_shrubs, lai_pine, lai_spruce, lai_decid = lad_profiles(
            grid, dbhfile, quantiles, hs, plot=False)

#    #define LAI other than derived from dbhfile (must be of len(quantiles))
#    lai_pine = []
#    lai_spruce = []
#    lai_decid = []

    gamma = 1.0  # adjust shoot light response
    gfact = 1.0  # coefficient for adjusting leaf gas-exchange parameters

    Pine = {
            'name': 'pine',
            'LAImax': lai_pine,
            'lad': lad_pine,
            'phenop': {
                'Xo': 0.0,
                'fmin': 0.1,
                'Tbase': -4.67,  # Kolari 2007
                'tau': 8.33,  # Kolari 2007
                'smax': 15.0  # Kolari 2014
                },
            'laip': {
                'lai_min': 0.8,
                'lai_ini': None,
                'DDsum0': 0.0,
                'Tbase': 5.0,
                'ddo': 45.0,
                'ddmat': 250.0,
                'ddur': 23.0,
                'sdl': 12.0,
                'sdur': 30.0
                },
            'photop': {
                'Vcmax': 50., #94.0,  # Tarvainen et al. 2018 Physiol. Plant.  # SAMULILLA 50.0
                'Jmax': 98., #143.0,  # SAMULILLA  90.0
                'Rd': 1.2, #1.3,  # SAMULILLA 1.2
                'tresp': {
                    'Vcmax': [72., 200., 649.],  #[78.3, 200.0, 650.1],
                    'Jmax': [50., 200., 646.], #[56.0, 200.0, 647.9],
                    'Rd': [33.0]
                    },
                'alpha': gamma * 0.2,
                'theta': 0.7,
                'La': 1600.0,
                'm': gfact * 2.5, # --- check value!! # SAMULILLA 2.0-2.3
                'g0': 4.0e-3,  # SAMULILLA 0.004
                'kn': 0.5,
                'beta': 0.95,
                'drp': 0.7
                },
            'leafp': {
                'lt': 0.02,  # SAMULILLA 0.05??
                'par_alb': 0.12,
                'nir_alb': 0.55,
                'emi': 0.98
                },
            'rootp': {
                'root_depth': 0.2,
                'beta': 0.943,
                'RAI_LAI_multiplier': 2.0,
                'fine_radius': 2.0e-3,
                'radial_K': 5.0e-8,
                }
            }

    Spruce = {
            'name': 'spruce',
            'LAImax': lai_spruce,
            'lad': lad_spruce,
            'phenop': {
                'Xo': 0.0,
                'fmin': 0.1,
                'Tbase': -4.67,  # Kolari 2007
                'tau': 8.33,  # Kolari 2007
                'smax': 15.0  # Kolari 2014
                },
            'laip': {
                'lai_min': 0.8,
                'lai_ini': None,
                'DDsum0': 0.0,
                'Tbase': 5.0,
                'ddo': 45.0,
                'ddmat': 250.0,
                'ddur': 23.0,
                'sdl': 12.0,
                'sdur': 30.0
                },
            'photop': {
                'Vcmax': 60., #69.7,  # Tarvainen et al. 2013 Oecologia # SAMULILLA 60.0
                'Jmax': 118., #130.2,  # SAMULILLA 114
                'Rd': 1.4, #1.3,  # SAMULILLA 1.5
                'tresp': {
                    'Vcmax': [72., 200., 649.],  #[53.2, 200.0, 640.0],
                    'Jmax': [50., 200., 646.], #[38.4, 200.0, 655.5],
                    'Rd': [33.0]
                    },
                'alpha': gamma * 0.2,
                'theta': 0.7,
                'La': 1600.0,
                'm': gfact * 2.5, # SAMULILLA 2.0-2.3
                'g0': 4.0e-3,  # SAMULILLA 0.004
                'kn': 0.5,
                'beta': 0.95,
                'drp': 0.7
                },
            'leafp': {
                'lt': 0.02,  # SAMULILLA 0.05??
                'par_alb': 0.12,
                'nir_alb': 0.55,
                'emi': 0.98
                },
            'rootp': {
                'root_depth': 0.2,
                'beta': 0.943,
                'RAI_LAI_multiplier': 2.0,
                'fine_radius': 2.0e-3,
                'radial_K': 5.0e-8,
                }
            }

    Decidious = {
            'name': 'decidious',
            'LAImax': lai_decid,
            'lad': lad_decid,
            'phenop': {
                'Xo': 0.0,
                'fmin': 0.01,
                'Tbase': -4.67,  # Kolari 2007
                'tau': 8.33,  # Kolari 2007
                'smax': 15.0  # Kolari 2014
                },
            'laip': {
                'lai_min': 0.1,
                'lai_ini': None,
                'DDsum0': 0.0,
                'Tbase': 5.0,
                'ddo': 45.0,
                'ddmat': 250.0,
                'ddur': 23.0,
                'sdl': 12.0,
                'sdur': 30.0
                },
            'photop': {
                'Vcmax': 45., #69.1,  # Medlyn et al 2002. Plant, Cell and Env.   # SAMULILLA 45.0
                'Jmax': 89., #116.3,  # SAMULILLA 80.0
                'Rd': 1.0,
                'tresp': {
                    'Vcmax': [72., 200., 649.],  #[77.0, 200.0, 636.4],
                    'Jmax': [50., 200., 646.], #[42.8, 200.0, 636.6],
                    'Rd': [33.0]
                    },
                'alpha': gamma * 0.2,
                'theta': 0.7,
                'La': 600.0,
                'm': gfact * 4.5,  # SAMULILLA 4.2 – 4.3
                'g0': 1.0e-2,  # SAMULILLA 0.01
                'kn': 0.5,  # SAMULILLA 0.2
                'beta': 0.95,
                'drp': 0.7
                },
            'leafp': {
                'lt': 0.05,
                'par_alb': 0.12,
                'nir_alb': 0.55,
                'emi': 0.98
                },
            'rootp': {
                'root_depth': 0.2,
                'beta': 0.943,
                'RAI_LAI_multiplier': 2.0,
                'fine_radius': 2.0e-3,
                'radial_K': 5.0e-8,
                }
            }

    Shrubs = {
            'name': 'shrubs',
            'LAImax': lai_shrubs,
            'lad': lad_shrubs,
            'phenop': {
                'Xo': 0.0,
                'fmin': 0.01,
                'Tbase': -4.67,  # Kolari 2007
                'tau': 8.33,  # Kolari 2007
                'smax': 15.0  # Kolari 2014
                },
            'laip': {
                'lai_min': 0.5,
                'lai_ini': None,
                'DDsum0': 0.0,
                'Tbase': 5.0,
                'ddo': 45.0,
                'ddmat': 250.0,
                'ddur': 23.0,
                'sdl': 12.0,
                'sdur': 30.0
                },
            'photop': {
                'Vcmax': 45., #69.1,  # Medlyn et al 2002. Plant, Cell and Env.  # SAMULILLA 40.0
                'Jmax': 89., #116.3,  # SAMULILLA 76.0
                'Rd': 1.0, #1.3,  # SAMULILLA 0.7
                'tresp': {
                    'Vcmax': [72., 200., 649.],  #[77.0, 200.0, 636.4],
                    'Jmax': [50., 200., 646.], #[42.8, 200.0, 636.6],
                    'Rd': [33.0]
                    },
                'alpha': gamma * 0.2,
                'theta': 0.7,
                'La': 600.0,
                'm': gfact * 4.5,  # SAMULILLA ??
                'g0': 1.0e-2,  # SAMULILLA ??
                'kn': 0.0,
                'beta': 0.95,
                'drp': 0.7
                },
            'leafp': {
                'lt': 0.05,
                'par_alb': 0.12,
                'nir_alb': 0.55,
                'emi': 0.98
                },
            'rootp': {
                'root_depth': 0.2,
                'beta': 0.943,
                'RAI_LAI_multiplier': 2.0,
                'fine_radius': 2.0e-3,
                'radial_K': 5.0e-8,
                }
            }

    planttypes = {'pine': Pine, 'spruce': Spruce, 'decidious': Decidious, 'shrubs': Shrubs}

    return planttypes