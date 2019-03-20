# -*- coding: utf-8 -*-
"""
Planttype parameters - default parameterization
"""

import numpy as np

"""
planttypes (list):
    i. planttype_i (dict):
        'name' (str): name of planttype
        'LAImax' (list): leaf area index of planttype groups
        'lad' (array of lists): normalized leaf area density profiles of planttype groups
        'phenop' (dict): parameters for seasonal cycle of photosynthetic activity
            'Xo': initial delayed temperature [degC]
            'fmin': minimum photocapacity [-]
            'Tbase': base temperature [degC]
            'tau': time constant [days]
            'smax': threshold temperature for full acclimation [degC]
        'laip' (dict): parameters for LAI seasonal dynamics
            'lai_min': minimum LAI, fraction of annual maximum [-]
            'lai_ini': initial LAI fraction, if None lai_ini = Lai_min * LAImax
            'DDsum0': degreedays at initial time [days]
            'Tbase': base temperature [degC]
            'ddo': degreedays at bud burst [days]
            'sdl':  daylength for senescence start [h]
            'sdur': duration of decreasing period [days]
        'photop' (dict): leaf gas-exchange parameters
            'Vcmax': maximum carboxylation velocity [umolm-2s-1]
            'Jmax': maximum rate of electron transport [umolm-2s-1]
            'Rd': dark respiration rate [umolm-2s-1]
            'alpha': quantum yield parameter [mol/mol]
            'theta': co-limitation parameter of Farquhar-model
            'La': stomatal parameter (Lambda, g1, ...) depending on model
            'g1':
            'g0': residual conductance for CO2 [molm-2s-1]
            'kn': used to scale photosynthetic capacity (vertical N gradient)
            'beta':  co-limitation parameter of Farquhar-model
            'drp'(dict): drought response parameters
            'tresp' (dict): temperature sensitivity parameters
                'Vcmax': [Ha, Hd, dS]; activation energy [kJmol-1], deactivation energy [kJmol-1],  entropy factor [J mol-1]
                'Jmax': [Ha, Hd, dS];
                'Rd': [Ha]; activation energy [kJmol-1)]
        'leafp' (dict): leaf properties
            'lt': leaf lengthscale [m]
        'rootp' (dict): root zone properties
            'root_depth': root depth [m]
            'beta': shape parameter for root distribution model
            'RAI_LAI_multiplier': multiplier for total fine root area index (RAI = 2*LAImax)
            'fine_radius': fine root radius [m]
            'radial_K': maximum bulk root membrane conductance in radial direction [s-1]
"""

# normalized leaf area density
# note! len(lad_xx) == number of canopy layers
lad_pine = np.array([0., 0., 7.53e-06, 2.48e-05, 4.44e-05, 6.12e-05, 7.44e-05, 8.23e-05, 8.67e-05, 9.15e-05, 9.65e-05, 1.07e-04, 1.35e-04, 1.63e-04, 2.25e-04, 3.33e-04, 4.70e-04, 6.39e-04, 8.32e-04, 1.08e-03, 1.38e-03, 1.76e-03, 2.20e-03, 2.73e-03, 3.33e-03, 4.01e-03, 4.81e-03, 5.96e-03, 7.52e-03, 9.55e-03, 1.22e-02, 1.59e-02, 2.05e-02, 2.60e-02, 3.31e-02, 4.20e-02, 5.27e-02, 6.63e-02, 8.22e-02, 9.99e-02, 1.19e-01, 1.38e-01, 1.55e-01, 1.70e-01, 1.80e-01, 1.87e-01, 1.88e-01, 1.87e-01, 1.81e-01, 1.73e-01, 1.63e-01, 1.50e-01, 1.38e-01, 1.26e-01, 1.10e-01, 9.73e-02, 8.28e-02, 6.90e-02, 5.69e-02, 4.37e-02, 3.17e-02, 2.22e-02, 1.50e-02, 9.22e-03, 4.41e-03, 2.62e-03, 1.40e-03, 5.44e-04, 2.01e-04, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0])
lad_spruce = np.array([0., 0., 0.0252, 0.0492, 0.0636, 0.0694, 0.0651, 0.064, 0.0548, 0.053, 0.0541, 0.0458, 0.046, 0.0476, 0.0447, 0.0412, 0.0428, 0.044, 0.0393, 0.0396, 0.041, 0.0412, 0.037, 0.0385, 0.0398, 0.0379, 0.0364, 0.0383, 0.0403, 0.039, 0.0415, 0.0443, 0.0451, 0.0465, 0.0496, 0.0518, 0.0525, 0.0555, 0.0577, 0.0578, 0.0604, 0.0617, 0.0614, 0.063, 0.0632, 0.0631, 0.0639, 0.0623, 0.0624, 0.0613, 0.0597, 0.0592, 0.057, 0.0563, 0.0523, 0.0511, 0.0483, 0.0469, 0.0418, 0.04, 0.0354, 0.0339, 0.0314, 0.0294, 0.0274, 0.0259, 0.024, 0.0217, 0.0186, 0.0168, 0.0152, 0.0139, 0.0134, 0.0113, 0.00906, 0.00736, 0.0066, 0.00636, 0.00611, 0.00587, 0.00563, 0.00539, 0.00516, 0.00423, 0.0032, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
lad_decid = np.array([0., 0., 6.88e-03, 1.31e-02, 1.62e-02, 2.09e-02, 2.88e-02, 3.54e-02, 3.86e-02, 4.24e-02, 4.63e-02, 4.86e-02, 5.10e-02, 5.17e-02, 5.43e-02, 5.71e-02, 5.90e-02, 6.03e-02, 6.12e-02, 6.19e-02, 5.89e-02, 5.99e-02, 6.18e-02, 6.39e-02, 6.57e-02, 6.74e-02, 6.76e-02, 6.79e-02, 7.12e-02, 7.48e-02, 7.87e-02, 8.28e-02, 8.38e-02, 8.78e-02, 9.15e-02, 9.40e-02, 9.48e-02, 9.32e-02, 9.26e-02, 9.10e-02, 8.81e-02, 8.34e-02, 7.99e-02, 7.59e-02, 7.02e-02, 6.54e-02, 6.10e-02, 5.54e-02, 5.08e-02, 4.68e-02, 4.15e-02, 3.81e-02, 3.43e-02, 3.12e-02, 2.75e-02, 2.47e-02, 2.12e-02, 1.85e-02, 1.54e-02, 1.36e-02, 1.15e-02, 9.40e-03, 7.71e-03, 6.31e-03, 5.26e-03, 3.69e-03, 2.34e-03, 1.84e-03, 1.46e-03, 1.13e-03, 1.03e-03, 9.35e-04, 8.51e-04, 7.73e-04, 8.26e-05, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])
lad_shrubs = np.array([0., 3.3, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

Pine = {
        'name': 'pine',
        'LAImax': 2.1,
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
            'sdl': 12.0,
            'sdur': 30.0
            },
        'photop': {
            'Vcmax': 50.,
            'Jmax': 98.,  # 1.97*Vcmax (Kattge and Knorr, 2007)
            'Rd': 1.2,  # 0.023*Vcmax
            'tresp': {
                'Vcmax': [72., 200., 649.],  # (Kattge and Knorr, 2007)
                'Jmax': [50., 200., 646.],  # (Kattge and Knorr, 2007)
                'Rd': [33.0]
                },
            'alpha': 0.2,
            'theta': 0.7,
            'g1': 2.5,
            'g0': 4.0e-3,
            'kn': 0.5,
            'beta': 0.95,
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.02,
            },
        'rootp': {
            'root_depth': 0.5,
            'beta': 0.943,
            'RAI_LAI_multiplier': 2.0,
            'fine_radius': 2.0e-3,
            'root_cond': 5e8,  # [s]
            }
        }

Spruce = {
        'name': 'spruce',
        'LAImax': 1.0,
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
            'sdl': 12.0,
            'sdur': 30.0
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
            'alpha': 0.2,
            'theta': 0.7,
            'g1': 2.5,
            'g0': 4.0e-3,
            'kn': 0.5,
            'beta': 0.95,
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.02,
            },
        'rootp': {
            'root_depth': 0.5,
            'beta': 0.943,
            'RAI_LAI_multiplier': 2.0,
            'fine_radius': 2.0e-3,
            'root_cond': 5e8,  # [s]
            }
        }

Decidious = {
        'name': 'decidious',
        'LAImax': 1.2,
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
            'sdl': 12.0,
            'sdur': 30.0
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
            'alpha': 0.2,
            'theta': 0.7,
            'g1': 4.5,
            'g0': 1.0e-2,
            'kn': 0.5,
            'beta': 0.95,
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.05,
            },
        'rootp': {
            'root_depth': 0.5,
            'beta': 0.943,
            'RAI_LAI_multiplier': 2.0,
            'fine_radius': 2.0e-3,
            'root_cond': 5e8,  # [s]
            }
        }

Shrubs = {
        'name': 'shrubs',
        'LAImax': 0.7,
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
            'sdl': 12.0,
            'sdur': 30.0
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
            'alpha': 0.2,
            'theta': 0.7,
            'g1': 4.5,
            'g0': 1.0e-2,
            'kn': 0.0,
            'beta': 0.95,
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.05,
            },
        'rootp': {
            'root_depth': 0.3,
            'beta': 0.943,
            'RAI_LAI_multiplier': 2.0,
            'fine_radius': 2.0e-3,
            'root_cond': 5e8,  # [s]
            }
        }

planttypes = {'pine': Pine,
              'spruce': Spruce,
              'decidious': Decidious,
              'shrubs': Shrubs
              }
