# -*- coding: utf-8 -*-
"""
Planttype parameters
"""

import numpy as np

"""
planttypes (list):
    i. planttype_i (dict):
        'name' (str): name of planttype
        'LAImax' (list): leaf area index of planttype groups
        'lad' (array of lists): normalized leaf area density profiles of planttype groups
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
        'rootp' (dict): root zone properties
            'root_depth': root depth [m]
            'beta': shape parameter for root distribution model
            'RAI_LAI_multiplier': multiplier for total fine root area index (RAI = 2*LAImax)
            'fine_radius': fine root radius [m]
            'radial_K': maximum bulk root membrane conductance in radial direction [s-1]
"""

# normalized leaf area density
# note! len(lad_xx) == number of canopy layers
# and len(lad_xx[i]) == len(LAImax) == number of groups in planttype
lad_pine = np.array([[0.00E+00], [0.00E+00], [7.53E-06], [2.48E-05], [4.44E-05], [6.12E-05], [7.44E-05], [8.23E-05], [8.67E-05], [9.15E-05], [9.65E-05], [1.07E-04], [1.35E-04], [1.63E-04], [2.25E-04], [3.33E-04], [4.70E-04], [6.39E-04], [8.32E-04], [1.08E-03], [1.38E-03], [1.76E-03], [2.20E-03], [2.73E-03], [3.33E-03], [4.01E-03], [4.81E-03], [5.96E-03], [7.52E-03], [9.55E-03], [1.22E-02], [1.59E-02], [2.05E-02], [2.60E-02], [3.31E-02], [4.20E-02], [5.27E-02], [6.63E-02], [8.22E-02], [9.99E-02], [1.19E-01], [1.38E-01], [1.55E-01], [1.70E-01], [1.80E-01], [1.87E-01], [1.88E-01], [1.87E-01], [1.81E-01], [1.73E-01], [1.63E-01], [1.50E-01], [1.38E-01], [1.26E-01], [1.10E-01], [9.73E-02], [8.28E-02], [6.90E-02], [5.69E-02], [4.37E-02], [3.17E-02], [2.22E-02], [1.50E-02], [9.22E-03], [4.41E-03], [2.62E-03], [1.40E-03], [5.44E-04], [2.01E-04], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00]])
lad_spruce = np.array([[0.00E+00], [0.00E+00], [2.52E-02], [4.92E-02], [6.36E-02], [6.94E-02], [6.51E-02], [6.40E-02], [5.48E-02], [5.30E-02], [5.41E-02], [4.58E-02], [4.60E-02], [4.76E-02], [4.47E-02], [4.12E-02], [4.28E-02], [4.40E-02], [3.93E-02], [3.96E-02], [4.10E-02], [4.12E-02], [3.70E-02], [3.85E-02], [3.98E-02], [3.79E-02], [3.64E-02], [3.83E-02], [4.03E-02], [3.90E-02], [4.15E-02], [4.43E-02], [4.51E-02], [4.65E-02], [4.96E-02], [5.18E-02], [5.25E-02], [5.55E-02], [5.77E-02], [5.78E-02], [6.04E-02], [6.17E-02], [6.14E-02], [6.30E-02], [6.32E-02], [6.31E-02], [6.39E-02], [6.23E-02], [6.24E-02], [6.13E-02], [5.97E-02], [5.92E-02], [5.70E-02], [5.63E-02], [5.23E-02], [5.11E-02], [4.83E-02], [4.69E-02], [4.18E-02], [4.00E-02], [3.54E-02], [3.39E-02], [3.14E-02], [2.94E-02], [2.74E-02], [2.59E-02], [2.40E-02], [2.17E-02], [1.86E-02], [1.68E-02], [1.52E-02], [1.39E-02], [1.34E-02], [1.13E-02], [9.06E-03], [7.36E-03], [6.60E-03], [6.36E-03], [6.11E-03], [5.87E-03], [5.63E-03], [5.39E-03], [5.16E-03], [4.23E-03], [3.20E-03], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00]])
lad_decid = np.array([[0.00E+00], [0.00E+00], [6.88E-03], [1.31E-02], [1.62E-02], [2.09E-02], [2.88E-02], [3.54E-02], [3.86E-02], [4.24E-02], [4.63E-02], [4.86E-02], [5.10E-02], [5.17E-02], [5.43E-02], [5.71E-02], [5.90E-02], [6.03E-02], [6.12E-02], [6.19E-02], [5.89E-02], [5.99E-02], [6.18E-02], [6.39E-02], [6.57E-02], [6.74E-02], [6.76E-02], [6.79E-02], [7.12E-02], [7.48E-02], [7.87E-02], [8.28E-02], [8.38E-02], [8.78E-02], [9.15E-02], [9.40E-02], [9.48E-02], [9.32E-02], [9.26E-02], [9.10E-02], [8.81E-02], [8.34E-02], [7.99E-02], [7.59E-02], [7.02E-02], [6.54E-02], [6.10E-02], [5.54E-02], [5.08E-02], [4.68E-02], [4.15E-02], [3.81E-02], [3.43E-02], [3.12E-02], [2.75E-02], [2.47E-02], [2.12E-02], [1.85E-02], [1.54E-02], [1.36E-02], [1.15E-02], [9.40E-03], [7.71E-03], [6.31E-03], [5.26E-03], [3.69E-03], [2.34E-03], [1.84E-03], [1.46E-03], [1.13E-03], [1.03E-03], [9.35E-04], [8.51E-04], [7.73E-04], [8.26E-05], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00]])
lad_shrubs = np.array([[0.00E+00], [3.30E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00], [0.00E+00]])

Pine = {
        'name': 'pine',
        'LAImax': [2.1],
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
            'La': 1600.0,
            'm': 2.5,
            'g0': 4.0e-3,
            'kn': 0.5,
            'beta': 0.95,
            'drp': 0.7
            },
        'leafp': {
            'lt': 0.02,
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
        'LAImax': [1.0],
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
            'La': 1600.0,
            'm': 2.5,
            'g0': 4.0e-3,
            'kn': 0.5,
            'beta': 0.95,
            'drp': 0.7
            },
        'leafp': {
            'lt': 0.02,
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
        'LAImax': [1.2],
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
            'La': 600.0,
            'm': 4.5,
            'g0': 1.0e-2,
            'kn': 0.2,
            'beta': 0.95,
            'drp': 0.7
            },
        'leafp': {
            'lt': 0.05,
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
        'LAImax': [0.7],
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
            'La': 600.0,
            'm': 4.5,
            'g0': 1.0e-2,
            'kn': 0.0,
            'beta': 0.95,
            'drp': 0.7
            },
        'leafp': {
            'lt': 0.05,
            },
        'rootp': {
            'root_depth': 0.2,
            'beta': 0.943,
            'RAI_LAI_multiplier': 2.0,
            'fine_radius': 2.0e-3,
            'radial_K': 5.0e-8,
            }
        }

planttypes = {'pine': Pine,
              'spruce': Spruce,
              'decidious': Decidious,
              'shrubs': Shrubs}