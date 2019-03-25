# -*- coding: utf-8 -*-
"""
CANOPY MODEL PARAMETERS
"""

from .planttype import planttypes
from .forestfloor import forestfloor

# site location
loc = {'lat': 61.51,  # latitude
       'lon': 24.0  # longitude
       }

# grid
grid = {'zmax': 30.0,  # heigth of grid from ground surface [m]
        'Nlayers': 100  # number of layers in grid [-]
        }

# --- control flags (True/False) ---
ctr = {'Eflow': True,  # ensemble flow
       'WMA': False, # well-mixed assumption
       'Ebal': True,  # computes leaf temperature by solving energy balance
       'WaterStress': None,  #'PsiL',  # Rew or PsiL or None
       'seasonal_LAI': True,  # account for seasonal LAI dynamics
       'pheno_cycle': True  # account for phenological cycle
       }

# --- micrometeo ---
micromet = {'zos': 0.01,  # forest floor roughness length [m]  -- not used?
            'dPdx': 0.01,  # horizontal pressure gradient
            'Cd': 0.15,  # drag coefficient
            'Utop': 5.0,  # ensemble U/ustar
            'Ubot': 0.0,  # lower boundary
            'Sc': {'T': 2.0, 'H2O': 2.0, 'CO2': 2.0}  # Schmidt numbers
            }

# --- radiation ---
radiation = {'clump': 0.7,  # clumping index [-]
             'leaf_angle': 1.0,  # leaf-angle distribution [-]
             'Par_alb': 0.12,  # shoot Par-albedo [-]
             'Nir_alb': 0.55,  # shoot NIR-albedo [-]
             'leaf_emi': 0.98
             }

# --- interception ---  SADANNAN KORJAUSKERTOIMET?
interception = {'wmax': 0.2e-03,  # maximum interception storage capacity for rain [m per unit of LAI]  - Watanabe & Mizunani coniferous trees
                'wmaxsnow': 1.6e-03,  # maximum interception storage capacity for snow [m per unit of LAI]
                'w_ini': 0.0,  # initial canopy storage [m]
                'Tmin': 0.0,  # temperature below which all is snow [degC]
                'Tmax': 1.0,  # temperature above which all is water [degC]
                'leaf_orientation': 0.5 # leaf orientation factor for randomdly oriented leaves
                }

cpara = {'loc': loc,
         'ctr': ctr,
         'grid': grid,
         'radiation': radiation,
         'micromet': micromet,
         'interception': interception,
         'planttypes': planttypes,
         'forestfloor': forestfloor
         }
