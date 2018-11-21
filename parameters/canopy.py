# -*- coding: utf-8 -*-
"""
CANOPY MODEL PARAMETERS
"""

from .planttype import get_planttypes
from .forestfloor import forestfloor

def get_cpara(dbhfile):

    dbhfile = "parameters/runkolukusarjat/" + dbhfile

    # initialize dictionary to store parameters
    cpara = {}

    # grid
    grid = {'zmax': 30.0,  # heigth of grid from ground surface [m]
            'Nlayers': 100  # number of layers in grid [-]
            }

    # --- control flags (True/False) ---
    ctr = {'Eflow': True,  # ensemble flow
           'WMA': True, # well-mixed assumption
           'StomaModel': 'MEDLYN_FARQUHAR',  # stomatal model
           'Ebal': False,  # computes leaf temperature by solving energy balance
           'SwModel': 'ZhaoQualls',
           'LwModel': 'ZhaoQualls',  #'Flerchinger'},  #
           'WaterStress': False,  # TRUE NOT SUPPORTED YET!
           'seasonal_LAI': True,  # account for seasonal LAI dynamics
           'pheno_cycle': True  # account for phenological cycle
           }

    # --- micrometeo ---
    micromet = {'zos': 0.01,  # forest floor roughness length [m]
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
    interception = {'wmax': 0.15e-03, #0.5e-03,  # maximum interception storage capacity for rain [m per unit of LAI]
                    'wmaxsnow': 1.2e-03, #4.0e-03,  # maximum interception storage capacity for snow [m per unit of LAI]
                    'w_ini': 0.0,  # initial canopy storage [m]
                    'Tmin': 0.0,  # temperature below which all is snow [degC]
                    'Tmax': 1.0,  # temperature above which all is water [degC]
                    }

    # --- forest floor ---
    # defined in parameters.forestfloor.py
    ffloor = forestfloor

    # --- plant types ---
    # defined in parameters.planttype.py
    planttypes = get_planttypes(dbhfile, grid)

    cpara.update({'ctr': ctr, 'grid': grid, 'radiation': radiation, 'micromet': micromet,
                  'interception': interception, 'planttypes': planttypes, 'forestfloor': ffloor})

    return cpara