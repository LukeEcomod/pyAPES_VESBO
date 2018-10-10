# -*- coding: utf-8 -*-
"""
CANOPY MODEL PARAMETERS
"""
from copy import deepcopy
from parameters.utilities import lad_profiles
from forestfloor import forestfloor

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
           'WMA': False, # well-mixed assumption
           'StomaModel': 'MEDLYN_FARQUHAR',  # stomatal model
           'Ebal': True,  # computes leaf temperature by solving energy balance
           'SwModel': 'ZhaoQualls',
           'LwModel': 'ZhaoQualls',  #'Flerchinger'},  #
           'WaterStress': False,  # TRUE NOT SUPPORTED YET!
           'seasonal_LAI': True,  # account for seasonal LAI dynamics
           'pheno_cycle': True  # account for phenological cycle
           }

    # --- site location ---
    loc = {'lat': 60.63,  # latitude
           'lon': 23.95  # longitude
           }

    # --- aerodynamics ---
    aero = {'zos': 0.01,  # forest floor roughness length [m]
            'dPdx': 0.01,  # horizontal pressure gradient
            'Cd': 0.15,  # drag coefficient
            'Utop': 5.0,  # ensemble U/ustar
            'Ubot': 0.0,  # lower boundary
            'Sc': {'T': 2.0, 'H2O': 2.0, 'CO2': 2.0}  # Schmidt numbers
            }

    # --- radiation ---
    radi = {'clump': 0.7,  # clumping index [-]
            'leaf_angle': 1.0,  # leaf-angle distribution [-]
            'Par_alb': 0.12,  # shoot Par-albedo [-]
            'soil_Par_alb': 0.05,  # soil (moss) Par-albedo [-]
            'Nir_alb': 0.55,  # shoot NIR-albedo [-]
            'leaf_emi': 0.98
            }

    # --- interception and snowmodel ---  SADANNAN KORJAUSKERTOIMET?
    interc_snow = {'wmax': 0.15e-03, #0.5e-03,  # maximum interception storage capacity for rain [m per unit of LAI]
                   'wmaxsnow': 1.2e-03, #4.0e-03,  # maximum interception storage capacity for snow [m per unit of LAI]
                   'kmelt': 2.31e-08,  # Melting coefficient [m degC-1 s-1] (=2.0 mm/C/d)
                   'kfreeze': 5.79e-09,  # Freezing  coefficient [m degC-1 s-1] (=0.5 mm/C/d)
                   'retention': 0.05,  # max fraction of liquid water in snow [-]
                   'Tmin': 0.0,  # temperature below which all is snow [degC]
                   'Tmax': 1.0,  # temperature above which all is water [degC]
                   'Tmelt': 0.0,  # temperature when melting starts [degC]
                   'w_ini': 0.0,  # initial canopy storage [m]
                   'swe_ini': 0.0,  # initial snow water equivalent [m]
                   }

    # --- forest floor ---
    # defined in parameters.forestfloor.py
    ffloor = forestfloor

    # --- default values for plant characteristics ---
    gamma = 1.5  # adjust shoot light response
    plant_default = {'phenop': {  # --- seasonal cycle of phenology ---
                                'Xo': 0.0,  # initial delayed temperature [degC]
                                'fmin': 0.01,  # minimum photocapacity [-]
                                'Tbase': -4.67,  # base temperature [degC]
                                'tau': 8.33,  # time constant [days]
                                'smax': 18.5  #threshold for full acclimation [degC]
                                },
                    'laip': {  # --- leaf-area seasonal dynamics ---
                             'lai_min': 0.1,  # minimum LAI, fraction of annual maximum [-]
                             'lai_ini': None,  # initial LAI fraction, if None lai_ini = Lai_min * LAImax
                             'DDsum0': 0.0,  # degreedays at initial time [days]
                             'Tbase': 5.0,  # base temperature [degC]
                             'ddo': 45.0,  # degreedays at bud burst [days]
                             'ddur': 23.0,  # duration of recovery period [days]
                             'sso': 240.0,  # start doy of decrease, based on daylength
                             'sdur': 30.0, # duration of decreasing period [days]
                             },
                    'photop': {  # --- leaf gas-exchange parameters ---
                            'Vcmax': None,  # maximum carboxylation velocity [umolm-2s-1]
                            'Jmax': None,  # maximum rate of electron transport [umolm-2s-1]
                            'Rd': None,  # dark respiration rate [umolm-2s-1]
                            'alpha': gamma * 0.2,  # quantum yield parameter [mol/mol]
                            'theta': 0.7,  # co-limitation parameter of Farquhar-model
                            'La':None,  # stomatal parameter (Lambda, m, ...) depending on model
                            'm':None,
                            'g0': 1.0e-3,  # residual conductance for CO2 [molm-2s-1]
                            'kn': 0.6,
                            'beta': 0.95,  # co-limitation parameter of Farquhar-model
                            'drp': 0.7,
                            'tresp': {  # --- temperature sensitivity parameters ---
                                    'Vcmax': [],  # [Ha, Hd, Topt]; activation energy [kJmol-1], deactivation energy [kJmol-1], optimum temperature [degC]
                                    'Jmax': [],  # [Ha, Hd, Topt];
                                    'Rd': [33.0]}  #[Ha]; activation energy [kJmol-1)]
                            },
                    'leafp': {  # --- leaf properties ---
                            'lt': 0.02,  # leaf lengthscale [m]
                            'par_alb': 0.12,  # leaf Par albedo [-]
                            'nir_alb': 0.55,  # leaf Nir albedo [-]
                            'emi': 0.98  # leaf emissivity [-]
                            },
                    'rootp': {  # --- root zone properties ----
                            'root_depth': 0.2,  # root depth [m]
                            'beta': 0.943,  # shape parameter for root distribution model. Y=1-beta.^z_in_cm; Y=cumulative distribution (Gale & Grigal 1987)
                            'RAI_LAI_multiplier': 2.0,  # multiplier for total fine root area index (RAI = 2*LAImax)
                            'fine_radius': 2e-3,  # fine root radius [m]
                            'radial_K': 5e-8  # maximum bulk root membrane conductance in radial direction [s-1]
                            },
                    }
    pine = deepcopy(plant_default)  # initialize with default
    spruce = deepcopy(plant_default)
    decid = deepcopy(plant_default)
    shrubs = deepcopy(plant_default)

    # --- stand characteristics ---
    # specify name and maximum leaf-area index, LAImax [m2/m2], and
    # adjust values of 'phenop' and 'laip' if default values not suitable
    pine.update({'name': 'pine', 'LAImax': []}) #[2.1]})
    pine['phenop'].update({'fmin': 0.1})
    pine['laip'].update({'lai_min': 0.8})

    spruce.update({'name': 'spruce', 'LAImax': []}) #[1.0]})
    spruce['phenop'].update({'fmin': 0.1})
    spruce['laip'].update({'lai_min': 0.8})

    decid.update({'name': 'decid', 'LAImax': []}) #[1.0]})

    shrubs.update({'name': 'shrubs', 'LAImax': [0.7]})
    shrubs['laip'].update({'lai_min': 0.5})

    # normed leaf area density profiles
    quantiles = [1.0]  # quantiles used in creating species stand lad profiles
    hs = 0.5  # height of understory shrubs [m]
    pine['lad'], spruce['lad'], decid['lad'], shrubs['lad'], lai_p, lai_s, lai_d = lad_profiles(
            grid, dbhfile, quantiles, hs, plot=False)
    if pine['LAImax'] == []:
        pine['LAImax'] = lai_p
    if spruce['LAImax'] == []:
        spruce['LAImax'] = lai_s
    if decid['LAImax'] == []:
        decid['LAImax'] = lai_d

    # adjust leaf gas-exchange parameters
    gfact = 1.2  # coefficient for adjusting (?)
    pine['photop'].update({'Vcmax': 55.0, 'Jmax': 105.0, 'Rd': 1.3,
                           'La': 1600.0, 'm': gfact*2.5})
    pine['photop']['tresp'].update({'Vcmax': [78, 200.0, 650.0],
                                    'Jmax': [56, 200.0, 647.0]})
    spruce['photop'].update({'Vcmax': 60.0, 'Jmax': 114.0, 'Rd': 1.5,
                             'La': 1600.0, 'm': gfact*2.5})
    spruce['photop']['tresp'].update({'Vcmax': [53.2, 202.0, 640.3],  # Tarvainen et al. 2013 Oecologia
                                      'Jmax': [38.4, 202.0, 655.8]})
    decid['photop'].update({'Vcmax': 50.0, 'Jmax': 95.0, 'Rd': 1.3,
                            'La': 600.0, 'm': gfact*4.5})
    decid['photop']['tresp'].update({'Vcmax': [77.0, 200.0, 636.7],  # Medlyn et al 2002.
                                     'Jmax': [42.8, 200.0, 637.0]})
    decid['leafp'].update({'lt': 0.05})
    shrubs['photop'].update({'Vcmax': 50.0, 'Jmax': 95.0, 'Rd': 1.3,
                            'La': 600.0, 'm': gfact*4.5, 'kn': 0.3})
    shrubs['photop']['tresp'].update({'Vcmax': [77.0, 200.0, 636.7],
                                     'Jmax': [42.8, 200.0, 637.0]})
    plant_types = [pine, spruce, decid, shrubs]

    cpara.update({'ctr': ctr, 'loc': loc, 'radi': radi, 'aero': aero,'grid': grid,
                  'interc_snow': interc_snow, 'plant_types': plant_types, 'ffloor': ffloor})

    return cpara