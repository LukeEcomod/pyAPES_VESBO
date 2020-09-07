# -*- coding: utf-8 -*-
"""
DEFINES PARAMETERS FOR pyAPES SIMULATION at HYYTIALA SMEAR II -SITE

to import: from parameters.SmearII import gpara, cpara, spara

pyAPES output variables and logger config in parameters.outputs

Todo:
- currently testing that all works, all parameters are not properly selected
- implement reading canopy structure from tree inventory data
- for thinning-simulations and long-term trend paper, use parameters.parametersets to
    override nominal parameters defined here!

"""

import numpy as np
from tools.utilities import lad_weibul, lad_constant


gpara = {'dt' : 1800.0,  # timestep in forcing data file [s]
           'start_time' : "2005-06-01",  # start time of simulation [yyyy-mm-dd]
           'end_time' : "2005-06-10",  # end time of simulation [yyyy-mm-dd]
           'forc_filename' : "Hyytiala/FIHy_forcing_2005-2010.dat",  # forcing data file*
           'results_directory':'results/Hyytiala/'
         }

# --- control flags (True/False) ---
ctr = {'Eflow': True,  # ensemble flow
       'WMA': False,  #True,  #  well-mixed assumption
       'Ebal': True,  #False,  #  computes leaf temperature by solving energy balance
       'WaterStress': 'Rew',  #'PsiL',  # Rew or PsiL or None
       'seasonal_LAI': True,  # account for seasonal LAI dynamics
       'pheno_cycle': True  # account for phenological cycle
       }

# site location
loc = {'lat': 61.51,  # latitude
       'lon': 24.0  # longitude
       }

# grid
grid = {'zmax': 25.0,  # heigth of grid from ground surface [m]
        'Nlayers': 101  # number of layers in grid [-]
        }


# --- micrometeo ---
micromet = {'zos': 0.01,  # forest floor roughness length [m]  -- not used?
            'dPdx': 0.0,  # horizontal pressure gradient
            'Cd': 0.15,  # drag coefficient
            'Utop': 5.0,  # ensemble U/ustar
            'Ubot': 0.01,  # lower boundary
            'Sc': {'T': 2.0, 'H2O': 2.0, 'CO2': 2.0}  # Schmidt numbers
            }

# --- radiation ---
radiation = {'clump': 0.7,  # clumping index [-]
             'leaf_angle': 1.0,  # leaf-angle distribution [-]
             'Par_alb': 0.12,  # shoot Par-albedo [-]
             'Nir_alb': 0.55,  # shoot NIR-albedo [-]
             'leaf_emi': 0.98  # leaf emissivity [-] 
             }

# --- interception ---
interception = {'wmax': 0.2,  # maximum interception storage capacity for rain [kg m-2 per unit of LAI] 
                'wmaxsnow': 0.8,  # maximum interception storage capacity for snow [kg m-2 per unit of LAI]
                'w_ini': 0.0,  # initial canopy storage [kg m-2]
                'Tmin': 0.0,  # temperature below which all is snow [degC]
                'Tmax': 2.0,  # temperature above which all is water [degC]
                'leaf_orientation': 0.5, # leaf orientation factor for randomdly oriented leaves
                }

# --- define two planttypes ---

z = np.linspace(0, grid['zmax'], grid['Nlayers'])  # grid [m] above ground

pt1 = { 'name': 'pine',
        'LAImax': 2.1, # maximum annual LAI m2m-2
        'lad': lad_weibul(z, LAI=1.0, h=15.0, hb=3.0, species='pine'),  # leaf-area density m2m-3
        # cycle of photosynthetic activity
        'phenop': {
            'Xo': 0.0,
            'fmin': 0.1,
            'Tbase': -4.67,  # Kolari 2007
            'tau': 8.33,  # Kolari 2007
            'smax': 18.0  # Kolari 2014
            },
        # cycle of LAI
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
        # A-gs model
        'photop': {
            'Vcmax': 55.0,
            'Jmax': 105.0,  # 1.97*Vcmax (Kattge and Knorr, 2007)
            'Rd': 1.3,  # 0.023*Vcmax
            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                'Vcmax': [78., 200., 649.],
                'Jmax': [56., 200., 646.],
                'Rd': [33.0]
                },
            'alpha': 0.2,   # quantum efficiency parameter -
            'theta': 0.7,   # curvature parameter
            'g1': 2.3,      # stomatal slope kPa^(0.5)
            'g0': 5.0e-3,   # residual conductance mol m-2 s-1
            'kn': 0.5,      # nitrogen attenuation coefficient -
            'beta': 0.95,   # co-limitation parameter -
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.02,     # leaf length scale m
            },
        # root zone
        'rootp': {
            'root_depth': 0.5, # rooting depth [m]
            'beta': 0.943, # root distribution shape [-]
            'RAI_LAI_multiplier': 2.0, # fine-root to leaf-area ratio [-]
            'fine_radius': 2.0e-3, # [m]
            'root_cond': 5.0e8, # [s]
            }
        }


pt2 = { 'name': 'spruce',
        'LAImax': 1.0, # maximum annual LAI m2m-2
        'lad': lad_weibul(z, LAI=1.0, h=15.0, hb=0.5, species='spruce'),  # leaf-area density m2m-3
        # cycle of photosynthetic activity
        'phenop': {
            'Xo': 0.0,
            'fmin': 0.1,
            'Tbase': -4.67,  # Kolari 2007
            'tau': 8.33,  # Kolari 2007
            'smax': 18.0  # Kolari 2014
            },
        # cycle of LAI
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
        # A-gs model
        'photop': {
            'Vcmax': 60.0,
            'Jmax': 114.0,  # 1.97*Vcmax (Kattge and Knorr, 2007)
            'Rd': 1.5,  # 0.023*Vcmax
            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                'Vcmax': [53., 202., 640.],
                'Jmax': [56., 202., 656.],
                'Rd': [33.0]
                },
            'alpha': 0.2,   # quantum efficiency parameter -
            'theta': 0.7,   # curvature parameter
            'g1': 2.3,      # stomatal slope kPa^(0.5)
            'g0': 5.0e-3,   # residual conductance mol m-2 s-1
            'kn': 0.5,      # nitrogen attenuation coefficient -
            'beta': 0.95,   # co-limitation parameter -
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.02,     # leaf length scale m
            },
        # root zone
        'rootp': {
            'root_depth': 0.5, # rooting depth [m]
            'beta': 0.943, # root distribution shape [-]
            'RAI_LAI_multiplier': 2.0, # fine-root to leaf-area ratio [-]
            'fine_radius': 2.0e-3, # [m]
            'root_cond': 5.0e8, # [s]
            }
        }
        
pt3 = { 'name': 'decid',
        'LAImax': 1.2, # maximum annual LAI m2m-2
        'lad': lad_weibul(z, LAI=1.0, h=10.0, hb=0.5, species='birch'),  # leaf-area density m2m-3
        # cycle of photosynthetic activity
        'phenop': {
            'Xo': 0.0,
            'fmin': 0.1,
            'Tbase': -4.67,  # Kolari 2007
            'tau': 8.33,  # Kolari 2007
            'smax': 18.0  # Kolari 2014
            },
        # annual cycle of LAI
        'laip': {
            'lai_min': 0.1, # relative to LAImax
            'lai_ini': None,
            'DDsum0': 0.0,
            'Tbase': 5.0,
            'ddo': 45.0,
            'ddmat': 250.0,
            'sdl': 12.0,
            'sdur': 30.0
            },
        # A-gs model
        'photop': {
            'Vcmax': 50.0,
            'Jmax': 95.0,  # 1.97*Vcmax (Kattge and Knorr, 2007)
            'Rd': 1.3,  # 0.023*Vcmax
            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                'Vcmax': [77., 200., 637.],
                'Jmax': [43., 200., 637.],
                'Rd': [33.0]
                },
            'alpha': 0.2,   # quantum efficiency parameter -
            'theta': 0.7,   # curvature parameter
            'g1': 4.5,      # stomatal slope kPa^(0.5)
            'g0': 5.0e-3,   # residual conductance mol m-2 s-1
            'kn': 0.5,      # nitrogen attenuation coefficient -
            'beta': 0.95,   # co-limitation parameter -
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.05,     # leaf length scale m
            },
        # root zone
        'rootp': {
            'root_depth': 0.5, # rooting depth [m]
            'beta': 0.943, # root distribution shape [-]
            'RAI_LAI_multiplier': 2.0, # fine-root to leaf-area ratio [-]
            'fine_radius': 2.0e-3, # [m]
            'root_cond': 5.0e8, # [s]
            }
        }

pt4 = { 'name': 'shrubs',
        'LAImax': 0.7, # maximum annual LAI m2m-2
        'lad': lad_constant(z, LAI=1.0, h=0.5, hb=0.0),  # leaf-area density m2m-3
        # cycle of photosynthetic activity
        'phenop': {
            'Xo': 0.0,
            'fmin': 0.1,
            'Tbase': -4.67,  # Kolari 2007
            'tau': 8.33,  # Kolari 2007
            'smax': 18.0  # Kolari 2014
            },
        # annual cycle of LAI
        'laip': {
            'lai_min': 0.1, # relative to LAImax
            'lai_ini': None,
            'DDsum0': 0.0,
            'Tbase': 5.0,
            'ddo': 45.0,
            'ddmat': 250.0,
            'sdl': 12.0,
            'sdur': 30.0
            },
        # A-gs model
        'photop': {
            'Vcmax': 40.0,
            'Jmax': 76.0,  # 1.97*Vcmax (Kattge and Knorr, 2007)
            'Rd': 0.7,  # 0.023*Vcmax
            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                'Vcmax': [77., 200., 637.],
                'Jmax': [43., 200., 637.],
                'Rd': [33.0]
                },
            'alpha': 0.2,   # quantum efficiency parameter -
            'theta': 0.7,   # curvature parameter
            'g1': 4.5,      # stomatal slope kPa^(0.5)
            'g0': 5.0e-3,   # residual conductance mol m-2 s-1
            'kn': 0.5,      # nitrogen attenuation coefficient -
            'beta': 0.95,   # co-limitation parameter -
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.05,     # leaf length scale m
            },
        # root zone
        'rootp': {
            'root_depth': 0.5, # rooting depth [m]
            'beta': 0.943, # root distribution shape [-]
            'RAI_LAI_multiplier': 2.0, # fine-root to leaf-area ratio [-]
            'fine_radius': 2.0e-3, # [m]
            'root_cond': 5.0e8, # [s]
            }
        }
        
""" --- forestfloor --- """

snowpack = {
        'kmelt': 2.31e-5,  # Melting coefficient [kg m-2 s-1 degC-1] (=2.0 mm/C/d)
        'kfreeze': 5.79e-6,  # Freezing  coefficient [kg m-2 s-1 degC-1] (=0.5 mm/C/d)
        'retention': 0.2,  # max fraction of liquid water in snow [-]
        'Tmelt': 0.0,  # temperature when melting starts [degC]
        'optical_properties': {
                'emissivity': 0.97,
                'albedo': {'PAR': 0.8, 'NIR': 0.8}
                },
        'initial_conditions': {'temperature': 0.0,
                               'snow_water_equivalent': 0.0}
        }

soil_respiration = {
        'r10': 2.5, # [umol m-2 s-1]
        'q10': 2.0, # [-]
        'moisture_coeff': [3.83, 4.43, 1.25, 0.854]  # Skopp moisture function param [a ,b, d, g]}
        }

Forest_moss = {
    'name': 'forest mosses',  # Hylocomium splendens and Pleurozium schreberi
    'layer_type': 'bryophyte',
    'coverage': 1.0,
    'height': 0.057,  # range: [0.021, 0.10]
    'roughness_height': 0.01,
    #'dry_mass': 0.668,  # range: [0.484, 1.62]
    'bulk_density': 14.3,  # range: [7.3, 28.74]
    'max_water_content': 9.7,  # range: [7.91, 11.8], fitted value is 27.6
    'water_content_ratio': 0.25,  # max_symplast_water_content:max_water_content -ratio
    #'max_symplast_water_content': 2.4,  # based on fitted value of other mosses than Sphagnum
    'min_water_content': 0.1,
    'porosity': 0.98,

    'photosynthesis': { # farquhar-parameters
        'Vcmax': 15.0, 'Jmax': 28.5, 'Rd': 0.75, # umolm-2s-1
        'alpha': 0.3, 'theta': 0.8, 'beta': 0.9, # quantum yield, curvature, co-limitation
        'gmax': 0.02, 'wopt': 7.0, 'a0': 0.7, 'a1': -0.263, 'CAP_desic': [0.44, 7.0],
        'tresp': {
            'Vcmax': [69.83, 200.0, 27.56],
            'Jmax': [100.28, 147.92, 19.8],
            'Rd': [33.0]
        }
    },
    'optical_properties': {
        'emissivity': 0.98,
        'albedo': {'PAR': 0.11, 'NIR': 0.29} # albedoes when fully hydrated [-]
    },
    'water_retention': {
        # 'theta_s': 0.526,  # based on fitted value of other mosses than Sphagnum
        # 'theta_r': 0.07,  # based on fitted value of other mosses than Sphagnum
        'alpha': 0.166,  # based on fitted value of other mosses than Sphagnum
        'n': 1.679,  # based on fitted value of other mosses than Sphagnum
        'saturated_conductivity': 1.17e-8,  # [m s-1], based on fitted value of other mosses than Sphagnum
        'pore_connectivity': -2.30, # based on fitted value of other mosses than Sphagnum
    },
    'initial_conditions': {
        'temperature': 10.0,
        'water_content': 10.0
    }

}

# --- compile forestfloor parameter dictionary

forestfloor = {
    'bottom_layer_types': {'forest_moss': Forest_moss},
    'snowpack': snowpack,
    'soil_respiration': soil_respiration
}

# --- compile canopy-model parameter dictionary

cpara = {'loc': loc,
         'ctr': ctr,
         'grid': grid,
         'radiation': radiation,
         'micromet': micromet,
         'interception': interception,
         'planttypes': {'pine': pt1, 'spruce': pt2, 'decid': pt3, 'shrubs': pt4},
         'forestfloor': forestfloor
         }

""" --- Soil submodel parameters --- """

# grid and soil properties: pF and conductivity from Launiainen et al. 2015 Hyytiala

soil_grid = {#thickness of computational layers [m]
            'dz': [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
                   0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
                   0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                   0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            # bottom depth of layers with different characteristics [m]
            'zh': [-0.05, -0.11, -0.35, -10.0]
            }

soil_properties = {'pF': {  # vanGenuchten water retention parameters
                        'ThetaS': [0.80, 0.50, 0.50, 0.41],
                        'ThetaR': [0.01, 0.08, 0.08, 0.03],
                        'alpha': [0.70, 0.06, 0.06, 0.05],
                        'n': [1.25, 1.35, 1.35, 1.21]
                        },
                  'saturated_conductivity_vertical': [2.42E-05, 2.08e-06, 3.06e-06, 4.17e-06],  # saturated vertical hydraulic conductivity [m s-1]
                  'saturated_conductivity_horizontal': [2.42E-05, 2.08e-06, 3.06e-06, 4.17e-06],  # saturated horizontal hydraulic conductivity [m s-1]
                  'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
                  'solid_composition': {
                         'organic': [0.1611, 0.0714, 0.1091, 0.028],
                         'sand': [0.4743, 0.525, 0.5037, 0.5495],
                         'silt': [0.3429, 0.3796, 0.3641, 0.3973],
                         'clay': [0.0217, 0.0241, 0.0231, 0.0252]
                         },
                  'freezing_curve': [0.2, 0.5, 0.5, 0.5],  # freezing curve parameter
                  'bedrock': {
                              'solid_heat_capacity': 2.16e6,  # [J m-3 (solid) K-1]
                              'thermal_conductivity': 3.0  # thermal conductivity of non-porous bedrock [W m-1 K-1]
                              }
                  }

# --- water model specs
water_model = {'solve': True,
               'type': 'Richards',  #'Equilibrium', # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
               'pond_storage_max': 0.05,  #  maximum pond depth [m]
               'initial_condition': {  # (dict) initial conditions
                       'ground_water_level': -2.0,  # groundwater depth [m]
                       'pond_storage': 0.  # initial pond depth at surface [m]
                       },
               'lower_boundary': {  # lower boundary condition (type, value, depth)
                       'type': 'head_oneway',
                       'value': -0.0,
#                       'type': 'impermeable',
#                       'value': None,
#                       'depth': -2.0
                       },
               'drainage_equation': {  # drainage equation and drainage parameters
                       'type': None,  #
#                       'type': 'Hooghoudt',  #
#                       'depth': 1.0,  # drain depth [m]
#                       'spacing': 45.0,  # drain spacing [m]
#                       'width': 1.0,  # drain width [m]
                       }
                }

# --- heat model specs
heat_model = {'solve': True,
              'initial_condition': {
                      'temperature': 4.0,  # initial soil temperature [degC]
                      },
              'lower_boundary': {  # lower boundary condition (type, value)
                      'type': 'temperature',
                      'value': 4.0
                      },
              }

# --- soil model parameter dictionary
spara = {'grid': soil_grid,
         'soil_properties': soil_properties,
         'water_model': water_model,
         'heat_model': heat_model}
