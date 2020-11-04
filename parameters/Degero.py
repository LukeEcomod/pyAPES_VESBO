# -*- coding: utf-8 -*-
"""
DEFINES PARAMETERS FOR pyAPES SIMULATION at Degerö Stormyr -SITE in Sweden

to import: from parameters.Degero import gpara, cpara, spara

pyAPES output variables and logger config in parameters.outputs

Todo:
- currently testing that all works, all parameters are not properly selected
- implement reading canopy structure from tree inventory data
- use parameters.parametersets to override default parameters defined here!

"""

import numpy as np
from tools.utilities import lad_constant
from pyAPES_utilities.soiltypes.sphagnum_peat import soil_properties, zh

gpara = {
    'dt' : 1800.0,  # timestep in forcing data file [s]
    'start_time' : "2015-06-01",  # start time of simulation [yyyy-mm-dd]
    'end_time' : "2015-10-01",  # end time of simulation [yyyy-mm-dd]
    'forc_filename' : "Degero/Degero_forcing_2014-2016.dat",  # forcing data file*
    'results_directory':'results/Degero/'
}

# --- control flags (True/False) ---
ctr = {
    'Eflow': True,  # ensemble flow
    'WMA': False,  #True,  #  well-mixed assumption
    'Ebal': True,  #False,  #  computes leaf temperature by solving energy balance
    'WaterStress': 'Rew',  #'PsiL',  # Rew or PsiL or None
    'seasonal_LAI': True,  # account for seasonal LAI dynamics
    'pheno_cycle': True  # account for phenological cycle
}

# site location
loc = {
    'lat': 64.11,  # latitude
    'lon': 19.33  # longitude
}

# grid
grid = {
    'zmax': 25.0,  # heigth of grid from ground surface [m]
    'Nlayers': 101  # number of layers in grid [-]
}


# --- micrometeo ---
micromet = {
    'zos': 0.01,  # forest floor roughness length [m]  -- not used?
    'dPdx': 0.0,  # horizontal pressure gradient
    'Cd': 0.15,  # drag coefficient
    'Utop': 9.0,  # ensemble U/ustar
    'Ubot': 0.01,  # lower boundary
    'Sc': {
        'T': 2.0, 'H2O': 2.0, 'CO2': 2.0}  # Schmidt numbers
}

# --- radiation ---
radiation = {
    'clump': 0.7,  # clumping index [-]
    'leaf_angle': 1.0,  # leaf-angle distribution [-]
    'Par_alb': 0.12,  # shoot Par-albedo [-]
    'Nir_alb': 0.55,  # shoot NIR-albedo [-]
    'leaf_emi': 0.98  # leaf emissivity [-]
}

# --- interception ---
interception = {
    'wmax': 0.2,  # maximum interception storage capacity for rain [kg m-2 per unit of LAI]
    'wmaxsnow': 0.8,  # maximum interception storage capacity for snow [kg m-2 per unit of LAI]
    'w_ini': 0.0,  # initial canopy storage [kg m-2]
    'Tmin': 0.0,  # temperature below which all is snow [degC]
    'Tmax': 2.0,  # temperature above which all is water [degC]
    'leaf_orientation': 0.5, # leaf orientation factor for randomdly oriented leaves
    }

# --- define planttypes ---

z = np.linspace(0, grid['zmax'], grid['Nlayers'])  # grid [m] above ground

pt1 = {
    'name': 'sedges',
    'LAImax': 0.6, # maximum annual LAI m2m-2
    'lad': lad_constant(z, LAI=1.0, h=0.0),  # leaf-area density m2m-3
    # cycle of photosynthetic activity
    'phenop': {
        'Xo': 0.0,
        'fmin': 0.1,
        'Tbase': -4.67,  # Kolari 2007
        'tau': 8.33,  # Kolari 2007ik
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
        'Vcmax': 45.0,
        'Jmax': 85.0,  # 1.97*Vcmax (Kattge and Knorr, 2007)
        'Rd': 0.9,  # 0.023*Vcmax
        'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
            'Vcmax': [72., 200., 649.],
            'Jmax': [50., 200., 646.],
            'Rd': [33.0]
        },
        'alpha': 0.2,   # quantum efficiency parameter -
        'theta': 0.7,   # curvature parameter
        'g1': 2.1,      # stomatal slope kPa^(0.5)
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
    'initial_conditions': {
        'temperature': 0.0,
        'snow_water_equivalent': 0.0}
}

soil_respiration = {
    'r10': 2.5, # [umol m-2 s-1]
    'q10': 2.0, # [-]
    'moisture_coeff': [3.83, 4.43, 1.25, 0.854]  # Skopp moisture function param [a ,b, d, g]}
}

# Note: renewed bryophyte parameters

# Based on literature review Pleurozium schreberi and Hylocomium splendens does not differ from each other
# median and range [min, max]
Forest_moss = {
    'name': 'forest mosses',  # Hylocomium splendens and Pleurozium schreberi
    'layer_type': 'bryophyte',
    'coverage': 0.0,
    'height': 0.057,  # [0.021, 0.10]
    'roughness_height': 0.01,
    'bulk_density': 14.3,  # [7.3, 28.74]
    'max_water_content': 9.7,  # [7.91, 11.8], fitted value is 27.6
    'water_content_ratio': 0.25,  # max_symplast_water_content:max_water_content -ratio
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

# this is general Sphagnum parametrisation based on literature review
Sphagnum = {
    'name': 'Sphagnum sp.',
    'layer_type': 'bryophyte',
    'coverage': 1.0,
    'height': 0.06, # [0.044, 0.076]
    'roughness_height': 0.02,
    'bulk_density': 35.1,  # [9.28, 46.7]
    'max_water_content': 17.8,  # [15.6, 24.4]
    'water_content_ratio': 0.43,  # max_symplast_water_content:max_water_content -ratio
    'min_water_content': 0.1,
    'porosity': 0.98,

    'photosynthesis': { # farquhar-parameters
        'Vcmax': 45.0, 'Jmax': 85.5, 'Rd': 1.35, # umolm-2s-1
        'alpha': 0.3, 'theta': 0.8, 'beta': 0.9, # quantum yield, curvature, co-limitation
        'gmax': 0.04, 'wopt': 7.0, 'a0': 0.7, 'a1': -0.263, 'CAP_desic': [0.58, 10.0],
        'tresp': {
            'Vcmax': [69.83, 200.0, 27.56],
            'Jmax': [100.28, 147.92, 19.8],
            'Rd': [33.0]
        }
    },
    'optical_properties': { # moisture responses are hard-coded
        'emissivity': 0.98,
        'albedo': {'PAR': 0.10, 'NIR': 0.27} # albedos when fully hydrated [-]
    },
    'water_retention': {
        'alpha': 0.381,  # based on fitted value
        'n': 1.781,  # based on fitted value
        'saturated_conductivity': 3.4e-4,  # [m s-1], based on fitted value
        'pore_connectivity': -2.11  # based on fitted value
    },
    'initial_conditions': {
        'temperature': 10.0,
        'water_content': 20.0
    }
}

Litter = {
    'name': 'Litter',
    'layer_type': 'litter',
    'coverage': 0.0,  # [-]
    'height': 0.03,  # [m]
    'roughness_height': 0.01,  # [m]
    'bulk_density': 45.0,  # [kg m\ :sup:`-3`]
    'max_water_content': 4.0, # 4.0,  # [g g\ :sup:`-1`\ ]
    'water_content_ratio': 0.25,  # max_symplast_water_content:max_water_content -ratio
    #'max_symplast_water_content': 1.0, # [g g\ :sup:`-1`\ ]
    'min_water_content': 0.1,
    'porosity': 0.95,  # [m\ :sup:`3` m\ :sup:`-3`\ ]
    'respiration': {# Taken from baresoil!! per what?
        'q10': 1.6,  # base heterotrophic respiration rate [umolm-2s-1]
        'r10': 2.0,  # temperature sensitivity [-]
        #'moisture_coeff': [add here]
    },
    'optical_properties': {  # [0.1102, 0.2909, 0.98]
        'emissivity': 0.98,  # [-]
        'albedo': {'PAR': 0.11, 'NIR': 0.29} # albedos when fully hydrated [-]
    },
    'water_retention': {#'theta_s': 0.95,  # max_water_content / WATER_DENSITY * bulk_density
                        #'theta_r': 0.01,  # min_water_content /WATER_DENSITY * bulk_density
        'alpha': 0.13,
        'n': 2.17,
        'saturated_conductivity': 1.16e-8,  # [m s-1]
        'pore_connectivity': -2.37,
    },
    'initial_conditions': {
        'temperature': 10.0,
        'water_content': 4.0
    }
}

# --- compile forestfloor parameter dictionary

forestfloor = {
    'bottom_layer_types': {
        'litter': Litter,
        'forest_moss': Forest_moss,
        'sphagnum': Sphagnum,
    },
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
         'planttypes': {'shrubs': pt1},
         'forestfloor': forestfloor
         }


# # """ --- Soil submodel parameters: these profiles resemble Degerö Stormyr --- """

soil_grid = {
    'dz': [  # thickness of computational layers [m]
        0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
        0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
        0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
        0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    ],
    'zh': [  # bottom depth of layers with different characteristics [m]
        -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1., -1.5, -2.0
    ]
}

N = len(soil_grid['zh'])

soil_properties = {
   'pF': {  # vanGenuchten water retention parameters
           'ThetaS': [0.945, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918, 0.918],
           'ThetaR': [0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098, 0.098],
           'alpha': [0.338, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072, 0.072],
           'n': [1.402, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371, 1.371]
       },
       'saturated_conductivity_vertical': [8.99E-05, 2.98E-05, 9.86E-06, 3.27E-06, 1.08E-06, 3.58E-07, 1.19E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07],  # saturated vertical hydraulic conductivity [m s-1]
       'saturated_conductivity_horizontal': [8.99E-05, 2.98E-05, 9.86E-06, 3.27E-06, 1.08E-06, 3.58E-07, 1.19E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07, 1.16E-07],  # saturated horizontal hydraulic conductivity [m s-1]
       'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
       'solid_composition': {
           'organic': [1.0 for i in range(N)],
           'sand': [0.0 for i in range(N)],
           'silt': [0.0 for i in range(N)],
           'clay': [0.0 for i in range(N)]
       },
       'freezing_curve': [0.5 for i in range(N)],  # freezing curve parameter
       'bedrock': {
           'solid_heat_capacity': 2.16e6,  # [J m-3 (solid) K-1]
           'thermal_conductivity': 3.0  # thermal conductivity of non-porous bedrock [W m-1 K-1]
           }
       }

# # --- water model specs
water_model = {
    'solve': True,
    'type': 'Equilibrium', #'Richards', # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
    'pond_storage_max': 0.002,  #  maximum pond depth [m]
    'initial_condition': {  # (dict) initial conditions
        'ground_water_level': -0.05,  # groundwater depth [m]
        'pond_storage': 0.0  # initial pond depth at surface [m]
    },
    'lower_boundary': {  # lower boundary condition (type, value, depth)
        'type': 'impermeable',
        'value': None,
        'depth': -2.0
    },
    'drainage_equation': {  # drainage equation and drainage parameters
        'type': 'Hooghoudt',  #
        'depth': 0.0,  # drain depth [m]
        'spacing': 100.0,  # drain spacing [m]
        'width': 1.0,  # drain width [m]
    }
}

# # --- heat model specs
heat_model = {
    'solve': True,
    'initial_condition': {
        'temperature': [
            9.6, 9.6, 9.6, 9.6, 9.6, 9.6,
            9.08175, 9.08175, 9.08175, 9.08175, 9.08175,
            8.109, 8.109,
            3.45675, 3.45675, 3.45675, 3.45675, 3.45675, 3.45675, 3.45675, 3.45675,
            2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122,
            2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122, 2.122
        ],  # initial soil temperature [degC]
    },
    'lower_boundary': {  # lower boundary condition (type, value)
        'type': 'temperature',
        'value': 4.0
    },
}

# # --- soil model parameter dictionary
spara = {'grid': soil_grid,
          'soil_properties': soil_properties,
          'water_model': water_model,
          'heat_model': heat_model}

# """ --- Soil submodel parameters for 1-layer soil where state is updated from
# forcing file in pyAPES.driver """

# soil_grid = {#thickness of computational layers [m]
#             'dz': [0.10],
#             # bottom depth of layers with different characteristics [m]
#             'zh': [-0.10]
#             }

# soil_properties = {
#     'pF': {  # vanGenuchten water retention parameters
#         'ThetaS': [0.945],
#         'ThetaR': [0.098],
#         'alpha': [0.338],
#         'n': [1.402]
#     },
#     'saturated_conductivity_vertical': [30*8.99E-05],  # saturated vertical hydraulic conductivity [m s-1]
#     'saturated_conductivity_horizontal': [30*8.99E-05],  # saturated horizontal hydraulic conductivity [m s-1]
#     'solid_heat_capacity': None,  # [J m-3 (solid) K-1] - if None, estimated from organic/mineral composition
#     'solid_composition': {
#         'organic': [1.0],
#         'sand': [0.0],
#         'silt': [0.0],
#         'clay': [0.0]
#     },
#     'freezing_curve': [0.5],  # freezing curve parameter
#     'bedrock': {
#         'solid_heat_capacity': 2.16e6,  # [J m-3 (solid) K-1]
#         'thermal_conductivity': 3.0  # thermal conductivity of non-porous bedrock [W m-1 K-1]
#     }
# }

# # --- water model specs
# water_model = {
#     'solve': False,
#     'type': 'Equilibrium', #'Richards', # solution approach 'Equilibrium' for equilibrium approach else solves flow using Richards equation
#     'pond_storage_max': 0.0,  #  maximum pond depth [m]
#     'initial_condition': {  # (dict) initial conditions
#         'volumetric_water_content': 0.5, #[m3m-3]
#         'ground_water_level': -0.05,  # groundwater depth [m]
#         'pond_storage': 0.0  # initial pond depth at surface [m]
#     },
#     'lower_boundary': {  # lower boundary condition (type, value, depth)
#         'type': 'head_oneway',
#         'value': -0.0,
#         #'type': 'impermeable',
#         #'value': None,
#         #'depth': -2.0
#     },
#     'drainage_equation': {  # drainage equation and drainage parameters
#         'type': 'Hooghoudt',  #
#         'depth': 0.0,  # drain depth [m]
#         'spacing': 100.0,  # drain spacing [m]
#         'width': 1.0,  # drain width [m]
#     }
# }

# # --- heat model specs
# heat_model = {
#     'solve': False,
#     'initial_condition': {
#         'temperature': 4.0,  # initial soil temperature [degC]
#     },
#     'lower_boundary': {  # lower boundary condition (type, value)
#         'type': 'temperature',
#         'value': 4.0
#     },
# }

