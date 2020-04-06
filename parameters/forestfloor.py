#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 10:29:44 2017

@author: ajkieloaho

Parameters and initial conditions for ForestFloor

Note: 12.3.2020
Based on literature review, bryophyte parameters have been revised:
- hylocomium and pleurozium are combined to forest mosses
- sphagnum remains untacked
- bryophyte characteristics: height, dry mass, bulk density, and max water content are revised
- bryophyte water retention: theta_s, theta_r, alpha, n, K_sat, and l are revised
- eventhough, dry mass, theta_s and theta_r are derived from other parameters, those are included but commented out from parameters
- min water content is fitted value from Mualem-VanGenuchten water retention curve (pressure head vs. gravimetric water content) and it is equivalent of theta_r

"""

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

Litter = {
    'name': 'Litter',
    'layer_type': 'litter',
    'coverage': 0.25,  # [-]
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

# Note: renewed bryophyte parameters

# Based on literature review Pleurozium schreberi and Hylocomium splendens does not differ from each other
# median and range [min, max]
Forest_moss = {
    'name': 'forest mosses',  # Hylocomium splendens and Pleurozium schreberi
    'layer_type': 'bryophyte',
    'coverage': 0.5,
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

# this is general Sphagnum parametrisation based on literature review
Sphagnum = {
    'name': 'Sphagnum sp.',
    'layer_type': 'bryophyte',
    'coverage': 0.5,
    'height': 0.06, #0.06,  # range: [0.044, 0.076]
    'roughness_height': 0.02,
    # 'dry_mass': 1.49,  # range: [0.592, 2.43]
    'bulk_density': 35.1,  # range: [9.28, 46.7]
    'max_water_content': 17.8,  # range: [15.6, 24.4]
    'water_content_ratio': 0.43,  # max_symplast_water_content:max_water_content -ratio
    #'max_symplast_water_content': 7.64,  # based on fitted value of Sphagnum
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
        # 'theta_s': 0.679,  # based on fitted value
        # 'theta_r': 0.176,  # based on fitted value
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

S_fuscum = {
    'name': 'Sphagnum fuscum',
    'layer_type': 'bryophyte',
    'coverage': 0.25,
    'height': 0.045,    # Soudziloskaia et al (2013)
    'roughness_height': 0.01,
    'bulk_density': 35.1,  # Soudziloskaia et al (2013)
    'max_water_content': 20.0,
    'min_water_content': 1.5,
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
        # 'theta_s': 0.63,  # max_water_content / WATER_DENSITY * bulk_density
        # 'theta_r': 0.10,  # min_water_content / WATER_DENSITY * bulk_density
        'alpha': 0.10,
        'n': 1.4,
        'saturated_conductivity': 3.5e-4,  # [m s-1]
        'pore_connectivity': -2.37
    },
    'initial_conditions': {
        'temperature': 10.0,
        'water_content': 20.0
    }
}


forestfloor = {
    'bottom_layer_types': {
        'litter': Litter,
        'forest_moss': Forest_moss,
        'sphagnum': Sphagnum,
    },
    'snowpack': snowpack,
    'soil_respiration': soil_respiration
}

#%% for empirical photosynthetic light-response - case

#Pleurozium = {
#        'name': 'Pleurozium schreberi',  # (Brid.) Mitt.
#        'layer_type': 'bryophyte',
#        'coverage': 0.0,
#        'height': 0.095,  # Soudziloskaia et al (2013)
#        'roughness_height': 0.01,
#        'bulk_density': 17.1,  # Soudziloskaia et al (2013)
#        'max_water_content': 10.0,
#        'min_water_content': 1.5,
#        'porosity': 0.98,
#        'photosynthesis': { # per what?
#            'amax': 1.8,  # check from Rice et al. 2011 Bryologist Table 1
#            'b': 150.0,  # check from Rice et al. 2011 Bryologist Table 1
#            'moisture_coeff': [6.4355, -14.0605, 9.1867, -0.8720],
#            'temperature_coeff': [-4.3e-5, -8.3e-4, 0.08, 0.1]
#        },
#        'respiration': {# per what?
#            'q10': 2.0,  # check from Rice et al. 2011 Bryologist Table 1
#            'r10': 0.7  # check from Rice et al. 2011 Bryologist Table 1
#        },
#        'optical_properties': {  # [0.1102, 0.2909, 0.98]
#            'emissivity': 0.98,
#            'albedo': {'PAR': 0.11, 'NIR': 0.29} # albedos when fully hydrated [-]
#        },
#        'water_retention': {
#            # 'theta_s': 0.17,
#            # 'theta_r': 0.026,
#            'alpha': 0.13,
#            'n': 2.17,
#            'saturated_conductivity': 1.16e-8,  # [m s-1]
#            'pore_connectivity': -2.37
#        },
#        'initial_conditions': {'temperature': 10.0,
#                               'water_content': 10.0}
#        }

#        'water_retention_parameters': [0.445, 0.02, 0.103, 2.229, 1.17e-4, -2.487],
#Hylocomium = {
#        'name': 'Hylocomium splendens',  # (Hedw.) B.S.G.
#        'layer_type': 'bryophyte',
#        'coverage': 0.8,
#        'height': 0.06,  # Soudziloskaia et al (2013)
#        'roughness_height': 0.01,
#        #'leaf_area_index': 1.212,
#        #'specific_leaf_area': 145.0,  # Bond-Lamberty and Gower (2007)
#        #'dry_mass': bulk_density x height
#        'bulk_density': 14.3,  # kg m-3 Soudziloskaia et al (2013)
#        'max_water_content': 10.0, # g g-1
#        'min_water_content': 1.5,
#        'porosity': 0.98,
#        'photosynthesis': {
#                'amax': 1.8,  # check from Rice et al. 2011 Bryologist Table 1
#                'b': 150.0,  # check from Rice et al. 2011 Bryologist Table 1
#                'moisture_coeff': [6.4355, -14.0605, 9.1867, -0.8720],
#                'temperature_coeff': [-4.3e-5, -8.3e-4, 0.08, 0.1]
#                },
#        'respiration': { #  [2.0, 1.1]
#                'q10': 2.0,  # check from Rice et al. 2011 Bryologist Table 1
#                'r10': 0.7  # check from Rice et al. 2011 Bryologist Table 1
#                },
#        'optical_properties': { # [0.1102, 0.2909, 0.98],
#                'emissivity': 0.98,
#                'albedo': {'PAR': 0.11, 'NIR': 0.29} # albedos when fully hydrated [-]
#                },
#        'water_retention': {
#            # 'theta_s': 0.17,  # max_water_content * 1e-3 * bulk_density
#            # 'theta_r': 0.026,  # min_water_content * 1e-3 * bulk_density
#            'alpha': 0.13,
#            'n': 2.17,
#            'saturated_conductivity': 1.16e-8,  # [m s-1]
#            'pore_connectivity': -2.37
#            },
#        'initial_conditions': {'temperature': 10.0,
#                               'water_content': 10.0},
#    }

#Sphagnum = {
#    'name': 'Sphagnum fuscum',
#    'layer_type': 'bryophyte',
#    'coverage': 0.0,
#    'height': 0.045,    # Soudziloskaia et al (2013)
#    'roughness_height': 0.01,
#    #'leaf_area_index': 2.136,
#    #'specific_leaf_area': 356.0,  # Bond-Lamberty and Gower (2007)
#    #'dry_mass': 1.58,  # Soudziloskaia et al (2013)
#    'bulk_density': 35.1,  # Soudziloskaia et al (2013)
#    'max_water_content': 18.0,
#    'min_water_content': 1.5,
#    'porosity': 0.98,
#    'photosynthesis': {  # [4.0, 175.0],  # Amax, b (half-saturation rate)
#        'amax': 4.0,  # Sphagnum Rice et al. 2008 Am. J. Bot. Table 3
#        'b': 175.0,  # Sphagnum Rice et al. 2008 Am. J. Bot. Table 3
#        'moisture_coeff': [6.4355, -14.0605, 9.1867, -0.8720],
#        'temperature_coeff': [-4.3e-5, -8.3e-4, 0.08, 0.1]
#    },
#    'respiration': {
#        'q10': 2.0,  # Rice et al. 2008 Am. J. Bot. Table 3
#        r10': 0.72  # Rice et al. 2008 Am. J. Bot. Table 3
#    },
#    'optical_properties': { #  [0.0975, 0.2674, 0.98]
#        'emissivity': 0.98,
#        'albedo': {'PAR': 0.10, 'NIR': 0.27} # albedos when fully hydrated [-]
#    },
#    'water_retention': {
#        # 'theta_s': 0.63,  # max_water_content / WATER_DENSITY * bulk_density
#        # 'theta_r': 0.10,  # min_water_content / WATER_DENSITY * bulk_density
#        'alpha': 0.10,
#        'n': 1.4,
#        'saturated_conductivity': 3.5e-4,  # [m s-1]
#        'pore_connectivity': -2.37
#    },
#    'initial_conditions': {
#        'temperature': 10.0,
#        'water_content': 4.0
#    }
#}


#Community-scale light-responses and respiation rates (umol m-2 s-1):

#    Sphagnum Rice et al. 2008 Am. J. Bot. Table 3; recomputed R10 assuming Q10 = 2.0
#    Amax 2.08...4.65
#    Rd10 0.14...0.22 x Amax
#    b = 150 ... 200 umolm-2s-1 NOTE I changed code so that b is now given directly as half-saturation rate

#    Pleurozium Rice et al. 2011 Bryologist Table 1
#    Amax = 4.97 (intro gives range 1...4.8 umolm-2s-1 and references)
#    Rd  0.14 x Amax


# -- Example of bryophyte parameters ---
# moss = {
#    'species': Scientific name
#    'ground_coverage': [0.0 1.0] sum of all forest floor elements have to bee 1.0!
#    'height': [m]
#    'roughness_height': [m]
#    'leaf_area_index': [m2 m-2] photosynthesis
#    'specific_leaf_area': [m3 m-3] photosynthesis
#    'dry_mass': [kg m-2] calculated from bulk_density * height NOT NEEDED
#    'bulk_density': [kg m-3]
#    'max_water_content': [g g-1] at field capacity! h = -0.01 m
#    'min_water_content': [g g-1] at air dry! h = -1000 m
#    'porosity': total pore volume [m3 m-3]
#    'photosynthesis': {
#        'Amax': [umol m-2 s-1]
#        'b': [mol mol-1]
#    },
#    'respiration': {
#        'Q10': [-]
#        'Rd10': [umol m-2 s-1]
#    },
#    'optical_properties': {
#        'emissivity':
#        'albedo_PAR':
#        'albedo_NIR':
#    },
#    'water_retention': {
#        # theta_s is max_water_content!
#        'theta_s': calculated from height, bulk and water densities NOT NEEDED
#        # theta_r is min_water_content!
#        'theta_r': calculated from height, bulk and water densities NOT NEEDED
#        'alpha': [cm-1]
#        'n': [-]
#        'saturated_conductivity': [m s-1]
#        'pore_connectivity': [-]
#    }
#}

#Pleurozium = {
#        'name': 'Pleurozium schreberi',  # (Brid.) Mitt.
#        'layer_type': 'bryophyte',
#        'coverage': 0.25,
#        'height': 0.095,  # Soudziloskaia et al (2013)
#        'roughness_height': 0.01,
#        'bulk_density': 17.1,  # Soudziloskaia et al (2013)
#        'max_water_content': 10.0,
#        'min_water_content': 1.5,
#        'porosity': 0.98,
#
#        'photosynthesis': { # farquhar-parameters
#            'Vcmax': 15.0, 'Jmax': 28.5, 'Rd': 0.75, # umolm-2s-1
#              'alpha': 0.3, 'theta': 0.8, 'beta': 0.9, # quantum yield, curvature, co-limitation
#              'gmax': 0.02, 'wopt': 7.0, 'a0': 0.7, 'a1': -0.263, 'CAP_desic': [0.44, 7.0],
#              'tresp': {'Vcmax': [69.83, 200.0, 27.56],
#                        'Jmax': [100.28, 147.92, 19.8],
#                        'Rd': [33.0]
#                       }
#              },
#        'optical_properties': {
#            'emissivity': 0.98,
#            'albedo': {'PAR': 0.11, 'NIR': 0.29} # albedos when fully hydrated [-]
#        },
#        'water_retention': {
#            # 'theta_s': 0.17,
#            # 'theta_r': 0.026,
#            'alpha': 0.13,
#            'n': 2.17,
#            'saturated_conductivity': 1.16e-8,  # [m s-1]
#            'pore_connectivity': -2.37
#        },
#        'initial_conditions': {'temperature': 10.0,
#                               'water_content': 10.0}
#        }

#Hylocomium = {
#        'name': 'Hylocomium splendens',  # (Hedw.) B.S.G.
#        'layer_type': 'bryophyte',
#        'coverage': 0.25,
#        'height': 0.06,  # Soudziloskaia et al (2013)
#        'roughness_height': 0.01,
#        #'leaf_area_index': 1.212,
#        #'specific_leaf_area': 145.0,  # Bond-Lamberty and Gower (2007)
#        #'dry_mass': bulk_density x height
#        'bulk_density': 14.3,  # kg m-3 Soudziloskaia et al (2013)
#        'max_water_content': 10.0, # g g-1
#        'min_water_content': 1.5,
#        'porosity': 0.98,
#
#        'photosynthesis': { # farquhar-parameters
#            'Vcmax': 15.0, 'Jmax': 28.5, 'Rd': 0.75, # umolm-2s-1
#              'alpha': 0.3, 'theta': 0.8, 'beta': 0.9, # quantum yield, curvature, co-limitation
#              'gmax': 0.02, 'wopt': 7.0, 'a0': 0.7, 'a1': -0.263, 'CAP_desic': [0.44, 7.0],
#              'tresp': {'Vcmax': [69.83, 200.0, 27.56],
#                        'Jmax': [100.28, 147.92, 19.8],
#                        'Rd': [33.0]
#                       }
#              },
#        'optical_properties': {
#                'emissivity': 0.98,
#                'albedo': {'PAR': 0.11, 'NIR': 0.29} # albedos when fully hydrated [-]
#                },
#        'water_retention': {
#            # 'theta_s': 0.17,  # max_water_content * 1e-3 * bulk_density
#            # 'theta_r': 0.026,  # min_water_content * 1e-3 * bulk_density
#            'alpha': 0.13,
#            'n': 2.17,
#            'saturated_conductivity': 1.16e-8,  # [m s-1]
#            'pore_connectivity': -2.37
#            },
#        'initial_conditions': {'temperature': 10.0,
#                               'water_content': 10.0},
#    }



# pf_xero = [0.445, 0.02, 0.103, 2.229, 1.17e-4, -2.487]
# pf_shpag = [0.95, 0.10, 0.34, 1.4, 3.5e-4, -4.38]