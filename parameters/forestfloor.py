#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 10:29:44 2017

@author: ajkieloaho
"""

initial_conditions = {
        'snowpack': {
                'snow_water_equivalent': 0.0
                },
        'baresoil': {
                'temperature': 10.0
                },
        'litter': {
                'temperature': 10.0,
                'water_content': 2.0
                },
        'bryophytes': {
                'hylocomium': {
                        'temperature': 10.0,
                        'water_content': 8.0
                        },
                'pleurozium': {
                        'temperature': 10.0,
                        'water_content': 8.0
                        },
                'sphagnum': {
                        'temperature': 10.0,
                        'water_content': 8.0
                        }
                }
        }

snowpack = {
        'kmelt': 2.31e-8,  # Melting coefficient [m degC-1 s-1] (=2.0 mm/C/d)
        'kfreeze': 5.79e-9,  # Freezing  coefficient [m degC-1 s-1] (=0.5 mm/C/d)
        'retention': 0.2,  # max fraction of liquid water in snow [-]
        'Tmelt': 0.0,  # temperature when melting starts [degC]
        'swe_ini': 0.0,  # initial snow water equivalent [m],
        'optical_properties': {
                'emissivity': 0.97,
                'albedo_PAR': 0.8,
                'albedo_NIR': 0.8,
                }
        }

baresoil = {
        'ground_coverage': 0.0,
        'roughness_length':  0.002, # check right value
        'optical_properties': {
                'emissivity': 0.98,
                'albedo_PAR': 0.05,  # [-]
                'albedo_NIR': 0.5,  # [-]
                },
        'respiration': {
                'R10': 2.5,  # base heterotrophic respiration rate [umolm-2s-1]
                'Q10': 2.0,  # temperature sensitivity [-]
                'limitpara': [3.83, 4.43, 1.25, 0.854]  # Skopp respiration function param [a ,b, d, g]
                }
        }

litter = {
        'ground_coverage': 0.0,  # [-]
        'height': 0.03,  # [m]
        'roughness_height': 0.01,  # [m]
        'bulk_density': 45.0,  # [kg m\ :sup:`-3`]
        'max_water_content': 10.0, # 4.0,  # [g g\ :sup:`-1`\ ]
        'min_water_content': 1.0,  # [g g\ :sup:`-1`\ ]
        'porosity': 0.95,  # [m\ :sup:`3` m\ :sup:`-3`\ ]
        'respiration': {  # Taken from baresoil!!
                'R10': 1.6,  # base heterotrophic respiration rate [umolm-2s-1]
                'Q10': 2.0,  # temperature sensitivity [-]
                'limitpara': [3.83, 4.43, 1.25, 0.854]  # Skopp respiration function param [a ,b, d, g]
                },
        'optical_properties': {  # [0.1102, 0.2909, 0.98]
                    'emissivity': 0.98,  # [-]
                    'albedo_PAR': 0.1102,  # [-]
                    'albedo_NIR': 0.2909,  # [-]
                    },
        'water_retention': {'theta_s': 0.95,  ### = porosity
                            'theta_r': 0.01,  ### määritellään min_water_content*1e-3*bulk_density
                            'alpha': 0.13,
                            'n': 2.17,
                            'saturated_conductivity': 1.16e-8,  # [m s-1]
                            'pore_connectivity': -2.37,
                            }
            }
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
#    "species": Scientific name
#    "ground_coverage": [0.0 1.0] sum of all forest floor elements have to bee 1.0!
#    "height": [m]
#    "roughness_height": [m]
#    "leaf_area_index": [m2 m-2] photosynthesis
#    "specific_leaf_area": [m3 m-3] photosynthesis
#    "dry_mass": [kg m-2] calculated from bulk_density * height NOT NEEDED
#    "bulk_density": [kg m-3]
#    "max_water_content": [g g-1] at field capacity! h = -0.01 m
#    "min_water_content": [g g-1] at air dry! h = -1000 m
#    "porosity": total pore volume [m3 m-3]
#    "photosynthesis": {
#        'Amax': [umol m-2 s-1]
#        'b': [mol mol-1]
#    },
#    "respiration": {
#        'Q10': [-]
#        'Rd10': [umol m-2 s-1]
#    },
#    'optical_properties': {
#        'emissivity':
#        'albedo_PAR':
#        'albedo_NIR':
#    },
#    "water_retention": {
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


Pleurozium = {
        "species": "Pleurozium schreberi",  # (Brid.) Mitt.
        "ground_coverage": 0.0,
        "height": 0.095,  # Soudziloskaia et al (2013)
        "roughness_height": 0.01,
        "leaf_area_index": 1.212,
        "specific_leaf_area": 262.1,  # Bond-Lamberty and Gower (2007)
        "bulk_density": 17.1,  # Soudziloskaia et al (2013)
        "max_water_content": 10.0,
        "min_water_content": 1.5,
        "porosity": 0.98,
        "photosynthesis": {
            'Amax': 1.8,  # check from Rice et al. 2011 Bryologist Table 1
            'b': 150.  # check from Rice et al. 2011 Bryologist Table 1
        },
        "respiration": {  # [2.0, 1.1]
            'Q10': 2.0,  # check from Rice et al. 2011 Bryologist Table 1
            'Rd10': 0.7  # check from Rice et al. 2011 Bryologist Table 1
        },
        'optical_properties': {  # [0.1102, 0.2909, 0.98]
            'emissivity': 0.98,
            'albedo_PAR': 0.1102,  # [-]
            'albedo_NIR': 0.2909,  # [-]
        },
        "water_retention": {
            # 'theta_s': 0.17,  # max_water_content * 1e-3 * bulk_density
            # 'theta_r': 0.026,  # min_water_content * 1e-3 * bulk_density
            'alpha': 0.13,
            'n': 2.17,
            'saturated_conductivity': 1.16e-8,  # [m s-1]
            'pore_connectivity': -2.37
        }
}
#        "water_retention_parameters": [0.445, 0.02, 0.103, 2.229, 1.17e-4, -2.487],
Hylocomium = {
        "species": "Hylocomium splendens",  # (Hedw.) B.S.G.
        "ground_coverage":1.0,
        "height": 0.06,  # Soudziloskaia et al (2013)
        "roughness_height": 0.01,
        "leaf_area_index": 1.212,
        "specific_leaf_area": 145.0,  # Bond-Lamberty and Gower (2007)
        "dry_mass": 0.86,  # Soudziloskaia et al (2013)
        "bulk_density": 14.3,  # Soudziloskaia et al (2013)
        "max_water_content": 10.0,
        "min_water_content": 1.5,
        "porosity": 0.98,
        "photosynthesis": {
                'Amax': 1.8,  # check from Rice et al. 2011 Bryologist Table 1
                'b': 150.  # check from Rice et al. 2011 Bryologist Table 1
                },
        "respiration": { #  [2.0, 1.1]
                'Q10': 2.0,  # check from Rice et al. 2011 Bryologist Table 1
                'Rd10': 0.7  # check from Rice et al. 2011 Bryologist Table 1
                },
        "optical_properties": { # [0.1102, 0.2909, 0.98],
                'emissivity': 0.98,
                'albedo_PAR': 0.1102,  # [-]
                'albedo_NIR': 0.2909,  # [-]
                },
        "water_retention": {
            # 'theta_s': 0.17,  # max_water_content * 1e-3 * bulk_density
            # 'theta_r': 0.026,  # min_water_content * 1e-3 * bulk_density
            'alpha': 0.13,
            'n': 2.17,
            'saturated_conductivity': 1.16e-8,  # [m s-1]
            'pore_connectivity': -2.37
            }
    }

Sphagnum = {
        "species": "Sphagnum fuscum",
        "ground_coverage": 0.0,
        "height": 0.045,    # Soudziloskaia et al (2013)
        "roughness_height": 0.01,
        "leaf_area_index": 2.136,
        "specific_leaf_area": 356.0,  # Bond-Lamberty and Gower (2007)
        "dry_mass": 1.58,  # Soudziloskaia et al (2013)
        "bulk_density": 35.1,  # Soudziloskaia et al (2013)
        "max_water_content": 18.0,
        "min_water_content": 1.5,
        "porosity": 0.98,
        "photosynthesis": {  # [4.0, 175.0],  # Amax, b (half-saturation rate)
            'Amax': 4.0,  # Sphagnum Rice et al. 2008 Am. J. Bot. Table 3
            'b': 175.  # Sphagnum Rice et al. 2008 Am. J. Bot. Table 3
            },
        "respiration": {
            'Q10': 2.0,  # Rice et al. 2008 Am. J. Bot. Table 3
            'Rd10': 0.72  # Rice et al. 2008 Am. J. Bot. Table 3
            },
        "optical_properties": { #  [0.0975, 0.2674, 0.98]
            'emissivity': 0.98,
            'albedo_PAR': 0.0975,  # [-]
            'albedo_NIR': 0.2674,  # [-]
            },
#        "water_retention_parameters": [0.95, 0.10, 0.34, 1.4, 3.5e-4, -4.38]
#        "water_retention_parameters": [0.2, 0.01, 0.13, 2.17, 2.07e-4, -2.37]
        "water_retention": {
            # 'theta_s': 0.63,  # max_water_content * 1e-3 * bulkd_density
            # 'theta_r': 0.10,  # min_water_content * 1e-3
            'alpha': 0.10,
            'n': 1.4,
            'saturated_conductivity': 3.5e-4,  # [m s-1]
            'pore_connectivity': -2.37
            }
    }

# check that ground coverages sum to 1.0
forestfloor = {
        'bryophytes': {
                'hylocomium': Hylocomium,
                'pleurozium': Pleurozium,
                'sphagnum': Sphagnum},
        'litter': litter,
        'baresoil': baresoil,
        'snowpack': snowpack,
        'initial_conditions': initial_conditions
        }

# pf_xero = [0.445, 0.02, 0.103, 2.229, 1.17e-4, -2.487]
# pf_shpag = [0.95, 0.10, 0.34, 1.4, 3.5e-4, -4.38]
