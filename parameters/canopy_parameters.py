# -*- coding: utf-8 -*-
"""
CANOPY MODEL PARAMETERS
"""

# initialize dictionary to store parameters
cpara = {}

# --- control flags ---
ctr = {'multilayer_model': False,  # compute in multilayer mode
       'seasonal_LAI': True,  # account for seasonal LAI dynamics
       'pheno_cylcle': True  # account for phenological cycle
       }

# --- site location ---
loc = {'lat': 61.4,  # latitude
       'lon': 23.7  # longitude
       }

# --- aerodynamics ---
aero = {'w': 0.01,  # leaf length scale [m]
        'zmeas': 2.0,  # wind speed measurement height above canopy [m]
        'zg': 0.5,  # height above ground where Ug is computed [m]
        'zos': 0.01  # forest floor roughness length [m]
        }

# --- interception and snowmodel ---  SADANNAN KORJAUSKERTOIMET?
interc_snow = {'wmax': 0.0005,  # maximum interception storage capacity for rain [m per unit of LAI]
               'wmaxsnow': 0.004,  # maximum interception storage capacity for snow [m per unit of LAI]
               'kmelt': 2.8934e-08,  # Melting coefficient [m degC-1 s-1]
               'kfreeze': 5.79e-09,  # Freezing  coefficient [m degC-1 s-1]
               'retention': 0.05,  # max fraction of liquid water in snow [-]
               'Tmin': 0.0,  # temperature below which all is snow [degC]
               'Tmax': 1.0,  # temperature above which all is water [degC]
               'Tmelt': 0.0,  # temperature when melting starts [degC]
               'w_ini': 0.0,  # initial canopy storage [m]
               'swe_ini': 0.0  # initial snow water equivalent [m]
               }

# --- default values for plant phenological cycle and lai seasonality ---
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
                         }
                }
pine, spruce, decid = plant_default, plant_default, plant_default  # initialize with default

# --- stand characteristics ---
# specify name and maximum leaf-area index, LAImax [m2/m2], and 
# adjust values of 'phenop' and 'laip' if default values not suitable
pine.update({'name': 'pine', 'LAImax': [2.1]})
pine['phenop'].update({'fmin': 0.1})
pine['laip'].update({'lai_min': 0.8})

spruce.update({'name': 'spruce', 'LAImax': [1.0]})
spruce['phenop'].update({'fmin': 0.1})
spruce['laip'].update({'lai_min': 0.8})

decid.update({'name': 'decid', 'LAImax': [1.0]})

if ctr['multilayer_model'] is False:
    """parameters for simple (not multilayer) model"""
    # --- physiology for calculation of transpiration and Efloor---
    phys_para = {'ga': 40.0,  # [m s-1]          MISSÄ KÄYTETÄÄN?
                 'q50': 50.0,  # half-sat. of leaf light response [W m-2]
                 'kp': 0.6,  # parameter to calculate fraction of Rnet at ground [-]
                 'f': 0.8,  # fraction of local Rnet available for evaporation at ground [-]
                 'rw': 0.20,  # transpiration moisture response parameter
                 'rwmin': 0.02}  # transpiration moisture response parameter
    pine['gsref'] = 1.6e-6  # [m s-1]
    spruce['gsref'] = 1.6e-6
    decid['gsref'] = 2.6e-6
    # --- parameters describing canopy ---
    canopy_para = {'hc': 16.0,  # canopy height [m]
                   'cf': 0.6}  # canopy closure [-]
    # add to parameter dictionary
    cpara.update({'phys_para':phys_para,'canopy_para':canopy_para})
    # plant types to list
    plant_types = [pine, spruce, decid]
else:
     """parameters for multilayer model"""
    # define shrubs, 'photop', 'leafp', 'lad', grid

cpara.update({'ctr': ctr, 'loc': loc, 'aero': aero,
              'interc_snow': interc_snow, 'plant_types': plant_types})

#        for computing aerodynamic resistances  -- yksiköt?
#        self.zmeas = 2.0
#        self.zground = 0.5
#        self.zo_ground = 0.01
#        self.soilrp = 300.0