# -*- coding: utf-8 -*-
"""
CANOPY MODEL PARAMETERS
"""

cpara = {
        # --- site location ---
        'lat': 61.4,  # latitude
        'lon': 23.7,  # longitude
        # --- stand characteristics ---
        'lai_conif': 4.0,  # LAI conifers [m2 m-2]
        'lai_decid': 0.5,  # LAI deciduous [m2 m-2]
        'hc': 16.0,  # canopy height [m]
        'cf': 0.6,  # canopy closure [-]
        # --- interception and snowmodel ---
        'wmax': 0.0005,  # maximum interception storage capacity for rain [m per unit of LAI]
        'wmaxsnow': 0.004,  # maximum interception storage capacity for snow [m per unit of LAI]
        'kmelt': 2.8934e-08,  # Melting coefficient [m degC-1 s-1]
        'kfreeze': 5.79e-09,  # Freezing  coefficient [m degC-1 s-1]
        'retention': 0.05,  # max fraction of liquid water in snow [-]
        'Tmin': 0.0,  # temperature below which all is snow [degC]
        'Tmax': 1.0,  # temperature above which all is water [degC]
        'Tmelt': 0.0,  # temperature when melting starts [degC]
        # --- tree physiology ---
        'ga': 40.0,  # [m s-1]          MISSÄ KÄYTETÄÄN?
        'q50': 50.0,  # half-sat. of leaf light response [W m-2]
        'gsref_conif': 1.6e-6,  # [m s-1]
        'gsref_decid': 2.6e-6,  # [m s-1]
        'kp': 0.6,  # parameter to calculate fraction of Rnet at ground [-]
        'f': 0.8,  # fraction of local Rnet available for evaporation at ground [-]
        'rw': 0.20,  # transpiration moisture response parameter
        'rwmin': 0.02,  # transpiration moisture response parameter
        # --- seasonal cycle of physiology ---
        'smax': 18.5,  # [degC]
        'tau': 13.0,  # [d]
        'xo': -4.0,  # [degC]
        'fmin': 0.05,  # residual photocapasity [-]
        # --- deciduous leaf-area dynamics ---
        'lai_min': 0.1,  # [-]
        'ddo': 45.0,  # degreedays for bud burst [days]
        'ddur': 23.0,  # [days]
        'sdl': 9.0,  # senescence, critical daylength [hours]
        'sdur': 30.0,  # [days]
        # --- initial state ---
        'w': 0.0,  # canopy storage [m]
        'swe': 0.0  # snow water equivalent [m]
        }

#        PARAMETREIHIN MYÖS??

#        SATEEN LAATU
#        Tmin = 0.0  # 'C, below all is snow
#        Tmax = 1.0  # 'C, above all is water
#        Tmelt = 0.0  # 'C, T when melting starts

#        SADANNAN KORJAUSKERTOIMET?

#        for computing aerodynamic resistances  -- yksiköt?
#        self.zmeas = 2.0
#        self.zground = 0.5
#        self.zo_ground = 0.01
#        self.soilrp = 300.0