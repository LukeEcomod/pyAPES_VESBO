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
        # --- interception and snowmodel --- Kmelt, Kfreeze [mm s-1]
        'wmax': 1.3,  # maximum interception storage capacity for rain [kg / unit of LAI]
        'wmaxsnow': 6.0,  # maximum interception storage capacity for snow [kg / unit of LAI]
        'kmelt': 2.8934e-05,  # Melting coefficient [mm s-1]
        'kfreeze': 5.79e-6,  # Freezing  coefficient [mm s-1]
        # --- tree physiology: ga ,q50 [Wm-2 (global)], gsref , kp [-]
        'ga': 40.0,  # [ms-1]          MISSÄ KÄYTETÄÄN?
        'q50': 50.0,  # half-sat. of leaf light response [Wm-2]
        'gsref_conif': 2.5e-3,  # [mm s-1]
        'gsref_decid': 3.75e-3,  # [mm s-1]
        'kp': 0.6,  # parameter to calculate fraction of Rnet at ground [-]
        'f': 0.8,  # fraction of local Rnet available for evaporation at ground [-]
        'rw': 0.20,  # transpiration moisture response parameter
        'rwmin': 0.02,  # transpiration moisture response parameter
        # --- seasonal cycle of physiology --- smax , , xo,fmin[-](residual photocapasity)
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
        'w': 0.0,  # canopy storage [mm]
        'swe': 0.0  # snow water equivalent [mm]
        }


#        # quality of precipitation            MIKSEI OLE PARAMETREJA?
#        Tmin = 0.0  # 'C, below all is snow
#        Tmax = 1.0  # 'C, above all is water
#        Tmelt = 0.0  # 'C, T when melting starts

#        ENTÄ SADANNAN KORJAUSKERTOIMET?
