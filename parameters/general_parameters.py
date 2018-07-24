# -*- coding: utf-8 -*-
"""
GENERAL PARAMETERS
"""

gpara = {
        'dt' : 1800.0,  # timestep in forcing data file [s]
        'start_time' : "2011-01-01",  # start time of simulation [yyyy-mm-dd]
        'end_time' : "2011-07-04",  #"2018-01-01",  # end time of simulation [yyyy-mm-dd]
        'forc_filename' : "Lettosuo_forcing_2010_2018.csv",  # forcing data file*
        'variables': [['forcing_air_temperature', 'above canopy air temperature [degC]', ('date', 'simulation')],
                      ['forcing_precipitation', 'precipitation [m s-1]', ('date', 'simulation')],
                      ['forcing_h2o','H2O concentration [mol mol-1]', ('date', 'simulation')],
                      ['forcing_co2','CO2 concentration [ppm]', ('date', 'simulation')],
                      ['forcing_wind_speed','wind speed [m s-1]', ('date', 'simulation')],
#                      ['canopy_h2o','H2O concentration [mol mol-1]', ('date', 'simulation', 'canopy')],
#                      ['canopy_co2','CO2 concentration [ppm]', ('date', 'simulation', 'canopy')],
                      ['canopy_wind_speed','canopy wind speed [m s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_PAR_sunlit','incident PAR sunlit leaves [umol m-2 s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_PAR_shaded','incident PAR shaded leaves [umol m-2 s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_lad','leaf area density [m3 m-2]', ('date', 'simulation', 'canopy')],
                      ['canopy_sunlit_fraction','fraction of sunlit leafs [-]', ('date', 'simulation', 'canopy')],
                      ['canopy_LAI','canopy LAI [m2 m-2]', ('date', 'simulation')],
                      ['canopy_phenostate','canopy phenological state [-]', ('date', 'simulation')],
                      ['canopy_interception', 'canopy interception [m s-1]', ('date', 'simulation')],
                      ['canopy_interception_storage', 'canopy interception storage [m]', ('date', 'simulation')],
                      ['canopy_evaporation', 'evaporation from interception storage [m s-1]', ('date', 'simulation')],
                      ['canopy_transpiration','transpiration [m s-1]', ('date', 'simulation')],
                      ['canopy_pt_transpiration', 'transpiration [m s-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_pt_An', 'gross primary production [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_pt_Rd', 'dark respiration [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_throughfall', 'throughfall to moss or snow [m s-1]', ('date', 'simulation')],
                      ['canopy_potential_infiltration', 'potential infiltration to soil [m s-1]', ('date', 'simulation')],
                      ['canopy_snow_water_equivalent', 'snow water equivalent [m]', ('date', 'simulation')],
                      ['canopy_moss_evaporation', 'evaporation from moss layer [m s-1]', ('date', 'simulation')],
                      ['canopy_LE', 'latent heat flux [W m-2]', ('date', 'simulation')],
                      ['canopy_NEE', 'net ecosystem exchage [umol m-2 s-1]', ('date', 'simulation')],
                      ['canopy_GPP', 'ecosystem gross primary production [umol m-2 s-1]', ('date', 'simulation')],
                      ['canopy_Reco', 'ecosystem respiration [umol m-2 s-1]', ('date', 'simulation')],
#                      ['soil_water_potential','soil water potential [m]', ('date', 'simulation', 'soil')],
                      ['soil_pond_storage', 'pond storage [m]', ('date', 'simulation')],
                      ['soil_ground_water_level', 'ground water level [m]', ('date', 'simulation')],
                      ['soil_infiltration', 'infiltration [m s-1]', ('date', 'simulation')],
                      ['soil_surface_runoff', 'surface runoff [m s-1]', ('date', 'simulation')],
                      ['soil_evaporation', 'evaporation from soil surface [m s-1]', ('date', 'simulation')],
                      ['soil_subsurface_drainage', 'subsurface drainage [m s-1]', ('date', 'simulation')],
                      ['soil_total_runoff', 'total runoff [m s-1]', ('date', 'simulation')],
                      ['canopy_MBE1', 'interception model mass balance error [m]', ('date', 'simulation')],
                      ['canopy_MBE2', 'snow model mass balance error [m]', ('date', 'simulation')],
                      ['canopy_MBE3', 'moss model mass balance error [m]', ('date', 'simulation')],
                      ['soil_MBE', 'soil mass balance error [m]', ('date', 'simulation')],
                      ['canopy_z', 'canopy model grid node elevations [m]', ('canopy')],
                      ['soil_z', 'soil model grid node elevations [m]', ('soil')],
                      ['canopy_Rnet','[W m-2]', ('date', 'simulation')],
                      ['canopy_Rnet_ground', '[W m-2]', ('date', 'simulation')],
                      ['canopy_U_ground', '[m s-1]', ('date', 'simulation')]
                      ]
        }


#  output_folder: "results"
#  output_format: "netcdf"
