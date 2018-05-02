# -*- coding: utf-8 -*-
"""
GENERAL PARAMETERS
"""

gpara = {
        'dt' : 1800.0,  # timestep in forcing data file [s]
        'start_time' : "2010-01-01",  # start time of simulation [yyyy-mm-dd]
        'end_time' : "2016-12-31",  # end time of simulation [yyyy-mm-dd]
        'forc_filename' : "Hyde_2010_2016.csv",  #"FMI_jokioinen.csv"  ## forcing data file*
        'variables': [['forcing_air_temperature', 'above canopy air temperature [degC]', ('date', 'simulation')],
                      ['forcing_precipitation', 'precipitation [m]', ('date', 'simulation')],
                      ['forcing_h2o','H2O concentration [mol mol-1]', ('date', 'simulation')],
                      ['forcing_co2','CO2 concentration [ppm]', ('date', 'simulation')],
                      ['forcing_wind_speed','wind speed [m s-1]', ('date', 'simulation')],
#                      ['canopy_h2o','H2O concentration [mol mol-1]', ('date', 'simulation', 'canopy')],
#                      ['canopy_co2','CO2 concentration [ppm]', ('date', 'simulation', 'canopy')],
#                      ['canopy_wind_speed','canopy wind speed [m s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_LAI','canopy LAI [m2 m-2]', ('date', 'simulation')],
                      ['canopy_phenostate','canopy phenological state [-]', ('date', 'simulation')],
                      ['canopy_interception', 'canopy interception [m]', ('date', 'simulation')],
                      ['canopy_interception_storage', 'canopy interception storage [m]', ('date', 'simulation')],
                      ['canopy_evaporation', 'evaporation from interception storage [m]', ('date', 'simulation')],
                      ['canopy_transpiration','transpiration [m]', ('date', 'simulation')],
                      ['canopy_throughfall', 'throughfall to moss or snow [m]', ('date', 'simulation')],
                      ['canopy_potential_infiltration', 'potential infiltration to soil [m]', ('date', 'simulation')],
                      ['canopy_snow_water_equivalent', 'snow water equivalent [m]', ('date', 'simulation')],
                      ['canopy_moss_evaporation', 'evaporation from moss layer [m]', ('date', 'simulation')],
                      ['canopy_LE', 'latent heat flux [W m-2]]', ('date', 'simulation')],
                      ['canopy_NEE', 'net ecosystem exchage [umol m-2 s-1]', ('date', 'simulation')],
                      ['canopy_GPP', 'ecosystem gross primary production [umol m-2 s-1]', ('date', 'simulation')],
#                      ['soil_water_potential','soil water potential [m]', ('date', 'simulation', 'soil')],
                      ['soil_pond_storage', 'pond storage [m]', ('date', 'simulation')],
                      ['soil_ground_water_level', 'ground water level [m]', ('date', 'simulation')],
                      ['soil_infiltration', 'infiltration [m]', ('date', 'simulation')],
                      ['soil_surface_runoff', 'surface runoff [m]', ('date', 'simulation')],
                      ['soil_evaporation', 'evaporation from soil surface [m]', ('date', 'simulation')],
                      ['soil_subsurface_drainage', 'subsurface drainage [m]', ('date', 'simulation')],
                      ['soil_total_runoff', 'total runoff [m]', ('date', 'simulation')],
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

#  * Forcing data file (comma delimited, .csv)
#    Headers:
#        yyyy,mo,dd,hh,mm - year, month, day, hour, minute
#        Prec - Precipitation [mm/dt]
#        Tair - Air temperature [degC]
#        U - Wind speed 10 min avg. [m/s]
#        RH - Relative humidity [%]
#        Rg - Global radiation [W/m2]
