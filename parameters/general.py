# -*- coding: utf-8 -*-
"""
GENERAL PARAMETERS
"""

gpara = {
        'pyAPES_path': '/Users/ajkieloaho/Repositories/pyAPES/',
        'dt' : 1800.0,  # timestep in forcing data file [s]
        'start_time' : "2015-10-01",  #"2009-10-01",  # start time of simulation [yyyy-mm-dd]
        'end_time' : "2019-01-01",  #"2019-01-01",  # end time of simulation [yyyy-mm-dd]
        'forc_filename' : "Lettosuo_forcing_2010_2018.csv",  # forcing data file*
        'results_directory':'results/case_Lettosuo/',
        'variables': [['forcing_air_temperature', 'above canopy air temperature [degC]', ('date', 'simulation')],
                      ['forcing_precipitation', 'precipitation [m s-1]', ('date', 'simulation')],
                      ['forcing_pressure', 'ambient pressure [Pa]', ('date', 'simulation')],
                      ['forcing_h2o','H2O concentration [mol mol-1]', ('date', 'simulation')],
                      ['forcing_co2','CO2 concentration [ppm]', ('date', 'simulation')],
                      ['forcing_wind_speed','wind speed [m s-1]', ('date', 'simulation')],
                      ['forcing_friction_velocity','friction velocity [m s-1]', ('date', 'simulation')],
                      ['canopy_WMA_assumption','WMA assumed (1=True, 0=False)', ('date', 'simulation')],
#                      ['canopy_h2o','H2O concentration [mol mol-1]', ('date', 'simulation', 'canopy')],
#                      ['canopy_co2','CO2 concentration [ppm]', ('date', 'simulation', 'canopy')],
                      ['canopy_temperature','air temperature [degC]', ('date', 'simulation', 'canopy')],
#                      ['canopy_wind_speed','canopy wind speed [m s-1]', ('date', 'simulation', 'canopy')],
#                      ['canopy_friction_velocity','canopy friction velocity [m s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_lad','leaf area density [m3 m-2]', ('date', 'simulation', 'canopy')],
#                      ['canopy_sunlit_fraction','fraction of sunlit leafs [-]', ('date', 'simulation', 'canopy')],
                      ['canopy_LAI','canopy LAI [m2 m-2]', ('date', 'simulation')],
                      ['canopy_phenostate','canopy phenological state [-]', ('date', 'simulation')],
                      ['canopy_interception', 'canopy interception [m s-1]', ('date', 'simulation')],
                      ['canopy_interception_storage', 'canopy interception storage [m]', ('date', 'simulation')],
                      ['canopy_evaporation', 'evaporation from interception storage [m s-1]', ('date', 'simulation')],
                      ['canopy_condensation', 'condensation to canopy interception storage [m s-1]', ('date', 'simulation')],
                      ['canopy_condensation_drip', 'condensation to canopy that drips [m s-1]', ('date', 'simulation')],
                      ['canopy_transpiration','transpiration [m s-1]', ('date', 'simulation')],
                      ['canopy_pt_root_water_potential', 'root water potential [m]', ('date', 'simulation', 'planttype')],
                      ['canopy_pt_transpiration', 'transpiration [m s-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_pt_gpp', 'gross primary production [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_pt_respiration', 'dark respiration [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
#                      ['canopy_pt_stomatal_conductance_h2o', 'stomatal conductance for H2O [mol m-2 leaf s-1]', ('date', 'simulation', 'planttype')],
#                      ['canopy_pt_boundary_conductance_h2o', 'boundary layer conductance for H2O [mol m-2 leaf s-1]', ('date', 'simulation', 'planttype')],
#                      ['canopy_pt_leaf_internal_co2', 'leaf internal CO2 mixing ratio [mol mol-1]', ('date', 'simulation', 'planttype')],
#                      ['canopy_pt_leaf_surface_co2', 'leaf surface CO2 mixing ratio [mol mol-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_Tleaf', 'leaf temperature [degC]', ('date', 'simulation', 'canopy')],
#                      ['canopy_Tleaf_wet', 'wet leaf temperature [degC]', ('date', 'simulation', 'canopy')],
#                      ['canopy_Tleaf_sl', 'sunlit leaf temperature [degC]', ('date', 'simulation', 'canopy')],
#                      ['canopy_Tleaf_sh', 'shaded leaf temperature [degC]', ('date', 'simulation', 'canopy')],
#                      ['canopy_leaf_net_LW', 'net leaf longwave radiation [W m-2]', ('date', 'simulation', 'canopy')],
#                      ['canopy_leaf_SW_absorbed', 'leaf absorbed shortwave radiation [W m-2]', ('date', 'simulation', 'canopy')],
                      ['canopy_throughfall', 'throughfall to moss or snow [m s-1]', ('date', 'simulation')],
#                      ['canopy_evaporation_ml', 'evaporation from interception storage (condensation incl.) [m s-1]', ('date', 'simulation', 'canopy')],
#                      ['canopy_throughfall_ml', 'throughfall within canopy [m s-1]', ('date', 'simulation', 'canopy')],
#                      ['canopy_condensation_drip_ml', 'condensation drip within canopy [m s-1]', ('date', 'simulation', 'canopy')],
#                      ['canopy_co2_flux', 'co2 flux [umol m-2 s-1]', ('date', 'simulation', 'canopy')],
#                      ['canopy_latent_heat_flux', 'latent heat flux [W m-2]', ('date', 'simulation', 'canopy')],
#                      ['canopy_sensible_heat_flux', 'sensible heat flux [W m-2]', ('date', 'simulation', 'canopy')],
                      ['canopy_SH', 'sensible heat flux [W m-2]', ('date', 'simulation')],
                      ['canopy_LE', 'latent heat flux [W m-2]', ('date', 'simulation')],
                      ['canopy_SWnet', 'net shortwave radiation [W m-2]', ('date', 'simulation')],
                      ['canopy_LWnet', 'net longwave radiation [W m-2]', ('date', 'simulation')],
                      ['canopy_net_radiation', 'net radiation [W m-2]', ('date', 'simulation')],
                      ['canopy_NEE', 'net ecosystem exchage [umol m-2 s-1]', ('date', 'simulation')],
                      ['canopy_GPP', 'ecosystem gross primary production [umol m-2 s-1]', ('date', 'simulation')],
                      ['canopy_respiration', 'ecosystem respiration [umol m-2 s-1]', ('date', 'simulation')],
                      ['canopy_IterWMA', 'number of iterations [-]', ('date', 'simulation')],
                      ['soil_water_potential','soil water potential [m]', ('date', 'simulation', 'soil')],
                      ['soil_pond_storage', 'pond storage [m]', ('date', 'simulation')],
                      ['soil_ground_water_level', 'ground water level [m]', ('date', 'simulation')],
                      ['soil_infiltration', 'infiltration [m s-1]', ('date', 'simulation')],
                      ['soil_surface_runoff', 'surface runoff [m s-1]', ('date', 'simulation')],
                      ['soil_evaporation', 'evaporation from soil surface [m s-1]', ('date', 'simulation')],
                      ['soil_transpiration', 'transpiration from soil [m s-1]', ('date', 'simulation')],
                      ['soil_drainage', 'subsurface drainage [m s-1]', ('date', 'simulation')],
                      ['soil_temperature', 'soil temperature [degC]', ('date', 'simulation', 'soil')],
                      ['soil_volumetric_water_content', 'soil water content [m3/m3]', ('date', 'simulation', 'soil')],
                      ['soil_volumetric_ice_content', 'soil ice content [m3/m3]', ('date', 'simulation', 'soil')],
                      ['soil_heat_flux', 'soil heat flux [W m-2]', ('date', 'simulation', 'soil')],
#                      ['soil_thermal_conductivity', 'thermal conductivity [W m-1 K-1]', ('date', 'simulation', 'soil')],
                      ['canopy_water_closure', 'interception model mass balance error [m s-1]', ('date', 'simulation')],
                      ['canopy_energy_closure', 'energy closure in canopy [W m-2]', ('date', 'simulation')],
                      ['canopy_fr_source', 'Frsource in canopy [W m-2]', ('date', 'simulation')],
                      ['soil_water_closure', 'soil water balance error [m s-1]', ('date', 'simulation')],
                      ['soil_energy_closure', 'soil heat balance error [W m-2]', ('date', 'simulation')],
                      ['canopy_z', 'canopy model grid node elevations [m]', ('canopy')],
                      ['soil_z', 'soil model grid node elevations [m]', ('soil')],
                      ['ffloor_potential_infiltration', 'potential infiltration to soil [m s-1]', ('date', 'simulation')],
                      ['ffloor_snow_water_equivalent', 'snow water equivalent [m]', ('date', 'simulation')],
                      ['ffloor_ground_heat', 'ground heat flux (forest floor) [W m-2]', ('date', 'simulation')],
                      ['ffloor_sensible_heat', 'sensible heat flux (forest floor) [W m-2]', ('date', 'simulation')],
                      ['ffloor_latent_heat', 'latent heat flux (forest floor) [W m-2]', ('date', 'simulation')],
                      ['ffloor_snow_water_closure', "water balance error (snowcover) [m s-1]", ('date', 'simulation')],
                      ['ffloor_bryo_water_closure', "water balance error (bryophytes) [m s-1]", ('date', 'simulation')],
                      ['ffloor_bryo_energy_closure', "energy balance error (bryophytes) [W m-2]", ('date', 'simulation')],
                      ['ffloor_soil_energy_closure', "energy balance error (soil) [W m-2]", ('date', 'simulation')],
                      ['ffloor_bryo_carbon_pool', 'carbon pool (bryophyte) [kg C m-2]', ('date', 'simulation')],
                      ['ffloor_bryo_photosynthesis', 'photosynthesis rate (bryophyte) [umol m-2(ground) s-1]', ('date', 'simulation')],
                      ['ffloor_bryo_respiration', 'respiration rate (bryophyte) [umol m-2(ground) s-1]', ('date', 'simulation')],
                      ['ffloor_litter_respiration', 'respiration rate (litter) [umol m-2(ground) s-1]', ('date', 'simulation')],
                      ['ffloor_respiration', 'respiration rate [umol m-2(ground) s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation', 'evaporation (forest floor) [m s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation_bryo', 'evaporation (bryophytes) [m s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation_litter', 'evaporation (litter) [m s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation_soil', 'evaporation (soil) [m s-1]', ('date', 'simulation')],
                      ['ffloor_temperature', 'temperature (forest floor) [degC]', ('date', 'simulation')],
                      ['ffloor_litter_temperature', 'temperature (litter) [degC]', ('date', 'simulation')],
                      ['ffloor_bryo_temperature', 'temperature (bryophyte) [degC]', ('date', 'simulation')],
                      ['ffloor_soil_temperature', 'temperature (soil) [degC]', ('date', 'simulation')],
                      ['ffloor_bryo_water_storage', 'water storage (bryophytes) [kg m-2]', ('date', 'simulation')],
                      ['ffloor_litter_water_storage', 'water storage (litter) [kg m-2]', ('date', 'simulation')],
                      ['ffloor_capillar_rise', 'capillary rise to bryophyte layer [m s-1]', ('date', 'simulation')],
                      ]}

logging_configuration = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
                'default': {'format': '%(asctime)s %(levelname)s %(name)s %(message)s'},
                'model': {'format': '%(levelname)s %(name)s %(funcName)s %(message)s'},
                },
        'handlers': {
                'console': {
                        'class' : 'logging.StreamHandler',
                        'formatter': 'model',
                        'level': 'INFO'  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        },
                'file': {
                        'class': 'logging.FileHandler',
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'formatter': 'model',
                        'filename': 'pyAPES.log',
                        'mode': 'w',  # a == append, w == overwrite
                        },
                },
        'loggers': {
                'pyAPES': {
                        'handlers': ['file', 'console'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagate': True,
                        },
                'canopy':{
                        'handlers': ['file', 'console'],
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagate': True,
                        },
                'soil':{
                        'handlers': ['file', 'console'],
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagate': True,
                        },
                },
        }

parallel_logging_configuration = {
        'version': 1,
        'formatters': {
                'default': {
                    'class': 'logging.Formatter',
                    'format': '%(asctime)s %(levelname)s %(name)s %(message)s'},
                'model': {
                    'class': 'logging.Formatter',
                    'format': '%(process)d %(levelname)s %(name)s %(funcName)s %(message)s'},
        },
        'handlers': {
                'console': {
                        'class' : 'logging.StreamHandler',
                        'formatter': 'model',
                        'level': 'INFO'  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                },
                'pyAPES_file': {
                        'class': 'logging.FileHandler',
                        'level': 'WARNING',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'formatter': 'model',
                        'filename': 'pyAPES.log',
                        'mode': 'w',  # a == append, w == overwrite
                },
                'parallelAPES_file': {
                        'class': 'logging.FileHandler',
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'formatter': 'default',
                        'filename': 'parallelAPES.log',
                        'mode': 'w',  # a == append, w == overwrite
                },
        },
        'loggers': {
                'pyAPES': {
                        #'handlers': ['file'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagate': True,
                        },
                'canopy':{
                        #'handlers': ['file'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagete': True,
                        },
                'soil':{
                        #'handlers': ['file'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagete': True,
                },
        },
        'root': {
                'level': 'INFO',
                'handlers': ['console', 'parallelAPES_file']
        }
    }

#  output_folder: "results"
#  output_format: "netcdf"
