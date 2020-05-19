# -*- coding: utf-8 -*-
"""
GENERAL PARAMETERS FOR RUNNING pyAPES
"""

output_variables = {'variables': [# variable name, description [units], (dimensions)

      # copy of forcing variables
      ['forcing_air_temperature', 'above canopy air temperature [degC]', ('date', 'simulation')],
      ['forcing_precipitation', 'precipitation [kg m-2 s-1]', ('date', 'simulation')],
      ['forcing_pressure', 'ambient pressure [Pa]', ('date', 'simulation')],
      ['forcing_h2o','H2O concentration [mol mol-1]', ('date', 'simulation')],
      ['forcing_co2','CO2 concentration [ppm]', ('date', 'simulation')],
      ['forcing_wind_speed','wind speed [m s-1]', ('date', 'simulation')],
      ['forcing_friction_velocity','friction velocity [m s-1]', ('date', 'simulation')],
      ['forcing_par','downward par (direct + diffuse) [W m-2]', ('date', 'simulation')],
      ['forcing_nir','downward nir (direct + diffuse) [W m-2]', ('date', 'simulation')],
      ['forcing_lw_in','downward lw [W m-2]', ('date', 'simulation')],

      # canopy state and model control statistics
      ['canopy_LAI','canopy LAI [m2 m-2]', ('date', 'simulation')],
      ['canopy_lad','leaf area density [m3 m-2]', ('date', 'simulation', 'canopy')],
      ['canopy_phenostate','canopy phenological state [-]', ('date', 'simulation')],
      ['canopy_z', 'canopy model grid node elevations [m]', ('canopy')],
      ['canopy_planttypes', 'canopy planttype names', ('planttype')],
      ['canopy_WMA_assumption','WMA assumed (1=True, 0=False)', ('date', 'simulation')],
      ['canopy_IterWMA', 'number of iterations [-]', ('date', 'simulation')],
      ['canopy_energy_closure', 'energy closure in canopy [W m-2]', ('date', 'simulation')],
      ['canopy_fr_source', 'Frsource in canopy [W m-2]', ('date', 'simulation')],  #error related to isothermal long-wave balance

      # micromet profiles and canopy-layer average leaf temperatures
      ['canopy_h2o','H2O concentration [mol mol-1]', ('date', 'simulation', 'canopy')],
      ['canopy_co2','CO2 concentration [ppm]', ('date', 'simulation', 'canopy')],
      ['canopy_temperature','air temperature [degC]', ('date', 'simulation', 'canopy')],
      ['canopy_wind_speed','canopy wind speed [m s-1]', ('date', 'simulation', 'canopy')],
      ['canopy_friction_velocity','canopy friction velocity [m s-1]', ('date', 'simulation', 'canopy')],
      ['canopy_Tleaf', 'leaf temperature [degC]', ('date', 'simulation', 'canopy')],
      ['canopy_Tleaf_wet', 'wet leaf temperature [degC]', ('date', 'simulation', 'canopy')],
#      ['canopy_Tleaf_sl', 'sunlit leaf temperature [degC]', ('date', 'simulation', 'canopy')],
#      ['canopy_Tleaf_sh', 'shaded leaf temperature [degC]', ('date', 'simulation', 'canopy')],

      # radiation
#      ['canopy_sunlit_fraction','fraction of sunlit leafs [-]', ('date', 'simulation', 'canopy')],
      ['canopy_SWnet', 'net shortwave radiation balance at canopy top [W m-2]', ('date', 'simulation')],
      ['canopy_LWnet', 'net longwave radiation balance at canopy top [W m-2]', ('date', 'simulation')],
      ['canopy_Rnet', 'net radiation balance at canopy top [W m-2]', ('date', 'simulation')],
      # leaf scale, per m-2 leaf
#      ['canopy_leaf_net_LW', 'net leaf longwave radiation [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_leaf_net_SW', 'net leaf shortwave radiation [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_par_absorbed_sunlit', 'absorbed PAR of sunlit leaves [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_par_absorbed_shaded', 'absorbed PAR of shaded leaves [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_nir_absorbed_sunlit', 'absorbed NIR of sunlit leaves [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_nir_absorbed_shaded', 'absorbed NIR of shaded leaves [W m-2]', ('date', 'simulation', 'canopy')],
      # vertical profiles, per m-2 ground
#      ['canopy_par_down', 'downward PAR [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_par_up', 'upward PAR [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_nir_down', 'downward NIR [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_nir_up', 'upward NIR [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_lw_down', 'downward LW [W m-2]', ('date', 'simulation', 'canopy')],
#      ['canopy_lw_up', 'upward LW [W m-2]', ('date', 'simulation', 'canopy')],

      # interception sub-model results
      ['canopy_interception', 'canopy interception [kg m-2 s-1]', ('date', 'simulation')],
      ['canopy_interception_storage', 'canopy interception storage [kg m-2]', ('date', 'simulation')],
      ['canopy_evaporation', 'evaporation from interception storage [kg m-2 s-1]', ('date', 'simulation')],
      ['canopy_condensation', 'condensation to canopy interception storage [kg m-2 s-1]', ('date', 'simulation')],
      ['canopy_condensation_drip', 'condensation to canopy that drips [kg m-2 s-1]', ('date', 'simulation')],
      ['canopy_throughfall', 'throughfall to moss or snow [kg m-2 s-1]', ('date', 'simulation')],
      ['canopy_evaporation_ml', 'evaporation from interception storage, profile (condensation incl.) [kg m-2 s-1]', ('date', 'simulation', 'canopy')],
#      ['canopy_throughfall_ml', 'throughfall within canopy, profile [kg m-2 s-1]', ('date', 'simulation', 'canopy')],
#      ['canopy_condensation_drip_ml', 'condensation drip within canopy, profile [kg m-2 s-1]', ('date', 'simulation', 'canopy')],
      ['canopy_water_closure', 'interception model mass balance error [kg m-2 s-1]', ('date', 'simulation')],

      # ecosystem-level fluxes (at highest gridpoint, per m-2 ground)
      ['canopy_SH', 'sensible heat flux [W m-2]', ('date', 'simulation')],
      ['canopy_LE', 'latent heat flux [W m-2]', ('date', 'simulation')],
      ['canopy_NEE', 'net ecosystem exchange [umol m-2 s-1]', ('date', 'simulation')],
      ['canopy_GPP', 'ecosystem gross primary production [umol m-2 s-1]', ('date', 'simulation')],
      ['canopy_Reco', 'ecosystem respiration [umol m-2 s-1]', ('date', 'simulation')],
      ['canopy_transpiration', 'transpiration of all planttypes [m s-1]', ('date', 'simulation')],

      # flux profiles within canopy
      ['canopy_co2_flux', 'co2 flux [umol m-2 s-1]', ('date', 'simulation', 'canopy')],
      ['canopy_latent_heat_flux', 'latent heat flux [W m-2]', ('date', 'simulation', 'canopy')],
      ['canopy_sensible_heat_flux', 'sensible heat flux [W m-2]', ('date', 'simulation', 'canopy')],

      # root sink profile: Kersti - this should be taken out from Soil? Now len(rootsink) != len(Soil.z)
      #['canopy_root_sink', 'root water uptake profile [m s-1]', ('date', 'simulation', 'soil')],

      # planttype -specific outputs: lists of length 'planttype'
      ['pt_total_gpp', 'gross-primary productivity [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
#      ['pt_total_dark_respiration', 'dark (or leaf + wood?) respiration [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
      ['pt_total_transpiration', 'transpiration [kg m-2 s-1]', ('date', 'simulation', 'planttype')],
#      ['pt_total_stomatal_conductance_h2o', 'stomatal conductance for H2O [mol m-2 s-1]', ('date', 'simulation', 'planttype')],
#      ['pt_total_boundary_conductance_h2o', 'leaf boundary layer conductance for H2O [mol m-2 s-1]', ('date', 'simulation', 'planttype')],
#      ['pt_root_water_potential', 'root water potential [m?]', ('date', 'simulation', 'planttype')], # CHECK UNITS!!!

      # vertical profiles: lists of length 'planttype'; layers where lad == 0 are set to np.NaN
#      ['pt_leaf_temperature', 'leaf temperature mean [degC]', ('date', 'simulation', 'planttype', 'canopy')],
      ['pt_leaf_temperature_sunlit', 'leaf temperature, sunlit leaves [degC]', ('date', 'simulation', 'planttype', 'canopy')],
      ['pt_leaf_temperature_shaded', 'leaf temperature, shaded leaves [degC]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_net_co2_sunlit', 'net co2 uptake, sunlit leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_net_co2_shaded', 'net co2 uptake, shaded leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_dark_respiration_sunlit', 'dark respiration, sunlit leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_dark_respiration_shaded', 'dark respiration, shaded leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_transpiration_sunlit', 'transpiration, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_transpiration_shaded', 'transpiration, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_latent_heat_sunlit', 'latent heat flux, sunlit leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_latent_heat_shaded', 'latent heat flux, shaded leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_sensible_heat_sunlit', 'sensible heat flux, sunlit leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_sensible_heat_shaded', 'sensible heat flux, shaded leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_stomatal_conductance_h2o_sunlit', 'stomatal conductance for H2O, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_stomatal_conductance_h2o_shaded', 'stomatal conductance for H2O, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_boundary_conductance_h2o_sunlit', 'boundary-layer conductance for H2O, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_boundary_conductance_h2o_shaded', 'boundary-layer conductance for H2O, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_leaf_internal_co2_sunlit', 'leaf internal CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_leaf_internal_co2_shaded', 'leaf internal CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_leaf_surface_co2_sunlit', 'leaf surface CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
#      ['pt_leaf_surface_co2_shaded', 'leaf surface CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],

      # soil model state and fluxes
      ['soil_z', 'soil model grid node elevations [m]', ('soil')],
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
#      ['soil_thermal_conductivity', 'thermal conductivity [W m-1 K-1]', ('date', 'simulation', 'soil')],
      ['soil_water_closure', 'soil water balance error [m s-1]', ('date', 'simulation')],
      ['soil_energy_closure', 'soil heat balance error [W m-2]', ('date', 'simulation')],

      # forest floor outputs
      ['ffloor_net_radiation', 'net radiation (forest floor) [W m-2]', ('date', 'simulation')],
      ['ffloor_sensible_heat', 'sensible heat flux (forest floor) [W m-2]', ('date', 'simulation')],
      ['ffloor_latent_heat', 'latent heat flux (forest floor) [W m-2]', ('date', 'simulation')],
      ['ffloor_ground_heat', 'ground heat flux (forest floor) [W m-2]', ('date', 'simulation')],
      ['ffloor_energy_closure', "energy balance error (forest floor) [W m-2]", ('date', 'simulation')],

      ['ffloor_evaporation', 'evaporation (forest floor) [kg m-2 s-1]', ('date', 'simulation')],
      ['ffloor_soil_evaporation', 'evaporation (soil) [kg m-2 s-1]', ('date', 'simulation')],
      ['ffloor_throughfall', 'throughfall, potential infiltration (forest floor)) [kg m-2 s-1]', ('date', 'simulation')],
      ['ffloor_interception', 'interception rate (forest floor) [kg m-2 s-1]', ('date', 'simulation')],
      ['ffloor_capillary_rise', 'capillary rise (forest floor) [kg m-2 s-1]', ('date', 'simulation')],
      ['ffloor_pond_recharge', 'recharge from pond storage (forest floor) [kg m-2 s-1]', ('date', 'simulation')],
      ['ffloor_water_closure', "water balance error (forest floor) [kg m-2 s-1]", ('date', 'simulation')],

      ['ffloor_net_co2', 'net co2 flux (forest floor) [umol m-2(ground) s-1]', ('date', 'simulation')],
      ['ffloor_photosynthesis', 'photosynthesis rate (forest floor) [umol m-2(ground) s-1]', ('date', 'simulation')],
      ['ffloor_respiration', 'respiration rate (forest floor) [umol m-2(ground) s-1]', ('date', 'simulation')],
      ['ffloor_soil_respiration', 'soil respiration rate (below forest floor) [umol m-2(ground) s-1]', ('date', 'simulation')],

      ['ffloor_surface_temperature', 'temperature (forest floor) [deg C]', ('date', 'simulation')],
      ['ffloor_water_storage', 'water storage (forest floor) [kg m-2]', ('date', 'simulation')],
      ['ffloor_snow_water_equivalent', 'snow water equivalent (forest floor) [kg m-2]', ('date', 'simulation')],
      ['ffloor_par_albedo', 'PAR albedo (forest floor) [-]', ('date', 'simulation')],
      ['ffloor_nir_albedo', 'NIR albedo (forest floor) [-]', ('date', 'simulation')],
      ['ffloor_groundtypes', 'forestfloor groundtype names', ('groundtype')],

      # ground-type specific outputs
#       ['gt_net_radiation', 'net radiation [W m-2]', ('date', 'simulation', 'groundtype')],
      ['gt_sensible_heat', 'sensible heat flux [W m-2]', ('date', 'simulation', 'groundtype')],
#       ['gt_latent_heat', 'latent heat flux [W m-2]', ('date', 'simulation', 'groundtype')],
      ['gt_ground_heat', 'ground heat flux  [W m-2]', ('date', 'simulation', 'groundtype')],
#       ['gt_conducted_heat', 'conducted heat flux to moss [W m-2]', ('date', 'simulation', 'groundtype')],
#       ['gt_heat_advection', 'advected heat  [W m-2]', ('date', 'simulation', 'groundtype')],
#       ['gt_energy_closure', "energy balance error [W m-2]", ('date', 'simulation', 'groundtype')],

#       ['gt_evaporation', 'evaporation [kg m-2 s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_interception', 'interception rate [kg m-2 s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_capillary_rise', 'capillary rise [kg m-2 s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_pond_recharge', 'recharge from pond storage  [kg m-2 s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_throughfall', 'throughfall  [kg m-2 s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_water_closure', 'water balance error  [kg m-2 s-1]', ('date', 'simulation', 'groundtype')],

#       ['gt_net_co2', 'net co2 flux  [umol m-2 s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_photosynthesis', 'photosynthesis_rate [umol m-2 s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_respiration', 'respiration_rate [umol m-2 s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_internal_co2', 'internal co2 mixing ratio [ppm]]', ('date', 'simulation', 'groundtype')],
#       ['gt_conductance_co2', 'air-moss conductance for co2 [mol m-2 s-1]', ('date', 'simulation', 'groundtype')],

#       ['gt_temperature', 'temperature [degC]', ('date', 'simulation', 'groundtype')],
#       ['gt_surface_temperature', 'surface temperature [degC]', ('date', 'simulation', 'groundtype')],
#       ['gt_water_content', 'water content [g g-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_volumetric_water', 'volumetric water content [m3 m-3]', ('date', 'simulation', 'groundtype')],
#       ['gt_water_storage', 'water storage [kg m-2]', ('date', 'simulation', 'groundtype')],
#       ['gt_water_potential', 'water potential [m]', ('date', 'simulation', 'groundtype')],
#       ['gt_hydraulic_conductivity', 'hydraulic_conductivity [m s-1]', ('date', 'simulation', 'groundtype')],
#       ['gt_thermal_conductivity', 'thermal_conductivity [W m-1 K-1]', ('date', 'simulation', 'groundtype')],
      ]

}

# --- logger configuration. Antti / Kersti: add option to define logger output file name?
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

# for parallel simulations
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
