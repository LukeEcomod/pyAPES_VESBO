# -*- coding: utf-8 -*-
"""
GENERAL PARAMETERS
"""

gpara = {
        'dt' : 1800.0,  # timestep in forcing data file [s]
        'start_time' : "2010-06-01",  # start time of simulation [yyyy-mm-dd]
        'end_time' : "2010-06-11",  #"2018-01-01",  # end time of simulation [yyyy-mm-dd]
        'forc_filename' : "Lettosuo_forcing_2010_2018.csv",  # forcing data file*
        'variables': [['forcing_air_temperature', 'above canopy air temperature [degC]', ('date', 'simulation')],
                      ['forcing_precipitation', 'precipitation [m s-1]', ('date', 'simulation')],
                      ['forcing_h2o','H2O concentration [mol mol-1]', ('date', 'simulation')],
                      ['forcing_co2','CO2 concentration [ppm]', ('date', 'simulation')],
                      ['forcing_wind_speed','wind speed [m s-1]', ('date', 'simulation')],
                      ['canopy_WMA_assumption','WMA assumed (1=True, 0=Flase)', ('date', 'simulation')],
                      ['canopy_h2o','H2O concentration [mol mol-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_co2','CO2 concentration [ppm]', ('date', 'simulation', 'canopy')],
                      ['canopy_temperature','air temperature []degC]', ('date', 'simulation', 'canopy')],
                      ['canopy_wind_speed','canopy wind speed [m s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_lad','leaf area density [m3 m-2]', ('date', 'simulation', 'canopy')],
                      ['canopy_sunlit_fraction','fraction of sunlit leafs [-]', ('date', 'simulation', 'canopy')],
                      ['canopy_LAI','canopy LAI [m2 m-2]', ('date', 'simulation')],
                      ['canopy_phenostate','canopy phenological state [-]', ('date', 'simulation')],
                      ['canopy_interception', 'canopy interception [m s-1]', ('date', 'simulation')],
                      ['canopy_interception_storage', 'canopy interception storage [m]', ('date', 'simulation')],
                      ['canopy_evaporation', 'evaporation from interception storage [m s-1]', ('date', 'simulation')],
                      ['canopy_condensation', 'condensation to canopy [m s-1]', ('date', 'simulation')],
                      ['canopy_transpiration','transpiration [m s-1]', ('date', 'simulation')],
                      ['canopy_pt_transpiration', 'transpiration [m s-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_pt_gpp', 'gross primary production [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_pt_respiration', 'dark respiration [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
                      ['canopy_Tleaf', 'leaf temperature [degC]', ('date', 'simulation', 'canopy')],
                      ['canopy_Tsurf', 'soil surface temperature [degC]', ('date', 'simulation')],
                      ['canopy_Tleaf_wet', 'wet leaf temperature [degC]', ('date', 'simulation', 'canopy')],
                      ['canopy_Tleaf_sl', 'sunlit leaf temperature [degC]', ('date', 'simulation', 'canopy')],
                      ['canopy_Tleaf_sh', 'shaded leaf temperature [degC]', ('date', 'simulation', 'canopy')],
                      ['canopy_net_leaf_LW', 'net leaf longwave radiation [W m-2]', ('date', 'simulation', 'canopy')],
                      ['canopy_net_leaf_radiation', 'leaf net radiation [W m-2]', ('date', 'simulation', 'canopy')],
                      ['canopy_throughfall', 'throughfall to moss or snow [m s-1]', ('date', 'simulation')],
                      ['canopy_co2_flux', 'co2 flux [umol m-2 s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_latent_heat_flux', 'latent heat flux [W m-2]', ('date', 'simulation', 'canopy')],
                      ['canopy_sensible_heat_flux', 'sensible heat flux [W m-2]', ('date', 'simulation', 'canopy')],
                      ['canopy_LE', 'latent heat flux [W m-2]', ('date', 'simulation')],
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
                      ['soil_drainage', 'subsurface drainage [m s-1]', ('date', 'simulation')],
                      ['soil_temperature', 'soil temperature [degC]', ('date', 'simulation', 'soil')],
                      ['soil_volumetric_water_content', 'soil ice content [m3/m3]', ('date', 'simulation', 'soil')],
                      ['soil_volumetric_ice_content', 'soil ice content [m3/m3]', ('date', 'simulation', 'soil')],
                      ['soil_thermal_conductivity', 'thermal conductivity [W m-1 K-1]', ('date', 'simulation', 'soil')],
                      ['canopy_water_closure_canopy', 'interception model mass balance error [m]', ('date', 'simulation')],
                      ['canopy_energy_closure_canopy', 'energy closure in canopy [W m-2]', ('date', 'simulation')],
                      ['canopy_fr_source', 'Frsource in canopy [W m-2]', ('date', 'simulation')],
                      ['soil_water_closure', 'soil mass balance error [m]', ('date', 'simulation')],
                      ['soil_energy_closure', 'soil heat balance error [W m2]', ('date', 'simulation')],
                      ['canopy_z', 'canopy model grid node elevations [m]', ('canopy')],
                      ['soil_z', 'soil model grid node elevations [m]', ('soil')],
                      ['ffloor_potential_infiltration', 'potential infiltration to soil [m s-1]', ('date', 'simulation')],
                      ['ffloor_snow_water_equivalent', 'snow water equivalent [m]', ('date', 'simulation')],
#                      ['ffloor_baresoil_evaporation', 'evaporation from baresoil [m s-1]', ('date', 'simulation')],
                      ['ffloor_ground_heat_flux', 'ground heat flux [W m-2]', ('date', 'simulation')],
                      ['ffloor_sensible_heat_flux', 'sensible heat flux from forest floor [W m-2]', ('date', 'simulation')],
                      ['ffloor_water_closure', 'forest floor mass balance error [m]', ('date', 'simulation')],
                      ['ffloor_water_closure_snow', 'snow model mass balance error [m]', ('date', 'simulation')],
                      ['ffloor_bryo_water_storage', 'water storage of bryophyte layer [kg m-2]', ('date', 'simulation')],
                      ['ffloor_bryo_temperature', 'bryophyte layer temperature [degC]', ('date', 'simulation')],
                      ['ffloor_bryo_evaporation', 'evaporation from bryophyte layer [m s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation', 'evaporation from forest floor [m s-1]', ('date', 'simulation')],
                      ['ffloor_soil_evaporation', 'soil evaporation [m s-1]', ('date', 'simulation')],
                      ['ffloor_temperature', 'forest floor temperature [degC]', ('date', 'simulation')],
                      ['ffloor_capillar_rise', 'capillary rise to bryophyte layer [m s-1]', ('date', 'simulation')]
                      ]
        }


#  output_folder: "results"
#  output_format: "netcdf"
