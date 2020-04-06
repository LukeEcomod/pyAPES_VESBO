# -*- coding: utf-8 -*-
"""
Created on Tue Oct 09 16:31:25 2018

@author: L1656
"""

from pyAPES import driver
from parameters.parameter_tools import get_parameter_list
from tools.iotools import read_results
from matplotlib import pyplot as plt
from pyAPES_utilities.plotting import plot_timeseries_xr, plot_timeseries_df, plot_fluxes
from tools.iotools import read_forcing

# Get parameters and forcing
# NOTE: start_time and end_time have to be same as in paramters.general.gpara
# this is because of netCDF that is initilized based on dates in general parameters
# it raises error msgs if it is not so. 
degero = get_parameter_list('Degero')

# %% Run model
outputfile=driver(parameters=degero, create_ncf=True)

#%%
results = read_results(outputfile)

Data = read_forcing("Degero_EC_2014_2016_ICOS.csv", cols='all',
                    start_time=results.date[0].values, end_time=results.date[-1].values)

# groud heat flux at 5 cm
results['ground_heat_flux'] = results['soil_heat_flux'].isel(soil=1).copy()

plot_fluxes(results, Data, norain=True,
            res_var=['canopy_Rnet','canopy_LWnet','canopy_SWnet', 'canopy_SH','canopy_LE','ground_heat_flux','canopy_NEE'],
            Data_var=['Rnet','LWnet','SWnet','SH','LE','G_5cm','NEE'],fmonth=5, lmonth=9)

# soil temperatures
results['soil_temperature_5cm'] = results['soil_temperature'].isel(soil=4).copy()
results['soil_temperature_10cm'] = results['soil_temperature'].isel(soil=9).copy()
results['soil_temperature_15cm'] = results['soil_temperature'].isel(soil=12).copy()
results['soil_temperature_30cm'] = results['soil_temperature'].isel(soil=19).copy()
results['soil_temperature_50cm'] = results['soil_temperature'].isel(soil=24).copy()

plot_fluxes(results, Data, norain=False,
            res_var=['soil_temperature_5cm','soil_temperature_10cm','soil_temperature_15cm','soil_temperature_30cm','soil_temperature_50cm'],
            Data_var=['Tsoil_5cm','Tsoil_10cm','Tsoil_15cm', 'Tsoil_30cm','Tsoil_50cm'],fmonth=5, lmonth=9)

plt.figure()
plot_timeseries_df(Data, 'GWL', colors=['k'])
plot_timeseries_xr(results, 'soil_ground_water_level')

plt.figure()
plot_timeseries_xr(results, 'ffloor_water_storage')
