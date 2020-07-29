# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 18:06:46 2020

@author: Samuli Launiainen

SCRIPT TO TEST RUNNING SCENARIO S2 FOR HYDE TIMESERIES PAPER
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pyAPES import driver
from parameters.parametersets_hydetrends import get_parameter_list
from tools.iotools import read_results
from tools.iotools import  read_data
from pyAPES_utilities.plotting import plot_fluxes

# Get parametersets

params = get_parameter_list('S2', years=[2005, 2005])
result_file = 'S2_2005.nc'

#%%
outputfile, Model = driver(parameters=params, create_ncf=True, result_file=result_file)

# read results from NetCDF-file to xarray-dataset: xarray documentation here:
# http://xarray.pydata.org/en/stable/index.html

results = read_results(outputfile)

# import fluxdata and meteorological datafiles into pd.dataframes: pandas documentation here:
# https://pandas.pydata.org/pandas-docs/stable/index.html

flxdata = read_data("forcing/Hyytiala/FIHy_flx_2005-2010.dat", sep=';',
                       start_time=results.date[0].values, end_time=results.date[-1].values)
metdata = read_data("forcing/Hyytiala/FIHy_met_2005-2010.dat", sep=';',
                       start_time=results.date[0].values, end_time=results.date[-1].values)

#%%
plot_fluxes(results, flxdata, norain=True,
            res_var=['canopy_Rnet','canopy_SH','canopy_LE',
                      'canopy_NEE','canopy_GPP','canopy_Reco'],
            Data_var=['Rnet','H','LE','NEE','GPP','Reco'],
            fmonth=5, lmonth=9, sim_idx=5)

plt.figure()
results['soil_temperature'].isel(soil=0).plot.line(x='date')

plt.figure()
results['soil_volumetric_water_content'].isel(soil=0).plot.line(x='date')