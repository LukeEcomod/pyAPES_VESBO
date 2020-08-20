# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:13:08 2020

@author: 03081268
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pyAPES import driver
from parameters.parametersets_S2 import get_parameter_list_S2
from tools.iotools import read_results
from tools.iotools import read_forcing, read_data
from pyAPES_utilities.plotting import plot_fluxes

scen = 'S2'
year = 2008

rfile = 'results/Scenarios/' + scen + '/' + scen + '_%4d.nc' %year

# Get parameterlist and forcing
params, p, pnames = get_parameter_list_S2('S2', years=[year, year], listout=True)

p = [list(i) for i in p] # list of lists
p = np.array(p)

pnames = [list(i) for i in pnames] # list of lists

rLAI = p[:,0]
rVmax = p[:,2]
rCO2 = p[:,1]

sRef = 0
sL = np.where((rVmax==1) & (rCO2 == 1))
sC = np.where((rLAI == 1) & (rVmax == 1))
sV = np.where((rLAI==1) & (rCO2 == 1))
sLC = np.where(rVmax==1)
sLV = np.where(rCO2==1)
sCV = np.where(rLAI==1)
# read results from NetCDF-file to xarray-dataset: xarray documentation here:
# http://xarray.pydata.org/en/stable/index.html

results = read_results(rfile)

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
            fmonth=5, lmonth=9, sim_idx=0)

#%%
cres, resu = scen_differences(results, pnames)



