# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 10:58:33 2020

@author: 03081268
"""

# import matplotlib.pyplot
# %matplotlib qt

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pyAPES import driver
from parameters.parameter_tools import get_parameter_list
from tools.iotools import read_results
from tools.iotools import read_forcing, read_data
from pyAPES_utilities.plotting import plot_fluxes

# Get parameters and forcing for SMEAR II -site

from parameters.SmearII import gpara, cpara, spara

forcing = read_forcing(
    forc_filename=gpara['forc_filename'],
    start_time=gpara['start_time'],
    end_time=gpara['end_time'],
    dt=gpara['dt']
)

#  wrap parameters and forcing in dictionary
params = {
    'general': gpara,
    'canopy': cpara,
    'soil': spara,
    'forcing': forcing
}

# to run multiple simulation with some variable parameters
params = get_parameter_list(params, 'test')

outputfile, Model = driver(parameters=params, create_ncf=True)

# read results from NetCDF-file to xarray-dataset: xarray documentation here:
# http://xarray.pydata.org/en/stable/index.html

results = read_results(outputfile)

# import fluxdata and meteorological datafiles into pd.dataframes: pandas documentation here:
# https://pandas.pydata.org/pandas-docs/stable/index.html

flxdata = read_data("forcing/Hyytiala/FIHy_flx_2005-2010.dat", sep=';',
                       start_time=results.date[0].values, end_time=results.date[-1].values)
metdata = read_data("forcing/Hyytiala/FIHy_met_2005-2010.dat", sep=';',
                       start_time=results.date[0].values, end_time=results.date[-1].values)