#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""
import numpy as np
from parameters.parameter_tools import iterate_parameters
from tools.iotools import read_forcing

import itertools
from copy import deepcopy as copy

# Modifications to override default parameters in parameters.SmearII

""" Simulation set 3: combinations """

# default parameters for planttypes

from parameters.generic_forest import gpara, cpara, spara

default_params = {'general': gpara,
                  'canopy': cpara,
                  'soil': spara
                 }
        
pts = cpara['planttypes']

co2_mean = 385.0 #ppm

rjv = pts['trees']['photop']['Jmax'] / pts['trees']['photop']['Vcmax']
rrv = pts['trees']['photop']['Rd'] / pts['trees']['photop']['Vcmax']
del pts

# parameter ranges

slai = ['LAI_r','LAI_l','LAI_h']
sCO2 = ['CO2_r','CO2_l','CO2_h']
sVcmax =['Vmax_r','Vmax_l','Vmax_h']
sdec = ['df_r']

lai = [2.0, 4.0, 6.0] #np.linspace(1.0, 7.0, 7)
co2 = co2_mean * np.linspace(1.0, 1.6, 10)
vcmax = [60.0] #, 80.0]
m = [2.0, 4.0]

# make all combinations
#combnames = itertools.product(slai, sCO2, sVcmax, sdec)
#combnames = list(combnames)
comb = itertools.product(lai, co2, vcmax, m) #, Vcmax, g1)
comb = list(comb)
# number of combinations
count = len(comb)
print(count)

# LAI per planttype
LAI = tuple(comb[i][0] for i in range(len(comb)))
CO2 = tuple(comb[i][1] for i in range(len(comb)))
Vcmax = tuple(comb[i][2] for i in range(len(comb)))
g1 = tuple(comb[i][3] for i in range(len(comb)))

Jmax = tuple(rjv * comb[i][2] for i in range(len(comb)))
Rd = tuple(rrv * comb[i][2] for i in range(len(comb)))

def get_parameter_list_S3(year, listout=False):
        
    parameters = {
        'count': count,
        'scenario': 'S2',
        'CO2': CO2,
        
        'canopy': {
            'planttypes': {
                'trees': {'LAImax': LAI,
                         'photop': {
                             'Vcmax': Vcmax,
                             'Jmax': Jmax,
                             'Rd': Rd}
                         },
                }
            }
        }
    parameters['general'] = {
                'start_time' : '%4d-07-01' %year,
                'end_time' : '%4d-07-15' %year,
                'forc_filename' : "Hyytiala/FIHy_forcing_2005-2010.dat"
        }

    # make parameter lists
    param_list = [iterate_parameters(parameters, copy(default_params), c) for c in range(count)]
    
    # read forcing, replace CO2
    start_time = param_list[0]['general']['start_time']
    end_time = param_list[0]['general']['end_time']
    dt = param_list[0]['general']['dt']
    forc_filename = param_list[0]['general']['forc_filename']

    
    forcing = read_forcing(
        forc_filename=forc_filename,
        start_time=start_time,
        end_time=end_time,
        dt=dt
    )
    
    for k in range(count):
        #print(k)
        forcing['CO2'] = param_list[k]['CO2']
        param_list[k]['forcing'] = forcing.copy()
        param_list[k]['nsim'] = k
    if listout:
        return param_list, comb
    else:
        return param_list


# EOF