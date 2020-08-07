#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""
from parameters.parameter_tools import iterate_parameters
from tools.iotools import read_forcing

import itertools
from copy import deepcopy as copy

# Modifications to override default parameters in parameters.SmearII

""" Simulation set 2: combinations """

# default parameters for planttypes:
# Get parameters and forcing for SMEAR II -site

from parameters.SmearII import gpara, cpara, spara

default_params = {'general': gpara,
                  'canopy': cpara,
                  'soil': spara
                 }
        
pts = cpara['planttypes']

co2_mean = 385.0 #ppm
rspruce = 0.25 #spruce contribution to total LAI stays roughly constant

lai_mean =  pts['pine']['LAImax'] + pts['spruce']['LAImax'] + pts['decid']['LAImax']

vmax_pine = pts['pine']['photop']['Vcmax']
jmax_pine = pts['pine']['photop']['Jmax']
rd_pine = pts['pine']['photop']['Rd']
del pts

# relative changes in parameters
rlai = [1.0, 0.78, 1.125]
rCO2 = [1.0, 0.93, 1.065]
rVcmax = [1.0, 0.95, 1.05]
rdec = [0.225]
#rdec = [0.18, 0.225, 0.27]

# make all combinations
comb = itertools.product(rlai, rCO2, rVcmax, rdec)
comb = list(comb)
# number of combinations
count = len(comb)
print(count)
# LAI per planttype
LAIpine = tuple((lai_mean * (1 - rspruce - comb[i][3]) * comb[i][0]) for i in range(len(comb)))
LAIdecid = tuple(lai_mean * comb[i][3] * comb[i][0] for i in range(len(comb)))
LAIspruce = tuple(lai_mean * rspruce * comb[i][0] for i in range(len(comb)))

# Vcmax, Jmax and Rd in pine are altered by Nleaf; spruce and decid kept constant
Vcmaxpine = tuple(vmax_pine * comb[i][2] for i in range(len(comb)))
Jmaxpine = tuple(jmax_pine * comb[i][2] for i in range(len(comb)))
Rdpine = tuple(rd_pine * comb[i][2] for i in range(len(comb)))

# CO2
CO2 = tuple(co2_mean * comb[i][1] for i in range(len(comb)))

def get_parameter_list(scenario, years=None):
    if scenario.upper()== 'S2':
        
        parameters = {
            'count': count,
            'scenario': 'S2',
            'CO2': CO2,
            
            'canopy': {
                'planttypes': {
                    'pine': {'LAImax': LAIpine,
                             'photop': {
                                 'Vcmax': Vcmaxpine,
                                 'Jmax': Jmaxpine,
                                 'Rd': Rdpine}
                             },
                    'spruce': {'LAImax': LAIspruce},
                    'decid': {'LAImax': LAIdecid}
                    }
                }
            }
        if years:
            parameters['general'] = {
                'start_time' : '%4d-06-01' %years[0],
                'end_time' : '%4d-06-02' %years[1],
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
            print(k)
            forcing['CO2'] = param_list[k]['CO2']
            param_list[k]['forcing'] = forcing.copy()
            
        return param_list
    
    else:
        raise ValueError("Unknown parameterset!")

# EOF