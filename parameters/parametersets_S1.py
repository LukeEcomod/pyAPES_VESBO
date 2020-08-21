#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""
from parameters.parameter_tools import iterate_parameters
from tools.iotools import read_forcing
import pandas as pd
import numpy as np
import itertools
from copy import deepcopy as copy

# Modifications to override default parameters in parameters.SmearII
laifile=r'forcing/Hyytiala/Hyde_LAI_annual.csv'

laidata = pd.read_csv(laifile, sep=';', header='infer', index_col='year')
co2_ave = pd.DataFrame({'CO2': np.mean(laidata['CO2'])}, index=laidata.index)
co2_mean = np.mean(laidata['CO2'])

""" Simulation set 1: scenarios """

# default parameters for planttypes:
# Get parameters and forcing for SMEAR II -site

from parameters.SmearII import gpara, cpara, spara

default_params = {'general': gpara,
                  'canopy': cpara,
                  'soil': spara
                 }

# nominal parameters:
        
pts = cpara['planttypes']

lai_ref =  pts['pine']['LAImax'] + pts['spruce']['LAImax'] + pts['decid']['LAImax']
vmax_pine = pts['pine']['photop']['Vcmax']
jmax_pine = pts['pine']['photop']['Jmax']
rd_pine = pts['pine']['photop']['Rd']
del pts

co2_mean = 386.0 #ppm
rspruce = 0.25 #let's assume plant contributions are constant over time
rpine = 0.50
rdecid = 0.25

def get_parameter_list_S1(year=None, listout=False):
        
        ix = np.where(laidata.index == year)
#        laipine = laidata['LAIpine'].loc[year]
#        laispruce = laidata['LAIspruce'].loc[year]
#        laidecid = laidata['LAIdecid'].loc[year]
#        lai = laipine + laispruce + laidecid
#        fpine = laipine / (lai + EPS)
#        fspruce = laispruce / (lai + EPS)
#        fdec = laidecid / (lai + EPS)
        
        lai = laidata['LAItot'].loc[year]        
        co2 = laidata['CO2'].loc[year]

        fvmax = laidata['f_vcmax'].loc[year]
        
        print(ix, lai, co2, fvmax)
        
        comb = itertools.product([lai, lai_ref],[co2, co2_mean], [fvmax, 1.0])
        comb = list(comb)
        count = len(comb)
        
        # LAI per planttype
        LAIpine = tuple((comb[i][0] * rpine) for i in range(len(comb)))
        LAIdecid = tuple((comb[i][0] * rdecid) for i in range(len(comb)))
        LAIspruce = tuple((comb[i][0] * rspruce) for i in range(len(comb)))

        # Vcmax, Jmax and Rd in pine are altered by Nleaf; spruce and decid kept constant
        Vcmaxpine = tuple(vmax_pine * comb[i][2] for i in range(len(comb)))
        Jmaxpine = tuple(jmax_pine * comb[i][2] for i in range(len(comb)))
        Rdpine = tuple(rd_pine * comb[i][2] for i in range(len(comb)))

        # CO2
        CO2 = tuple(comb[i][1] for i in range(len(comb)))

        parameters = {
            'count': count,
            'scenario': 'S1',
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
        
        parameters['general'] = {
                    'start_time' : '%4d-04-01' %year,
                    'end_time' : '%4d-10-31' %year,
                    'forc_filename' : "Hyytiala/FIHy_forcing_1997-2019.dat"
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
            #print(param_list[k])
            forcing['CO2'] = param_list[k]['CO2']
            param_list[k]['forcing'] = forcing.copy()
            param_list[k]['nsim'] = k
        if listout:
            combnames = itertools.product(['L', 'Lr'], ['C', 'Cr'],['V', 'Vr'])
            combnames = list(combnames)
            return param_list, comb, combnames
        else:
            return param_list


# EOF