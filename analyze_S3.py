# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:13:08 2020

@author: 03081268
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pyAPES import driver
from parameters.parametersets_S3 import get_parameter_list_S3
from tools.iotools import read_results
from pyAPES_utilities.scenario_figs import summarize_scenarios
year = 2008

rfile = r'results/Scenarios/S3/S3_%4d.nc' % year

# Get parameterlist and forcing
params, p = get_parameter_list_S3(year, listout=True)

p = [list(i) for i in p] # list of lists
p = np.array(p)

p = {'LAI': p[:,0], 'Ca': p[:,1], 'Vmax': p[:,2], 'g1': p[:,3]}

#levs = {'LAI': np.unique(LAI), 'Ca': np.unique(Ca), 'Vmax': np.unique(Vmax), 'g1': np.unique(g1)}

# read results from NetCDF-file to xarray-dataset: xarray documentation here:
# http://xarray.pydata.org/en/stable/index.html

results = read_results(rfile)

#%%
resu = summarize_scenarios(results, p)

#%%
#cres, resu = scen_differences(results, pnames)

fig, ax = plt.subplots(3,3)
for k in np.unique(p['LAI']):
    ix0 = (resu['LAI']==k) & (resu['Vmax']==60) & (resu['g1']==2.0)
    X = resu.loc[ix0].copy()
    
    ca = X['Ca'].values
    
    ax[0,0].plot(ca, X.GPP/X.GPP.iloc[0], label='LAI='+ str(k)); ax[0,0].set_ylabel('GPP')
    ax[0,1].plot(ca, X.ET/X.ET.iloc[0]); ax[0,1].set_ylabel('ET')
    ax[0,2].plot(ca, X.LUE / X.LUE.iloc[0]); ax[0,2].set_ylabel('LUE')
    ax[1,0].plot(ca, X.WUE / X.WUE.iloc[0]); ax[1,0].set_ylabel('WUE')
    ax[1,1].plot(ca, X.Gs / X.Gs.iloc[0]); ax[1,1].set_ylabel('Gs')
    ax[1,2].plot(ca, X.CiCa); ax[1,2].set_ylabel('Ci/Ca')
    ax[2,0].plot(ca, X.alpha); ax[2,0].set_ylabel('ET/ETeq')
    ax[2,1].plot(ca, X.EF); ax[2,1].set_ylabel('EF')
    ax[2,2].plot(ca, X.tr / X.tr.iloc[0]); ax[2,2].set_ylabel('Tr')

for k in np.unique(p['LAI']):
    ix1 = (resu['LAI']==k) & (resu['Vmax']==40) & (resu['g1']==4.0)
    X = resu.loc[ix1].copy()
    
    ca = X['Ca'].values
    
    ax[0,0].plot(ca, X.GPP/X.GPP.iloc[0], '--', label='LAI='+ str(k)); ax[0,0].set_ylabel('GPP')
    ax[0,1].plot(ca, X.ET/X.ET.iloc[0], '--'); ax[0,1].set_ylabel('ET')
    ax[0,2].plot(ca, X.LUE / X.LUE.iloc[0], '--'); ax[0,2].set_ylabel('LUE')
    ax[1,0].plot(ca, X.WUE / X.WUE.iloc[0], '--'); ax[1,0].set_ylabel('WUE')
    ax[1,1].plot(ca, X.Gs / X.Gs.iloc[0], '--'); ax[1,1].set_ylabel('Gs')
    ax[1,2].plot(ca, X.CiCa, '--'); ax[1,2].set_ylabel('Ci/Ca')
    ax[2,0].plot(ca, X.alpha, '--'); ax[2,0].set_ylabel('ET/ETeq')
    ax[2,1].plot(ca, X.EF, '--'); ax[2,1].set_ylabel('EF')
    ax[2,2].plot(ca, X.tr / X.tr.iloc[0], '--'); ax[2,2].set_ylabel('Tr')
    
ax[0,0].legend()