# -*- coding: utf-8 -*-
"""
Created on Thu May 14 13:35:04 2020

@author: 03110850
"""

from pyAPES_utilities.plotting import xarray_to_df, plot_timeseries_df
from matplotlib import pyplot as plt
prop_cycle = plt.rcParams['axes.prop_cycle']
default = prop_cycle.by_key()['color']
pal=default

def plot_daily(results, sim_idx=0):

    results['overstory_evaporation'] = results['canopy_evaporation_ml'][:,:,2:].sum(dim='canopy')
    results['understory_transpiration'] = results['pt_total_transpiration'][:,:,2]
    results['understory_evaporation'] = results['canopy_evaporation_ml'][:,:,:2].sum(dim='canopy')
    results['Tr_pine'] = results['pt_total_transpiration'][:,:,1]
    results['Tr_spruce'] = results['pt_total_transpiration'][:,:,3]
    results['Tr_birch'] = results['pt_total_transpiration'][:,:,0]
    results['gpp_pine'] = results['pt_total_gpp'][:,:,1]
    results['gpp_spruce'] = results['pt_total_gpp'][:,:,3]
    results['gpp_birch'] = results['pt_total_gpp'][:,:,0]
    results['gpp_understory'] = results['pt_total_gpp'][:,:,2]
    results['overstory_transpiration'] = results['Tr_pine'] + results['Tr_spruce'] + results['Tr_birch']
    results['overstory_gpp'] = results['gpp_pine'] + results['gpp_spruce'] + results['gpp_birch']

    # VPD
    from canopy.micromet import e_sat
    # vapor pressure
    esat, s = e_sat(results['forcing_air_temperature'].values)
    results['vpd'] = (esat - results['forcing_h2o'] * results['forcing_pressure']) / results['forcing_pressure']  # mol/mol

    results['evapotranspiration'] = (results['overstory_transpiration'] +
                         results['understory_transpiration'] +
                         results['overstory_evaporation'] +
                         results['understory_evaporation'] +
                         results['ffloor_evaporation'])
    results['evapotranspiration2'] = results['canopy_LE'] / 44100.0 * 18.015e-3

    variables = ['evapotranspiration2','evapotranspiration','canopy_transpiration', 'overstory_transpiration','canopy_evaporation','ffloor_evaporation',
                 'canopy_GPP','overstory_gpp','Tr_pine','Tr_spruce','Tr_birch','gpp_pine','gpp_spruce','gpp_birch','gpp_understory', 'vpd','ffloor_photosynthesis']

    Data = xarray_to_df(results, variables, sim_idx=sim_idx)

    for var in ['canopy_GPP','overstory_gpp','gpp_pine','gpp_spruce','gpp_birch','gpp_understory','ffloor_photosynthesis']:
        Data[var] *= 1e-6* 12.01e-3 * 1e3 * 3600 * 24  # umol m-2 s-1 - > gC m-2 d-1

    for var in ['evapotranspiration','evapotranspiration2','ffloor_evaporation','canopy_evaporation','canopy_transpiration', 'overstory_transpiration','Tr_pine','Tr_spruce','Tr_birch']:
        Data[var] *= 3600 * 24  # mm/s --> mm/d or kg m-2 d-1

    Data['canopy_transpiration'] *= 1000.0

    Data = Data.resample('D').mean()

    Data['iWUE_pine'] = (Data['gpp_pine'] * 1e3 / 12.01e-3) / (Data['Tr_pine'] / 18.015e-3 / Data['vpd'])
    Data['iWUE_spruce'] = (Data['gpp_spruce'] * 1e3 / 12.01e-3) / (Data['Tr_spruce'] / 18.015e-3 / Data['vpd'])
    Data['iWUE_birch'] = (Data['gpp_birch'] * 1e3 / 12.01e-3) / (Data['Tr_birch'] / 18.015e-3 / Data['vpd'])

    Data['gs_pine'] = (Data['Tr_pine'] / 18.015e-3 / Data['vpd'])
    Data['gs_spruce'] = (Data['Tr_spruce'] / 18.015e-3 / Data['vpd'])
    Data['gs_birch'] = (Data['Tr_birch'] / 18.015e-3 / Data['vpd'])

    variables = [['evapotranspiration','canopy_evaporation'],#,'ffloor_evaporation'],
                 ['canopy_transpiration', 'overstory_transpiration','Tr_pine','Tr_spruce','Tr_birch'],
                 ['canopy_GPP','overstory_gpp','gpp_pine','gpp_spruce','gpp_birch']]
                 # ['iWUE_pine','iWUE_spruce','iWUE_birch'],
                 # ['gs_pine','gs_spruce','gs_birch']]

    label=['mm d-1', 'mm d-1', 'gC m-2 d-1','','']

    plt.figure(figsize=(9,8))
    n=len(variables)
    for i in range(n):
        if i == 0:
            ax=plt.subplot(n,1,i+1)
        else:
            plt.subplot(n,1,i+1,sharex=ax)
        plot_timeseries_df(Data, variables[i], colors=pal, xticks=True,
                           labels=variables[i])
        plt.ylabel(label[i])
        plt.legend(loc="upper right", frameon=False, borderpad=0.0)

    plt.tight_layout()

    print(Data[['canopy_GPP','overstory_gpp','gpp_pine','gpp_spruce','gpp_birch','canopy_evaporation',
                'evapotranspiration','overstory_transpiration','Tr_pine','Tr_spruce','Tr_birch']].sum())

    return Data[['canopy_GPP','gpp_pine','gpp_spruce','gpp_birch','Tr_pine','Tr_spruce','Tr_birch']]