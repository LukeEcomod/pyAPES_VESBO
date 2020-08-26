# -*- coding: utf-8 -*-
"""
Script to read annual nc-files, pick variables and save as a ncf-file
Created on Mon Aug 17 09:57:38 2020

@author: 03081268
"""

import pickle
import xarray
import glob
#import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# utility functions
#from utils.functions import vpd_from_rh 
#from pyAPES_utilities.bigleaf_functions import light_use_efficiency, water_use_efficiency, \
#                            intrinsic_water_use_efficiency, canopy_conductance, eq_evap, \
#                            saturation_vapor_pressure
from pyAPES_utilities.plotting import xarray_to_df
                         
scens = ['LCN', 'LC','LN','L','CN','C','N','ref']
sims = 8

#mres = []
#keys = ['date', 'NEE', 'GPP', 'Reco', 'ET', 'LUE', 'WUE', 'IWUE', 'Gs', 'ETeq', 'CiCa']

outfile = r'results/Scenarios/S1/scens_250820.pk'
    
fpath = r'results/Scenarios/S1/*.nc'
files = glob.glob(fpath)

data = []
for k in range(1,len(files)):
    print('reading ' + str(k))
    data.append(xarray.open_dataset(files[k]))
print('done')

##df = xarray.merge([xarray.open_dataset(f) for f in files]
#df = xarray.concat(data, 'date')
#df.to_netcdf(r'results/Scenarios/S1/S1all.nc', format='NETCDF4')
#print('ncf written')

#%%

out = []
for k in range(0, sims):
    dd = []
    m = 0
    for d in data:        
        print('converting: scen '  + str(k) + ' file ' + str(m))
        nc = xarray_to_df(d, [
                             'forcing_air_temperature',
                             'forcing_precipitation',
                             'forcing_pressure',
                             'forcing_h2o',
                             'forcing_co2',
                             'forcing_wind_speed',
                             'forcing_friction_velocity',
                             'forcing_par',
                             'forcing_nir',
                             'forcing_lw_in',
                             'canopy_LAI',
                             #'canopy_energy_closure',
                             #'canopy_fr_source',
                             'canopy_SWnet',
                             'canopy_LWnet',
                             'canopy_Rnet',
                             'canopy_interception',
                             #'canopy_interception_storage',
                             'canopy_evaporation',
                             #'canopy_condensation',
                             #'canopy_condensation_drip',
                             'canopy_throughfall',
                             #'canopy_water_closure',
                             'canopy_SH',
                             'canopy_LE',
                             'canopy_NEE',
                             'canopy_GPP',
                             'canopy_Reco',
                             'canopy_transpiration',
                             #'pt_total_gpp',
                             #'pt_total_dark_respiration',
                             #'pt_total_transpiration',
                             #'pt_total_stomatal_conductance_h2o',
                             #'pt_total_boundary_conductance_h2o',
                             #'soil_temperature',
                             #'soil_volumetric_water_content',
                             'ffloor_net_radiation',
                             'ffloor_sensible_heat',
                             'ffloor_latent_heat',
                             'ffloor_ground_heat',
                             #'ffloor_energy_closure',
                             #'ffloor_evaporation',
                             #'ffloor_soil_evaporation',
                             #'ffloor_throughfall',
                             'ffloor_interception',
                             #'ffloor_capillary_rise',
                             #'ffloor_pond_recharge',
                             #'ffloor_water_closure',
                             'ffloor_net_co2',
                             'ffloor_photosynthesis',
                             'ffloor_respiration',
                             'ffloor_soil_respiration',
                             'ffloor_surface_temperature',
                             #'ffloor_water_storage',
                             'ffloor_snow_water_equivalent',
                             #'ffloor_par_albedo',
                             #'ffloor_nir_albedo'],
                             ], sim_idx=k)
        dd.append(nc)
        del nc
        m += 1
    df = pd.concat(dd)
    out.append(df)
    del df
    
    file = open(outfile, 'wb')
    pickle.dump([out, scens], file)
    file.close()
    