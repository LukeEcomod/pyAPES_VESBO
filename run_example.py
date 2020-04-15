# -*- coding: utf-8 -*-
"""
Created on Tue Oct 09 16:31:25 2018

@author: Samuli Launiainen
"""

# to show figures in qt (pop-up's!)
# write in console: 
# import matplotlib.pyplot
# %matplotlib qt

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pyAPES import driver
#from parameters.parameter_tools import get_parameter_list
from tools.iotools import read_results
from tools.iotools import read_forcing, read_data
#from pyAPES_utilities.plotting import plot_timeseries_xr, plot_timeseries_df, plot_fluxes


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

# %% Run model, results are saved into netcdf-file and log-file pyAPES.log. driver returns filepath to results

outputfile, Model = driver(parameters=params, create_ncf=True)

#%% from now on, we just play with results and SMEAR II -data; this is all external to the model

# read results from NetCDF-file to xarray-dataset: xarray documentation here:
#http://xarray.pydata.org/en/stable/index.html

results = read_results(outputfile)

# import fluxdata and meteorological datafiles into pd.dataframes: pandas documentation here: 
# https://pandas.pydata.org/pandas-docs/stable/index.html

flxdata = read_data("forcing/Hyytiala/FIHy_flx_2005-2010.dat", sep=';',
                       start_time=results.date[0].values, end_time=results.date[-1].values)
metdata = read_data("forcing/Hyytiala/FIHy_met_2005-2010.dat", sep=';',
                       start_time=results.date[0].values, end_time=results.date[-1].values)

#%%
# --- prints content of results -dataset
results
# variables and their dimensions defined in parameters.outputs; brief explanation here:
# forcing_ : model forcing variables
# canopy_ : canopy-model results and integrated ecosystem level fluxes
# pt_ : planttype -specific results
# soil_ : soil-model results
# ffloor_ : forestfloor-model results
# gt_ : groundtype -specific results (part of ffloor_)

#it has following dimesions
results.dims

# date = time
# simulation (simulation nr, if only one, this is 0)
# canopy (canopy layers, 0 at ground),
# planttype (planttypes),
# soil (nr soil layers, 0 at top)
# groundtype (groundtypes at forestfloor)

# depending on variable, the number of dimensions vary.

# --- documentation of variable can be accessed e.g. as:
results['canopy_NEE'].attrs
results['canopy_NEE'].dims

results['pt_leaf_temperature'].attrs
results['pt_leaf_temperature'].dims

#%% let's now plot some variables
# indexes for dimensions:
sim = 0 # simulation 0
# python indices start from 0, -1 refers to last element in array, : means 'all'

# -- properties
t = results.date # time
zc = results.canopy_z # height above ground [m]
zs = results.soil_z # depth in soil is shown negative [m]

# now plot some results using matplolib; Kersti has some more advanced codes in pyAPES_utilities.plotting

# --- canopy structure and planttypes

# leaf-area density profile at end of simulation: these are not realistic for SMEAR II yet.
plt.figure('LAD')
var = 'canopy_lad'
plt.plot(results[var][-1, sim, :], zc); plt.xlabel(results[var].attrs['units'])
plt.ylabel(zc.attrs['units'])
plt.title('canopy LAI = ' + (str(results['canopy_LAI'][-1, sim,].values)))

#%% --- micromet.Micromet -submodel solves momentum exchange and scalar profiles in canopy air-space

# matric of heights for fast plotting
zcm = np.ones(np.shape(results['canopy_wind_speed'][:, sim, :])) * zc.values

# plot ensemble flow profiles: compute time-averaged wind speed and ust at each height

U = np.mean(results['canopy_wind_speed'][:, sim, :], axis=0)
ust = np.mean(results['canopy_friction_velocity'][:, sim, :], axis=0)

plt.figure('flow')

plt.subplot(121);
plt.plot(results['canopy_wind_speed'][:, sim, :], zcm, 'k', alpha=0.1)
plt.plot(U, zc, 'r-'); 
plt.xlabel('wind speed [m/s]')
plt.ylabel(zc.attrs['units'])

plt.subplot(122);
plt.plot(results['canopy_friction_velocity'][:, sim, :], zcm, 'k', alpha=0.1)
plt.plot(ust, zc, 'r-'); 
plt.xlabel('mean ustar [m/s]')
plt.ylabel(zc.attrs['units'])

# plot ensemble concentration profiles as difference to top layer; separately for day and night
par = results['forcing_par'][:,sim] # use PAR as criteria for night / day

var = ['canopy_temperature', 'canopy_co2', 'canopy_h2o']

plt.figure('scalarprofiles')
n = 1
for v in var:
    plt.subplot(2,2,n)
    x = results[v][:, sim, :] - results[v][:, sim, -1] # s - s_ref
    xm = np.mean(x, axis=0)
    xmn = np.mean(x.loc[par < 20, :], axis=0) # night
    xmd = np.mean(x.loc[par > 200, :], axis=0) # day
    plt.plot(x, zcm, 'k', alpha=0.1)
    plt.plot(xmd, zc, 'r', label='day'); plt.plot(xmn, zc, 'b', label='night')
    plt.xlabel(results[v].attrs['units']); plt.ylabel(zc.attrs['units'])
    n += 1

plt.legend()

# plot timeseries of some micromet variables at above-canopy (forging) and in sub-canopy
zsub = 3.0 # sub-canopy height [m], next find gridnode closest to this
ixsub = np.where((zc - zsub) == min(abs(zc - zsub)))[0]
#print(ixsub, zc[ixsub])

plt.figure('micromet')
var = ['canopy_wind_speed', 'canopy_friction_velocity', 'canopy_temperature', 'canopy_co2', 'canopy_h2o']

n = 1
for v in var:
    plt.subplot(5,1,n)
    plt.plot(t, results[v][:, sim, -1], 'k-', label='ref')
    plt.plot(t, results[v][:, sim, ixsub], 'r-', label='sub')
    plt.ylabel(results[v].attrs['units'])
    n +=1

#%% --- canopy.radiation.Radiation -model solves short-wave and long-wave radiation:
    # it provides profiles of incident downward, and upward radiation components [W m-2 (ground)]
    # and radiation absorbed per unit leaf area  [W m-2 (leaf) at each layer, separately for sunlit and 
    # shaded faction of leaves

# daytime short-wave profiles
var = ['canopy_par_down', 'canopy_par_up', 'canopy_nir_down', 'canopy_nir_up']

plt.figure('shortwave')

for v in var:
    plt.subplot(2,2,1)
    x = results[v][:, sim, :]
    #xm = np.mean(x, axis=0)
    #xmn = np.mean(x.loc[par < 20, :], axis=0) # night
    xmd = np.mean(x.loc[par > 200, :], axis=0) # day
    #plt.plot(x, zcm, 'k', alpha=0.1)
    plt.plot(xmd, zc, '-', label=results[v].attrs['units'])
    plt.xlabel(results[v].attrs['units']); plt.ylabel(zc.attrs['units'])
plt.legend()

var = ['canopy_lw_down', 'canopy_lw_up']
for v in var:
    plt.subplot(2,2,2)
    x = results[v][:, sim, :]
    #xm = np.mean(x, axis=0)
    #xmn = np.mean(x.loc[par < 20, :], axis=0) # night
    xmd = np.mean(x.loc[par > 200, :], axis=0) # day
    #plt.plot(x, zcm, 'k', alpha=0.1)
    plt.plot(xmd, zc, '-', label=results[v].attrs['units'])
    plt.xlabel(results[v].attrs['units']); plt.ylabel(zc.attrs['units'])
plt.legend()

# ... plot and some timeseries above canopy

plt.subplot(2,2,3)
# append Rnet from data
plt.plot(t, metdata['Rnet'], 'ko', alpha=0.3)

plt.plot(t, results['canopy_Rnet'], 'k-', label='Rn')
plt.plot(t, results['canopy_SWnet'], 'r-', label='SWn')
plt.plot(t, results['canopy_LWnet'], 'b-', label='LWn')
plt.ylabel('W m-2')
plt.legend()

plt.subplot(2,2,4)
alb = 1.0 - results['canopy_SWnet'] / (results['forcing_par'] + results['forcing_nir'])
plt.plot(t, alb, 'k')
plt.ylabel('SW albedo [-]')

#%% -- now let's move to ecosystem-level fluxes and look how they form

#--- timeseries of NEE, GGP, RECO

#var = ['canopy_NEE', 'canopy_GPP', 'canopy_Reco']

plt.figure('co2 fluxes')

plt.subplot(311); 
plt.plot(t, flxdata['NEE'], 'ko', alpha=0.3)
plt.plot(t, results['canopy_NEE'], 'r-');
plt.ylabel(results['canopy_NEE'].attrs['units'])

plt.subplot(312); 
plt.plot(t, flxdata['GPP'], 'ko', alpha=0.3)
plt.plot(t, results['canopy_GPP'], 'r-');
plt.ylabel(results['canopy_GPP'].attrs['units'])

plt.subplot(313); 
plt.plot(t, flxdata['Reco'], 'ko', alpha=0.3)
plt.plot(t, results['canopy_Reco'], 'r-');
plt.ylabel(results['canopy_Reco'].attrs['units'])

plt.figure('energy fluxes')

plt.subplot(411); 
plt.plot(t, metdata['Rnet'], 'ko', alpha=0.3)
plt.plot(t, results['canopy_Rnet'], 'r-');
plt.ylabel(results['canopy_Rnet'].attrs['units'])

plt.subplot(412); 
plt.plot(t, flxdata['H'], 'ko', alpha=0.3)
plt.plot(t, results['canopy_SH'], 'r-');
plt.ylabel(results['canopy_SH'].attrs['units'])

plt.subplot(413); 
plt.plot(t, flxdata['LE'], 'ko', alpha=0.3)
plt.plot(t, results['canopy_LE'], 'r-');
plt.ylabel(results['canopy_LE'].attrs['units'])

plt.subplot(414); 
plt.plot(t, flxdata['Gflux'], 'ko', alpha=0.3)

# let's take modeled G at 10cm depth
zref = -0.10 # soil depth [m], next find gridnode closest to this
ixs = np.where((zs - zref) == min(abs(zs - zref)))[0]

G = results['soil_heat_flux'][:, sim, ixs]
plt.plot(t, G, 'r-');
plt.ylabel(results['soil_heat_flux'].attrs['units'])

#%% -- ecosystem fluxes are integrated exchange rates from ground to canopy top.
# Example for daytime NEE and H, LE and leaf-air temperature difference

plt.figure('canopy flux profiles daytime', figsize=(8,8))

var = ['canopy_co2_flux', 'canopy_sensible_heat_flux', 'canopy_latent_heat_flux']

n = 1
for v in var:
    plt.subplot(2,2,n)
    x = results[v][:, sim, :]
    xmd = np.mean(x.loc[par > 200, :], axis=0) # day
    plt.plot(xmd, zc, '-', label=x.attrs['units'])
    plt.xlabel(x.attrs['units']); plt.ylabel(zc.attrs['units'])
    #plt.legend()
    n += 1

# leaf-air temperature difference, average at canopy layers
plt.subplot(2,2,n)
x = results['canopy_Tleaf'][:, sim, :] - results['canopy_temperature'][:, sim, :]
xmd = np.mean(x.loc[par > 200, :], axis=0) # day
plt.plot(xmd, zc, '-', label='T_l - T_a [degC]')
plt.xlabel('T_{leaf} - T_{air} [degC]'); plt.ylabel(zc.attrs['units'])

#%% -- leaf gas-exchange and leaf energy balance is computed in canopy.planttype -submodel
# let's look how sunlit and shaded leaves differ in pines; let's plot ensemble profiles and look timeseries
# of photosynthesis, heat fluxes and leaf temperature
print(results['canopy_planttypes'])
ptype = np.where(results['canopy_planttypes'] == 'pine')[0][0] # at index 1 we have pine

plt.figure('planttype sunlit / shaded leaf fluxes', figsize=(8,12))

vsl = ['pt_net_co2_sunlit', 'pt_latent_heat_sunlit', 'pt_sensible_heat_sunlit',
       'pt_stomatal_conductance_h2o_sunlit', 'pt_leaf_temperature_sunlit', 'pt_boundary_conductance_h2o_sunlit']
vsh = ['pt_net_co2_shaded', 'pt_latent_heat_shaded', 'pt_sensible_heat_shaded', 'pt_stomatal_conductance_h2o_shaded',
       'pt_leaf_temperature_shaded', 'pt_boundary_conductance_h2o_shaded']

for k in range(6):
    plt.subplot(3,2,k+1)
    
    xsl = np.mean(results[vsl[k]][:, sim, ptype, :].loc[par > 200, :], axis=0) # day
    plt.plot(xsl, zc, '-', label='sunlit')
    xsh = np.mean(results[vsh[k]][:, sim, ptype, :].loc[par > 200, :], axis=0) # day
    plt.plot(xsh, zc, '--', label='shaded')    
    plt.xlabel(results[vsl[k]].attrs['units']); plt.ylabel(zc.attrs['units'])
plt.legend()

# timeseries

plt.figure('sunlit leaves timeseries at all heights', figsize=(12,8))

for k in range(6):
    plt.subplot(3,2,k+1)
    
    plt.plot(t, results[vsl[k]][:, sim, ptype, :], '-', label='sunlit') 
    plt.ylabel(results[vsl[k]].attrs['units'])

#%% --- canopy.forestfloor -submodel handles moss and litter layer CO2, water and energy exchange

var = ['ffloor_net_radiation', 'ffloor_sensible_heat', 'ffloor_latent_heat',
       'ffloor_ground_heat', 'ffloor_photosynthesis', 'ffloor_water_storage']

plt.figure('forestfloor fluxes', figsize=(12,8))

for k in range(6):
    plt.subplot(3,2,k+1)
    
    plt.plot(t, results[var[k]][:, sim], '-', label='sunlit') 
    plt.ylabel(results[var[k]].attrs['units'])
    
#%% --- and forestfloor could consist of different groundtypes...

# ...
    
#%% --- soil - submodule computes soil moisture and temperature

var = ['soil_temperature', 'soil_volumetric_water_content']
lyrs = [0, 1, 2, 3, 4] # five top layers 
plt.figure('soil')
k = 1
for v in var:
    plt.subplot(2,1,k)
    plt.plot(t, results[v][:,sim,lyrs])
    plt.ylabel(results[v].attrs['units'])
    k += 1
