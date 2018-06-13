# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:40:57 2018

@author: L1656
    """
import pandas as pd
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from tools.plotting import plot_results, plot_fluxes, plot_columns, plot_pt_results, plot_lad_profiles
from tools.iotools import read_results, read_forcing, save_df_to_csv
from tools.dataprocessing_scripts import read_lettosuo_data
import seaborn as sns

results = read_results([output_control2, output_partial])
results = read_results('results/201806130907_CCFPeat_results.nc')
results = read_results('results/201806131258_CCFPeat_results.nc')


plt.figure()
for i in range(2):
    results[i]['canopy_PAR_shaded'].sel(date=results[i]['date.month']==7).mean(dim='date').plot()
    results[i]['canopy_PAR_sunlit'].sel(date=results[i]['date.month']==7).mean(dim='date').plot()
plt.figure()
for i in range(2):
    results[i]['canopy_wind_speed'].sel(date=results[i]['date.month']==7).mean(dim='date').plot()
plt.figure()
for i in range(2):
    results[i]['canopy_sunlit_fraction'].sel(date=results[i]['date.month']==7).mean(dim='date').plot()
plt.figure()
for i in range(2):
    results[i]['canopy_lad'].sel(date=results[i]['date.month']==7).mean(dim='date').plot()


plot_results(results)
plot_fluxes(results)
plot_pt_results(results,'canopy_pt_transpiration')
plot_pt_results(results,'canopy_pt_An')
plot_lad_profiles("letto2014.txt", quantiles=[0.75, 1.0])

#output_control1 = driver(create_ncf=True, dbhfile="letto2009.txt")
output_control2 = driver(create_ncf=True, dbhfile="letto2014.txt")
output_partial = driver(create_ncf=True, dbhfile="letto2016_partial.txt")
#output_clearcut = driver(create_ncf=True, dbhfile="letto2016_clearcut.txt")
#output = driver(create_ncf=True, dbhfile="letto2014.txt", LAImax=[5, 0, 0, 0.7])
#output = driver(create_ncf=True, LAI_sensitivity=True, dbhfile="letto2014.txt")

variable = 'canopy_transpiration'

pal = sns.color_palette("hls", 13)

from parameters.sensitivity_sampling import LAIcombinations
plt.figure()
plt.subplot(211)
for i in range(len(LAIcombinations)):
    plotxarray2(results.isel(simulation=i), variable, colors=pal[i:], xticks=True, label=", stand LAI = " + str(sum(LAIcombinations[i])))
plt.subplot(212)
for i in range(len(LAIcombinations)):
    plotcumulative(results.isel(simulation=i), variable, colors=pal[i:], xticks=True, label=", stand LAI = " + str(sum(LAIcombinations[i])))
plt.tight_layout(rect=(0, 0, 0.8, 1))

"""-------- Colormesh plotting --------"""
from parameters.sensitivity_sampling import LAIcombinations
LAI=np.empty(len(LAIcombinations))
# pine, spruce, decid
string='psd'
species_f=np.empty([len(LAIcombinations),3])
species= ['('] * len(LAIcombinations)
for i in range(len(LAIcombinations)):
    LAI[i] = sum(LAIcombinations[i])
    species_f[i,:] = LAIcombinations[i]/LAI[i]
    for j in range(3):
        species[i] += str(int(species_f[i,j]*100))
        if j < 2:
            species[i] += ','
        else:
            species[i] += ')'
#        if species_f[i,j] >= 0.5:
#            species[i] += string[j].upper()
#        if species_f[i,j] == 0.25:
#            species[i] += string[j]
WTD = -results['soil_ground_water_level'].sel(date=results['date.season']=='JJA').groupby('date.year').mean(dim='date')
N_years = len(WTD['year'])
species_uniques = list(set(species))
LAI_uniques = list(set(LAI))

WTD_mesh = np.empty([len(species_uniques),len(LAI_uniques),N_years])
for i in range(len(species_uniques)):
    x = species_uniques[i]
    for j in range(len(LAI_uniques)):
        y = LAI_uniques[j]
        for k in range(len(WTD['simulation'])):
            if species[k] == x and LAI[k] == y:
                WTD_mesh[i,j,:]=WTD[:,k]


fig=plt.figure(figsize=(15,5))
for i in range(N_years):
    ax = plt.subplot(1,N_years+1,i+1)
    p = ax.pcolormesh(WTD_mesh[:,:,i],cmap=plt.get_cmap('coolwarm'),vmin=WTD.min(), vmax=WTD.max())
    ax.set_title(str(WTD['year'][i].values))
    yticks = np.arange(len(species_uniques))+0.5
    xticks = np.arange(len(LAI_uniques))+0.5
    ax.set_yticks(yticks)
    ax.set_xticks(xticks)
    ax.set_yticklabels([])
    ax.set_xticklabels(LAI_uniques)
#    ax.set_aspect('equal')
    if i == 0:
        ax.set_yticklabels(species_uniques)
        ax.set_ylabel('Species composition, share of total LAI (%)\n (pine, spruce, decidious)')
    if i == int(N_years/2):
        ax.set_xlabel('Total LAI of tree stand (m$^2$ m$^{-2}$)')

cax = plt.axes([0.94, 0.1, 0.015, 0.85])
cb = fig.colorbar(p, cax=cax)
cb.set_label('Mean water table depth during Jun-Aug (m)')
cb.ax.invert_yaxis()
fig.subplots_adjust(left=0.08, bottom=0.1, right=1.05, top=0.95, wspace=0.1)

fig.savefig('standstructure_WTD.png')

"""---------------------------------------------"""


"""plotting timeseries"""

lettosuo_data = read_lettosuo_data()

# RH
plot_columns(lettosuo_data,[44,57,63,81,92,100,112,118,124])

# Tair
plot_columns(lettosuo_data,[43,58,69,88,90,101,110,116,122])

# Prec
plot_columns(lettosuo_data,[89,107,108,109,64,115,121])
daily_data=lettosuo_data.resample('24H').sum()
plot_columns(daily_data,[85,100,101,102,60])

# wind
plot_columns(lettosuo_data,[91,111,117,123,19,23])
plot_columns(lettosuo_data,[19,23,7,24])
# ustar
plot_columns(lettosuo_data,[7,24])

# Global radiation
plot_columns(lettosuo_data,[95,97,52,53])

# Pressure
plot_columns(lettosuo_data,[59,93,119,125])
plot_columns(lettosuo_data,[59,-1])

# Soil moist
plot_columns(lettosuo_data,[47,48,49,50,65,66,67,68,84,85,86,87])

# NEE
plot_columns(lettosuo_data,[29,30,0])

# GPP
plot_columns(lettosuo_data,[31,32,2])

# Reco
plot_columns(lettosuo_data,[33,34,1])

# LE
plot_columns(lettosuo_data,[15,36])

""" lettosuo EC data """
lettosuo_data['L']=latent_heat(letto_data['Tair'].values)
lettosuo_data['ET']=lettosuo_data['Letto1_EC: LE (Wm-2)'] / lettosuo_data['L'] * 1800
lettosuo_data['ET'][lettosuo_data['ET'] < 0.0] = 0.0

letto_data = lettosuo_data[['Letto1_EC: LE (Wm-2)',
                            'ET',
                            'Letto1_EC: PaikNEE_N',
                            'Letto1_EC: GPP_N',
                            'Letto1_EC: Reco_N']]

letto_data=letto_data.rename(columns={'Letto1_EC: LE (Wm-2)':'LE',
                                      'Letto1_EC: PaikNEE_N':'NEE',
                                      'Letto1_EC: GPP_N':'GPP',
                                      'Letto1_EC: Reco_N':'Reco'})