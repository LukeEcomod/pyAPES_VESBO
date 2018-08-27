# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:40:57 2018
@author: L1656
"""
import pandas as pd
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from tools.plotting import plot_timeseries_xr,plot_results, plot_fluxes, plot_columns, plot_pt_results, plot_lad_profiles, plot_timeseries_df
from tools.iotools import read_results, read_forcing, save_df_to_csv
from tools.dataprocessing_scripts import read_lettosuo_data
import seaborn as sns
pal = sns.color_palette("hls", 5)


results = read_results(['results/201808271237_CCFPeat_results.nc',
                        'results/201808271241_CCFPeat_results.nc',
                        'results/201808271244_CCFPeat_results.nc',
                        'results/201808271317_CCFPeat_results.nc'])
#results = read_results('results/201808271237_CCFPeat_results.nc')
#results = read_results('results/201808231333_CCFPeat_results.nc')
results = read_results(outputfile)
plot_results(results)
gwl_meas = read_forcing("lettosuo_WTD_pred.csv", cols='all')
plt.figure()
plot_timeseries_df(gwl_meas, ['part','clear','ctrl'],colors=[pal[2],pal[3],pal[0]],xticks=True, limits=False)
results['soil_pond_storage'].plot()

results['canopy_Rabs'].mean(dim='date').plot()
results['canopy_LWleaf'].mean(dim='date').plot()
results['canopy_Tleaf'].mean(dim='date').plot()
results['canopy_T'].mean(dim='date').plot()
results['canopy_Tleaf_wet'].mean(dim='date').plot()
results['canopy_lad'].mean(dim='date').plot()

results['canopy_IterWMA'].plot()
results['forcing_precipitation'].plot()

labels=['Ebal & no WMA', 'no Ebal & no WMA', 'Ebal & WMA', 'no Ebal & WMA']
plt.figure()
plot_timeseries_xr(results, 'canopy_IterWMA', labels=labels)
plt.figure()
plt.subplot(211)
plot_timeseries_xr(results, 'canopy_transpiration', labels=labels,unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
plt.subplot(212)
plot_timeseries_xr(results, 'canopy_transpiration', labels=labels, cum=True, unit_conversion={'unit':'mm', 'conversion':1e3})
plt.figure()
plt.subplot(211)
plot_timeseries_xr(results, 'canopy_evaporation', labels=labels, unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
plt.subplot(212)
plot_timeseries_xr(results, 'canopy_evaporation', labels=labels, cum=True, unit_conversion={'unit':'mm', 'conversion':1e3})

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
    (results[i]['canopy_sunlit_fraction'].sel(date=results[i]['date.month']==7).mean(dim='date') * results[i]['canopy_lad'].sel(date=results[i]['date.month']==7).mean(dim='date')).plot()
LAI=[np.sum(results[i]['canopy_lad'].sel(date=results[i]['date.month']==7).mean(dim='date').values*0.3) for i in range(2)]
LAI_sunlit=[np.sum(results[i]['canopy_sunlit_fraction'].sel(date=results[i]['date.month']==7).mean(dim='date').values * results[i]['canopy_lad'].sel(date=results[i]['date.month']==7).mean(dim='date').values*0.3) for i in range(2)]
plt.figure()
for i in range(2):
    results[i]['canopy_lad'].sel(date=results[i]['date.month']==7).mean(dim='date').plot()


plot_results(results)
plot_fluxes(results)
plot_pt_results(results,'canopy_pt_Tleaf')
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

""" WTD """
from tools.plotting import plot_timeseries_df
import seaborn as sns
pal = sns.color_palette("hls", 5)

gwl_meas = read_forcing("lettosuo_WTD.csv", cols='all')
#
#plot_columns(gwl_meas[['ctrl4m','ctrl8m','ctrl12m','ctrl22.5m']])
#plot_columns(gwl_meas[['part4m','part8m','part12m','part22.5m']])
#plot_columns(gwl_meas[['clear4m','clear8m','clear12m','clear22.5m']])
#plot_columns(gwl_meas[['WT_E','WT_N','WT_S','WT_W']])
gwl_meas['ctrl'] = np.nanmean([gwl_meas['ctrl8m'].values, gwl_meas['ctrl22.5m'].values],axis=0)
gwl_meas['clear4m'][gwl_meas.index < '11-01-2015']=np.nan
gwl_meas['clear12m'][gwl_meas.index < '11-01-2015']=np.nan

gwl_calib = gwl_meas[(gwl_meas.index <= '03-15-2016')]
#plot_columns(gwl_calib[['ctrl4m','ctrl12m', 'ctrl']],slope=1.0)

#plot_columns(gwl_calib[['WT_E','WT_N','WT_S','WT_W', 'ctrl']],slope=1.0)
gwl_meas['WT_Ec'] = gwl_meas['WT_E'] - 7.57
gwl_meas['WT_Nc'] = gwl_meas['WT_N'] -7.21
gwl_meas['WT_Sc'] = gwl_meas['WT_S'] - 16.67
gwl_meas['WT_Wc'] = gwl_meas['WT_W'] - 13.10
#plt.figure()
#plot_timeseries_df(gwl_meas, ['WT_Ec','WT_Sc','WT_Wc','WT_Nc','ctrl'],colors=pal,xticks=False, limits=False)
#
#plot_columns(gwl_calib[['part4m','part8m','part12m','part22.5m', 'ctrl']],slope=1.0)
gwl_meas['part12mc'] = gwl_meas['part12m'] -1.46
gwl_meas['part22.5mc'] = gwl_meas['part22.5m'] - 9.97
#plt.figure()
#plot_timeseries_df(gwl_meas, ['part12mc','part22.5mc','ctrl'],colors=pal,xticks=False, limits=False)
#
#plot_columns(gwl_calib[['clear4m','clear8m','clear12m','clear22.5m', 'ctrl']],slope=1.0)
gwl_meas['clear4mc'] = gwl_meas['clear4m'] - 2.32
gwl_meas['clear12mc'] = gwl_meas['clear12m'] - 9.39
#plt.figure()
#plot_timeseries_df(gwl_meas, ['clear4mc','clear12mc','ctrl'],colors=pal,xticks=False, limits=False)

plt.figure()
plot_timeseries_df(gwl_meas, ['WT_Ec','WT_Sc','WT_Wc'],colors=[pal[2]],xticks=False, limits=False)
plot_timeseries_df(gwl_meas, ['part12mc','part22.5mc'],colors=[pal[1]],xticks=False, limits=False)
plot_timeseries_df(gwl_meas, ['clear4mc','clear12mc'],colors=[pal[3]],xticks=False, limits=False)
plot_timeseries_df(gwl_meas, ['ctrl8m','ctrl22.5m'],colors=[pal[0]],xticks=True, limits=False)

gwl_meas['ctrl'] = np.nanmean([gwl_meas['ctrl8m'].values, gwl_meas['ctrl22.5m'].values],axis=0)
gwl_meas['part'] = np.nanmean([gwl_meas['part12mc'].values, gwl_meas['part22.5mc'].values,gwl_meas['WT_Ec'].values, gwl_meas['WT_Sc'].values,gwl_meas['WT_Wc'].values],axis=0)
gwl_meas['part1'] = np.nanmean([gwl_meas['part12mc'].values, gwl_meas['part22.5mc'].values],axis=0)
gwl_meas['part2'] = np.nanmean([gwl_meas['WT_Ec'].values, gwl_meas['WT_Sc'].values,gwl_meas['WT_Wc'].values],axis=0)
gwl_meas['clear'] = np.nanmean([gwl_meas['clear4mc'].values, gwl_meas['clear12mc'].values],axis=0)

for col in gwl_meas:
    gwl_meas[col] = gwl_meas[col] / 100.0
plt.figure()
plot_timeseries_df(gwl_meas, ['part','clear','ctrl'],colors=[pal[2],pal[3],pal[0]],xticks=True, limits=False)
plt.figure()
plot_timeseries_df(gwl_meas, ['part2','part1','clear','ctrl'],colors=[pal[1],pal[2],pal[3],pal[0]],xticks=True, limits=False)

save_df_to_csv(gwl_meas[['WT_Ec','WT_Sc','WT_Wc','part12mc','part22.5mc','clear4mc','clear12mc','ctrl8m','ctrl22.5m','part','clear','ctrl']],
               'lettosuo_WTD_pred', readme=' - Check timezone!! \nWT_Ec, WT_Sc,...: selected WTD timeseries, all predicted to control site using period from 1-Nov-2015 to 15-Mar-2016\npart, clear, ctrl: mean values for treatments')