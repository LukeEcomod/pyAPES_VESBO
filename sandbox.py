# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:40:57 2018

@author: L1656
    """
import pandas as pd
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from tools.plotting import plotresults, plotxarray, plotxarray2, plotresultsMLM, plot_columns, plotcumulative
from tools.timeseries_tools import fill_gaps
import seaborn as sns
from tools.iotools import read_results, read_forcing, save_df_to_csv

plotresults(results.isel(simulation=0))

plotresultsMLM(results.isel(simulation=0))


results = read_results('201806041539_CCFPeat_results.nc')
results['canopy_GPP'].values = results['canopy_GPP'].values * 44.01 * 1e-3
results['canopy_Reco'].values = results['canopy_Reco'].values * 44.01 * 1e-3

dat = read_forcing("Lettosuo_data_2010_2018.csv", cols=['GPP','Reco','ET'])

dat.GPP=-dat.GPP

dat.GPP = dat.GPP * 1800 * 1e-6
dat.Reco = dat.Reco * 1800 * 1e-6


results['canopy_GPP'].values = -results['canopy_GPP'].values * 1800 * 44.01 * 1e-9
results['canopy_pt_An'].values = -results['canopy_pt_An'].values * 1800 * 44.01 * 1e-9
results['canopy_Reco'].values = results['canopy_Reco'].values * 1800 * 44.01 * 1e-9
results['canopy_pt_Rd'].values = results['canopy_pt_Rd'].values * 1800 * 44.01 * 1e-9

plt.figure(figsize=(9.5, 9.5))
pal = sns.color_palette("hls", 5)
if type(results) != list:
    results = [results]
for i in range(len(results)):
    result=results[i]
    plt.subplot(len(results),1,i+1)
    plotcumulative(result, 'canopy_Reco', colors=pal, xticks=True, cum=True,stack=True,label=['Soil'])
    plotcumulative([result.isel(planttype=i) for i in range(len(result.planttype))],
                'canopy_pt_Rd', colors=pal[1:], xticks=True,stack=True, cum=True,
                label=['Pine','Spruce','Birch','Shrubs'])
    plotcumulative(result, 'canopy_GPP', colors=[(0,0,0),(0,0,0)], xticks=True, cum=True,stack=True,label=['Moss'])
    plotcumulative([result.isel(planttype=i) for i in range(len(result.planttype))],
                'canopy_pt_An', colors=pal[1:], xticks=True,stack=True, cum=True)

plt.figure(figsize=(9.5, 9.5))
pal = sns.color_palette("hls", 5)
if type(results) != list:
    results = [results]
for i in range(len(results)):
    result=results[i]
    plt.subplot(len(results),1,i+1)
    plotcumulative([result.isel(planttype=i) for i in range(len(result.planttype))],
                'canopy_pt_transpiration', colors=pal[1:], xticks=True, m_to='mm',stack=True, cum=True,
                label=['pine','spruce','decid','shrubs'])
    plt.xlim('01-01-2016','01-01-2018')
    plt.ylim(0,0.22)
#    plotcumulative(result, 'canopy_Reco', colors=pal, xticks=True, cum=True,stack=True)
#    plotcumulative(result, 'canopy_GPP', colors=pal[2:], xticks=True, cum=True,stack=True)

#                   m_to='mm', label=', total')

yearly_cum = yearly_cumulative(dat, ['Reco'])
plt.plot(dat.index, yearly_cum[0], ':k')
yearly_cum = yearly_cumulative(dat, ['GPP'])
plt.plot(dat.index, yearly_cum[0],'k')


plt.tight_layout(rect=(0, 0, 0.8, 1))


output_control1 = driver(create_ncf=True, dbhfile="letto2009.txt")
output_control2 = driver(create_ncf=True, dbhfile="letto2014.txt")
output_partial = driver(create_ncf=True, dbhfile="letto2016_partial.txt")
output_clearcut = driver(create_ncf=True, dbhfile="letto2016_clearcut.txt")
output = driver(create_ncf=True, dbhfile="letto2014.txt", LAImax=[5, 0, 0, 0.7])
output = driver(create_ncf=True, LAI_sensitivity=True, dbhfile="letto2014.txt")

results = read_results(["results/201806050958_CCFPeat_results.nc",
                        "results/201806051022_CCFPeat_results.nc",
                        "results/201806051047_CCFPeat_results.nc"])

pal = sns.color_palette("hls", 4)
plt.figure(figsize=(9.5, 9.5))
plt.subplot(611); plotxarray2(results, 'soil_ground_water_level', colors=pal, xticks=False)
plt.subplot(612); plotxarray2(results, 'canopy_LAI', colors=pal, xticks=False)
plt.subplot(613); plotcumulative(results, 'canopy_transpiration', colors=pal, xticks=False, m_to='mm', cum=True)
plt.subplot(614); plotcumulative(results, 'canopy_evaporation', colors=pal, xticks=False, m_to='mm', cum=True)
plt.subplot(615); plotcumulative(results, 'canopy_moss_evaporation', colors=pal, xticks=False, m_to='mm', cum=True)
plt.subplot(616); plotcumulative(results, 'soil_total_runoff', colors=pal, xticks=True, m_to='mm', cum=True)
plt.tight_layout(rect=(0, 0, 0.8, 1))

plt.figure(figsize=(9.5, 9.5))
plt.subplot(511); plotxarray2(results, 'soil_ground_water_level', colors=pal, xticks=False)
plt.xlim('01-01-2016','01-01-2018')
plt.ylim(-0.8,0)
plt.subplot(512); plotcumulative(results, 'soil_total_runoff', colors=pal, xticks=True, m_to='mm', cum=True)
plt.xlim('01-01-2016','01-01-2018')
plt.subplot(513); plotcumulative(results, 'canopy_transpiration', colors=pal, xticks=False, m_to='mm', cum=True)
plt.xlim('01-01-2016','01-01-2018')
plt.subplot(514); plotcumulative(results, 'canopy_evaporation', colors=pal, xticks=False, m_to='mm', cum=True)
plt.xlim('01-01-2016','01-01-2018')
plt.subplot(515); plotcumulative(results, 'canopy_moss_evaporation', colors=pal, xticks=True, m_to='mm', cum=True)
plt.xlim('01-01-2016','01-01-2018')

plt.tight_layout(rect=(0, 0, 0.8, 1))

plt.figure(figsize=(9.5, 9.5))
plt.subplot(513); plotcumulative(results, 'canopy_GPP', colors=pal, xticks=False, m_to='mm', cum=True)
plt.xlim('01-01-2016','01-01-2018')
plt.subplot(514); plotcumulative(results, 'canopy_Reco', colors=pal, xticks=True, m_to='mm', cum=True)
plt.xlim('01-01-2016','01-01-2018')

plt.tight_layout(rect=(0, 0, 0.8, 1))

plt.figure(); 
plt.subplot(411); plotxarray2(results, 'canopy_GPP', colors=pal, xticks=False)
plt.subplot(412); plotcumulative(results, 'canopy_GPP', colors=pal, xticks=False)
plt.subplot(413); plotxarray2(results, 'canopy_LE', colors=pal, xticks=False)
plt.subplot(414); plotcumulative(results, 'canopy_LE', colors=pal, xticks=True)
plt.tight_layout(rect=(0, 0, 0.8, 1))

results=read_results(output)





variable = 'canopy_transpiration'
variable = 'canopy_GPP'
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

# RH
plot_columns(lettosuo_data,[44,57,63,81,92,100,112,118,124])

# Tair
plot_columns(lettosuo_data,[43,58,69,88,90,101,110,116,122])

# Prec
plot_columns(lettosuo_data,[89,107,108,109,64,115,121])
daily_data=lettosuo_data.resample('24H').sum()
plot_columns(daily_data,[85,100,101,102,60])

# wind
plot_columns(lettosuo_data,[123,127])
plot_columns(lettosuo_data,[91,111,117,123,19,23])
plot_columns(lettosuo_data,[19,23,7,24])
# ustar
plot_columns(lettosuo_data,[7,24])

# Global radiation
plot_columns(lettosuo_data,[95,97,52,53])

# Pressure
plot_columns(lettosuo_data,[59,93,119,125])
plot_columns(lettosuo_data,[59,-1])


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