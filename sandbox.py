# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:40:57 2018
@author: L1656
"""

import pandas as pd
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from tools.plotting import plot_timeseries_xr,plot_results, plot_fluxes, plot_columns, plot_pt_results, plot_lad_profiles, plot_timeseries_df, plot_lad_profiles
from tools.iotools import read_results, read_forcing, save_df_to_csv
from tools.dataprocessing_scripts import read_lettosuo_data
import seaborn as sns
pal = sns.color_palette("hls", 5)

results = read_results(outputfile)
plot_results(results)
plot_fluxes(results)

plt.figure()
plot_timeseries_xr(results.isel(soil=0, canopy=0), ['forcing_air_temperature','canopy_temperature','canopy_Tsurf','soil_temperature'])

plt.figure()
plot_timeseries_xr(results.isel(soil=0), ['canopy_ground_heat_flux','canopy_soil_sensible_heat_flux'])

plt.figure()
plot_timeseries_xr(results.isel(soil=0), ['soil_heat_be'])


T_surf, Hw, Frw, Gw, Ep, LEw, closure = baresoil_energybalance(z_can=0.3, U=0.002, T=8, H2O=0.005, P=98630,
                       T_ave=8, soil_alb={'Par':0.05, 'Nir': 0.5}, soil_emi=0.98, zr=0.01,
                       T_soil=7, h_soil=-0.41, z_soil=-0.005, Kh=5e-7, Kt=0.05,
                       Par_gr=20, Nir_gr=20, LWn=-10, Ebal=True)

#results = read_results(['results/201808271237_CCFPeat_results.nc',
#                        'results/201808271241_CCFPeat_results.nc',
#                        'results/201808271244_CCFPeat_results.nc',
#                        'results/201808271317_CCFPeat_results.nc'])
#results = read_results('results/201809131815_CCFPeat_results.nc')


control=driver(create_ncf=True, dbhfile="letto2014.txt")
partial=driver(create_ncf=True, dbhfile="letto2016_partial.txt")
clearcut=driver(create_ncf=True, dbhfile="letto2016_clearcut.txt")

results = read_results([control, partial, clearcut])

variable = 'canopy_transpiration'
plt.figure()
plt.subplot(211)
plot_timeseries_xr(results, variable, unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
plt.subplot(212)
plot_timeseries_xr(results, variable, cum=True, unit_conversion={'unit':'mm', 'conversion':1e3})

dbhfile="letto2016_partial"
plot_lad_profiles(dbhfile + ".txt")
plt.xlim(0,0.52)
plt.ylim(0,25)
plt.xlabel('leaf area density [m$^2$m$^{-3}$]')
plt.gca().axes.spines['top'].set_visible(False)
plt.gca().axes.spines['right'].set_visible(False)
#plt.legend(frameon=False, borderpad=0.0, labels=['Pine','Spruce', 'Birch'])
#plt.savefig(dbhfile +'.png',dpi=500)

""" WATER BALANCE FIGURE """
pal = np.array(sns.color_palette("hls", 7))
pal[1] = pal[0]
pal=pal[1:]
pal[1] = (pal[1]**1.5)
plt.figure(figsize=(7,2))
yyyy=2017
ax = plt.subplot(141)
i=0
plot_timeseries_xr(results[i].sel(date=results[i]['date.year']==yyyy), ['canopy_moss_evaporation', 'canopy_transpiration', 'canopy_evaporation',
                   'soil_subsurface_drainage'],
                   colors=pal, cum=True, stack=True, xticks=False,
                   unit_conversion={'unit':'mm', 'conversion':1e3}, legend=False)
plt.xticks([])
plt.ylabel('')
plt.yticks(np.linspace(0.,800.,3))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax1=plt.subplot(142, sharey=ax)
i=1
plot_timeseries_xr(results[i].sel(date=results[i]['date.year']==yyyy), ['canopy_moss_evaporation', 'canopy_transpiration', 'canopy_evaporation',
                   'soil_subsurface_drainage'],
                   colors=pal, cum=True, stack=True, xticks=False,
                   unit_conversion={'unit':'mm', 'conversion':1e3}, legend=False, limits=False)
plt.xticks([])
plt.ylabel('')
_=plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1=plt.subplot(143, sharey=ax)
i=2
plot_timeseries_xr(results[i].sel(date=results[i]['date.year']==yyyy), ['canopy_moss_evaporation', 'canopy_transpiration', 'canopy_evaporation',
                   'soil_subsurface_drainage'],
                   colors=pal, cum=True, stack=True, xticks=False,
                   unit_conversion={'unit':'mm', 'conversion':1e3}, legend=False, limits=False)
plt.xticks([])
plt.ylabel('')
_=plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
plt.tight_layout()
plt.legend(bbox_to_anchor=(1.02,0.5), loc="center left", frameon=False, borderpad=0.0, fontsize=12)
#plt.savefig('WB.png',dpi=500)

""" PLANT TYPE TRANSPIRATION """
pal = sns.color_palette("hls", 7)[1:5]
pal = np.flipud(pal)
pal[1] = (pal[2]**1.5)
pal[0] = (pal[2]**3)
variable='canopy_pt_transpiration'
plt.figure(figsize=(7,1.7))
yyyy=2017
ax = plt.subplot(141)
i=0
plot_timeseries_pt(results[i].sel(date=results[i]['date.year']==yyyy), variable, unit_conversion={'unit':'mm', 'conversion':1e3},
                    xticks=True, stack=True, cum=True,legend=False, colors=pal)
plt.xticks([])
plt.title('')
plt.ylabel('')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax1=plt.subplot(142, sharey=ax)
i=1
plot_timeseries_pt(results[i].sel(date=results[i]['date.year']==yyyy), variable, unit_conversion={'unit':'mm', 'conversion':1e3},
                    xticks=True, stack=True, cum=True,legend=False, limits=False, colors=pal)
plt.xticks([])
plt.title('')
plt.ylabel('')
_=plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1=plt.subplot(143, sharey=ax)
i=2
plot_timeseries_pt(results[i].sel(date=results[i]['date.year']==yyyy), variable, unit_conversion={'unit':'mm', 'conversion':1e3},
                    xticks=True, stack=True, cum=True,legend=False, limits=False, colors=pal)
plt.xticks([])
plt.title('')
plt.ylabel('')
_=plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
plt.tight_layout()
plt.legend(bbox_to_anchor=(1.02,0.5), loc="center left", frameon=False, borderpad=0.0, fontsize=12)
#plt.savefig('Tr.png',dpi=500)



plot_results(results[1])
plot_results(results[2])

results['soil_pond_storage'].plot()


plt.figure()
results['canopy_energy_closure_canopy'].plot()
results['canopy_fr_source'].plot()
plt.figure(1)
results['canopy_Rabs'].mean(dim='date').plot()
results['canopy_LWleaf'].mean(dim='date').plot()
plt.figure(5)
i=1
results['canopy_Tleaf'].mean(dim='date').plot()
results['canopy_Tleaf_wet'].mean(dim='date').plot()
results['canopy_Tleaf_sl'].mean(dim='date').plot()
results['canopy_Tleaf_sh'].mean(dim='date').plot()
results['canopy_temperature'].mean(dim='date').plot()
plt.legend(['Tleaf','Tleaf_wet','Tleaf_sl','Tleaf_sh','Tair'])
plt.figure(4)
idx=25
results['canopy_Tleaf'].isel(date=idx).plot()
#results['canopy_Tleaf_wet'].isel(date=idx).plot()
#results['canopy_Tleaf_sl'].isel(date=idx).plot()
#results['canopy_Tleaf_sh'].isel(date=idx).plot()
results['canopy_T'].isel(date=idx).plot()
plt.legend(['Tleaf','Tleaf_wet','Tleaf_sl','Tleaf_sh','Tair'])
plt.figure(3)
results['canopy_throughfall'].plot()
results['canopy_Tleaf_wet'].mean(dim='date').plot()
results['canopy_lad'].mean(dim='date').plot()
plt.figure(3)
results['canopy_T'].isel(date=3).plot()
results['canopy_Tleaf'].isel(date=20).plot()

plt.figure()
results['canopy_IterWMA'].plot()
#results['forcing_h2o'].plot()
results['canopy_WMA_assumption'].plot()

results['canopy_WMA'].sum()/len(results['canopy_WMA'])

results = read_results(['results/201809051734_CCFPeat_results.nc',
                        'results/201809060950_CCFPeat_results.nc',
                        'results/201809060954_CCFPeat_results.nc',
                        'results/201809061000_CCFPeat_results.nc'])
labels=['Ebal & no WMA', 'no Ebal & no WMA', 'Ebal & WMA', 'no Ebal & WMA']
labels=['Lv','Ls']
plt.figure()
plot_timeseries_xr(results, 'canopy_IterWMA', labels=labels)

variable = 'canopy_moss_evaporation'
plt.figure()
plt.subplot(211)
plot_timeseries_xr(results, variable, unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
plt.subplot(212)
plot_timeseries_xr(results, variable, cum=True, unit_conversion={'unit':'mm', 'conversion':1e3})

plt.figure()
plt.subplot(211)
plot_timeseries_xr(results, 'canopy_evaporation', unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
plt.subplot(212)
plot_timeseries_xr(results, 'canopy_evaporation', cum=True, unit_conversion={'unit':'mm', 'conversion':1e3})
plt.figure()
plt.subplot(211)
plot_timeseries_xr(results, 'canopy_condensation', labels=labels, unit_conversion={'unit':'mm h-1', 'conversion':1e3*3600})
plt.subplot(212)
plot_timeseries_xr(results, 'canopy_condensation', labels=labels, cum=True, unit_conversion={'unit':'mm', 'conversion':1e3})

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

results=read_results('results/201804241836_CCFPeat_results.nc')

variable = 'soil_ground_water_level'

pal = sns.color_palette("hls", 7)

from parameters.sensitivity_sampling import LAIcombinations
plt.figure()
for i in range(len(LAIcombinations)):
    plot_timeseries_xr(results, variable, sim_idx=i, colors=pal[int(sum(LAIcombinations[i])):], labels=", stand LAI = " + str(sum(LAIcombinations[i])))

plt.subplot(211)

for i in range(len(LAIcombinations)):
    plot_timeseries_xr(results, variable, sim_idx=i, colors=pal[int(sum(LAIcombinations[i])):], labels=", stand LAI = " + str(sum(LAIcombinations[i])))
plt.subplot(212)
for i in range(len(LAIcombinations)):
    plot_timeseries_xr(results, variable, sim_idx=i, colors=pal[i:], labels=", stand LAI = " + str(sum(LAIcombinations[i])))
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
from canopy.evapotranspiration import latent_heat
letto_meteo = read_forcing("Lettosuo_meteo_2010_2018.csv", cols=['Tair'])
lettosuo_data['L'] = latent_heat(letto_meteo['Tair'].values)
lettosuo_data['ET_control'] = lettosuo_data['Letto1_EC: LE (Wm-2)'] / lettosuo_data['L'] * 3600  # [mm h-1]
lettosuo_data['ET_clearcut'] = lettosuo_data['avohakkuu_EC: latent heat flux: FL [W m-2]'] / lettosuo_data['L'] * 3600  # [mm h-1]
lettosuo_data['ET_partial'] = lettosuo_data['energyfluxes_lettosuo: LE [W m-2]'] / lettosuo_data['L'] * 3600  # [mm h-1]
lettosuo_data['ET_control'][lettosuo_data['ET_control'] < 0.0] = np.nan
lettosuo_data['ET_clearcut'][lettosuo_data['ET_clearcut'] < 0.0] = np.nan
lettosuo_data['ET_partial'][lettosuo_data['ET_partial'] < 0.0] = np.nan

letto_EC = lettosuo_data[['ET_clearcut',
                          'ET_control',
                          'ET_partial',
                          'Letto1_EC: GPP_S',
                          'Letto1_EC: GPP_N',
                          'avohakkuu_EC: GPP [mg CO2 m-2 s-1]',
                          'Partial_EC_gapfilled_fluxes: GPP [mg CO2 m-2 s-1]']].copy()

letto_EC = letto_EC.rename(columns={'Letto1_EC: GPP_S': 'GPP_control_S',
                                    'Letto1_EC: GPP_N': 'GPP_control_N',
                                    'avohakkuu_EC: GPP [mg CO2 m-2 s-1]': 'GPP_clearcut',
                                    'Partial_EC_gapfilled_fluxes: GPP [mg CO2 m-2 s-1]': 'GPP_partial'})

letto_EC['GPP_control_S'] = -letto_EC['GPP_control_S']
letto_EC['GPP_control_N'] = -letto_EC['GPP_control_S']
letto_EC['GPP_clearcut'] = -letto_EC['GPP_clearcut']
letto_EC['GPP_partial'] = -letto_EC['GPP_partial']

""" PLOT DIURNAL GPP AND ET"""
from tools.plotting import plot_diurnal
cc=[(255./255.,130./255,0.),(0.,181./255,226./255),(1.,0.,0.),(0.85,0.85,0.85)]
plt.figure(figsize=(7,2.2))
ax = plt.subplot(131)
for yyyy in np.linspace(2010,2015,6):
    plot_diurnal(letto_EC.GPP_control_S[letto_EC.index.year==yyyy], color=cc[3], legend=False)
plot_diurnal(letto_EC.GPP_control_S[letto_EC.index.year<=2015], color=cc[0], legend=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax1=plt.subplot(132, sharey=ax)
yyyy=2016
plot_diurnal(letto_EC.GPP_partial[letto_EC.index.year==yyyy], color=cc[1], legend=False)
plot_diurnal(letto_EC.GPP_clearcut[letto_EC.index.year==yyyy], color=cc[2], legend=False)
_=plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1=plt.subplot(133, sharey=ax)
yyyy=2017
plot_diurnal(letto_EC.GPP_partial[letto_EC.index.year==yyyy], color=cc[1], legend=False)
plot_diurnal(letto_EC.GPP_clearcut[letto_EC.index.year==yyyy], color=cc[2], legend=False)
_=plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
plt.tight_layout()#rect=(0, 0, 0.8, 1))
#plt.savefig('diurnal_GPP.png',dpi=500)


letto_meteo = read_forcing("Lettosuo_meteo_2010_2018.csv", cols=['Prec'])
dates = letto_meteo.index
ix = pd.rolling_sum(letto_meteo.Prec.values, 48, 1)
dryc = np.ones(len(dates))
f = np.where(ix > 1.0)[0]  # wet canopy indices
dryc[f] = 0.0
months = letto_meteo.index.month
fmonth = 5
lmonth = 9
f = np.where((months >= fmonth) & (months <= lmonth) & (dryc == 1))[0]
ff = np.zeros(len(dates))
ff[f] = 1.0

#plt.figure()
##letto_meteo.Prec.plot()
##letto_EC.ET_control.plot()
##letto_EC.ET_control[(ff==1)].plot()
#letto_EC.ET_clearcut.plot()
##letto_EC.ET_clearcut[(ff==1)].plot()
#letto_EC.ET_partial.plot()
##letto_EC.ET_partial[(ff==1)].plot()



plt.figure(figsize=(7,2.2))
ax = plt.subplot(131)
for yyyy in np.linspace(2010,2015,6):
    plot_diurnal(letto_EC.ET_control[(letto_EC.index.year==yyyy) & (ff==1)], color=cc[3], legend=False)
plot_diurnal(letto_EC.ET_control[(letto_EC.index.year<=2015) & (ff==1)], color=cc[0], legend=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax1=plt.subplot(132, sharey=ax)
yyyy=2016
plot_diurnal(letto_EC.ET_partial[(letto_EC.index.year==yyyy) & (ff==1)], color=cc[1], legend=False)
plot_diurnal(letto_EC.ET_clearcut[(letto_EC.index.year==yyyy) & (ff==1)], color=cc[2], legend=False)
_=plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1=plt.subplot(133, sharey=ax)
yyyy=2017
plot_diurnal(letto_EC.ET_partial[(letto_EC.index.year==yyyy) & (ff==1)], color=cc[1], legend=False)
plot_diurnal(letto_EC.ET_clearcut[(letto_EC.index.year==yyyy) & (ff==1)], color=cc[2], legend=False)
_=plt.setp(plt.gca().axes.get_yticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
plt.tight_layout()#rect=(0, 0, 0.8, 1))
#plt.savefig('diurnal_ET.png',dpi=500)

""" PLOT WTD """
gwl_meas = read_forcing("lettosuo_WTD_pred.csv", cols='all')
plt.figure(figsize=(7,2))
ax=plt.subplot(111)
plot_timeseries_df(gwl_meas[gwl_meas.index >= "11.1.2015"], ['ctrl','part','clear',],colors=cc,xticks=True, 
                   limits=False, legend=False)
plot_timeseries_df(gwl_meas[gwl_meas.index < "11.1.2015"], ['part'],colors=[cc[0]],xticks=True, 
                   limits=False, legend=False)
plt.xlim("5.1.2010","5.1.2018")
plt.ylim(-0.85,0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.yticks(np.linspace(-0.8,0.0,3), np.linspace(0.8,0.0,3))
plt.tight_layout()#rect=(0, 0, 0.8, 1))
plt.savefig('WTD.png',dpi=500)

""" PLOT LAD PROFILES """
plot_lad_profiles()

""" WTD PROCESSING"""
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