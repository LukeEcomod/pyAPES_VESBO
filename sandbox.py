# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:40:57 2018

@author: L1656
    """
import pandas as pd
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from plotting import plotresults, plotxarray, plotresultsMLM

#outputfile=driver(create_ncf=True)
#outputfile = 'results/201804241836_CCFPeat_results.nc'

filepath='C:/Users/L1656/Documents/Git_repos/CCFPeat/' + outputfile
results=xr.open_dataset(filepath)

results.coords['simulation']=results.simulation.values
results.coords['soil']=results.soil_z.values
results.coords['canopy']=results.canopy_z.values

plotresults(results.isel(simulation=0))

plotresultsMLM(results.isel(simulation=0))

import seaborn as sns
pal = sns.color_palette("hls", 6)
from parameters.sensitivity_sampling import LAIcombinations
for i in range(len(LAIcombinations)):
    k=int(sum(LAIcombinations[i]))-1
    plotxarray(results.isel(simulation=i), ['soil_ground_water_level'], colors=pal[k:], xticks=True)

"""-------- Colormesh plotting --------"""

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


# readind lettosuo data
# filepaths
forc_fp = ["H:/Lettosuo/Forcing_data/Annalea1/avohakkuu_EC.csv",
           "H:/Lettosuo/Forcing_data/Annalea2/Letto1_EC.csv",
           "H:/Lettosuo/Forcing_data/Annalea1/avohakkuu_meteo.csv",
           "H:/Lettosuo/Forcing_data/Annalea1/Letto1_meteo.csv",
           "H:/Lettosuo/Forcing_data/Annalea1/Letto1_metsanpohja.csv",
           "H:/Lettosuo/Forcing_data/Annalea1/Letto2_metsanpohja.csv",
           "H:/Lettosuo/Forcing_data/FMI_jokioinen/jokioinen_meteo.csv",
           "H:/Lettosuo/Forcing_data/FMI_jokioinen/jokioinen_rad.csv",
           "H:/Lettosuo/Forcing_data/Annalea2/Letto1_meteo_gapfilled.csv"]

index=pd.date_range('01-01-2009','01-01-2018',freq='0.5H')
lettosuo_data=pd.DataFrame(index=index, columns=[])

for fp in forc_fp:
    dat = pd.read_csv(fp, sep=',', header='infer')
    if fp.split("/")[-2] == 'Annalea2':
        dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
        # period end
        dat.index = dat.index - pd.Timedelta(hours=0.5)
    else:
        dat.index = pd.to_datetime(dat.ix[:,0], yearfirst=True)
    if fp.split("/")[-2] == 'FMI_jokioinen':
        # UTC -> UTC + 2
        dat.index = dat.index + pd.Timedelta(hours=2)
    dat.index=dat.index.map(lambda x: x.replace(second=0))
    dat.ix[:,0]=dat.index
    dat = dat.drop_duplicates(subset=dat.columns[0])
    dat = dat.drop(dat.columns[0], axis=1)
#    dat.plot(kind='line',subplots=True,title=fp)
    dat.columns = fp.split("/")[-1].split(".")[0] + ': ' + dat.columns
    lettosuo_data=lettosuo_data.merge(dat, how='outer', left_index=True, right_index=True)

# divide hourly precipitation to half hour
lettosuo_data.ix[1:-1:2,89]=lettosuo_data.ix[0:-1:2,89].values
lettosuo_data.ix[:,89]=lettosuo_data.ix[:,89].values/2.0

def plot_columns(data, col_index):
    col_names=[]
    for i in col_index:
        col_names.append(data.columns[i])
    data[col_names].plot(kind='line',marker='o',markersize=1)
    plt.legend()
    pd.plotting.scatter_matrix(data[col_names], figsize=(10, 10), alpha=.2)

# RH
plot_columns(lettosuo_data,[44,57,63,81,92,100])

# Tair
plot_columns(lettosuo_data,[43,58,69,88,90,101])

# Prec
plot_columns(lettosuo_data,[89,64])

# wind
plot_columns(lettosuo_data,[19,23,91])

# ustar
plot_columns(lettosuo_data,[7,24])

# Global radiation
plot_columns(lettosuo_data,[97,52,95])

# Pressure
plot_columns(lettosuo_data,[59,93])

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

# lad profile test
lad_p, lad_s, lad_d, _, _, _, _, _, _ = model_trees(z, quantiles, normed=False,
    dbhfile=r"C:\Users\L1656\Documents\Git_repos\CCFPeat\parameters\runkolukusarjat\letto2016_partial.txt",
    plot=True)



#from canopy.evapotranspiration import e_sat
#
#esat, _, _ = e_sat(results.forcing_air_temperature.values)  # [Pa]
#RH = results.forcing_h2o.values * 101300 / esat
#vpd = (1 - RH) * esat * 1e-3  # [kPa]
#plt.plot(vpd)
##plotxarray(results, ['forcing_h2o'], colors=pal, xticks=False)
#plt.figure()
#plotxarray(results, ['canopy_Rnet'], colors=pal, xticks=False)
#plotxarray(results1, ['canopy_Rnet'], colors=pal[1:], xticks=False)
#
#plt.figure()
#plotxarray(results1, ['canopy_Rnet_ground'], colors=pal[1:], xticks=False)
#plotxarray(results, ['canopy_Rnet_ground'], colors=pal, xticks=False)
#
#plt.figure()
#plotxarray(results1, ['canopy_U_ground'], colors=pal[1:], xticks=False)
#plotxarray(results, ['canopy_U_ground'], colors=pal, xticks=False)
#
#plt.figure()
#plotxarray(results, ['forcing_h2o'], colors=pal, xticks=False)
#plotxarray(results1, ['forcing_h2o'], colors=pal[1:], xticks=False)

def create_hyde_forcingfile():

    """ read files in format Forcing_YYYY.dat """
    direc = "H:/Samulilta/Hyde/Forcing_corrected/"
    frames = []
    columns =['U','ust','Ta','RH','CO2','H2O','Prec','P',
              'dirPar','diffPar','dirNir','diffNir','Rnet',
              'LWin','LWout','LWnet','Tsh','Tsa','Tsc','Wh','Wa']
    for year in range(1997,2017):
        forc_fp = direc + "Forcing_" + str(year) + ".dat"
        dat = pd.read_csv(forc_fp, sep=',', header='infer')
        index = pd.date_range('01-01-' + str(year),'31-12-' + str(year) + ' 23:59:59',freq='0.5H')
        if len(index) != len(dat):
            print "Year " + str(year) + ": Missing values!"
        else:
            dat.index = index
            frames.append(dat[columns])

    hyde_data=pd.concat(frames)

    hyde_data.insert(0,'yyyy',hyde_data.index.year.values)
    hyde_data.insert(1,'mo',hyde_data.index.month.values)
    hyde_data.insert(2,'dd',hyde_data.index.day.values)
    hyde_data.insert(3,'hh',hyde_data.index.hour.values)
    hyde_data.insert(4,'mm',hyde_data.index.minute.values)

    """ read files in format FIHy_YYYY_pd.csv """
    direc = "H:/Samulilta/Hyde/FIHy_1997_2016_new/"
    frames = []
    columns =['NEE','GPP','LE','ET','fRg']
    for year in range(1997,2017):
        forc_fp = direc + "FIHy_" + str(year) + "_pd.csv"
        dat = pd.read_csv(forc_fp, sep=',', header='infer')
        index = pd.date_range('01-01-' + str(year),'31-12-' + str(year) + ' 23:59:59',freq='0.5H')
        if len(index) != len(dat):
            print "Year " + str(year) + ": Missing values!"
        else:
            dat.index = index
            frames.append(dat[columns])

    df = pd.concat(frames)
    hyde_data=hyde_data.merge(df, how='outer', left_index=True, right_index=True)
    hyde_data=hyde_data.rename(columns={'Ta':'Tair', 'ust':'Ustar', 'fRg':'Rg'})

    fp = "C:/Users/L1656/Documents/Git_repos/CCFPeat/forcing/Hyde_data_1997_2016.csv"
    hyde_data.to_csv(path_or_buf=fp, sep=',', na_rep='NaN',index=False)
