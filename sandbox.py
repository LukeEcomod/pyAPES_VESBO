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
#outputfile = 'results/201803201103_CCFPeat_results.nc'

filepath='C:/Users/L1656/Documents/Git_repos/CCFPeat/' + outputfile
results=xr.open_dataset(filepath)

results.coords['simulation']=results.simulation.values
results.coords['soil']=results.soil_z.values
results.coords['canopy']=results.canopy_z.values

plotresults(results)

plotresultsMLM(results)

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
lad_p, lad_s, lad_d, _, _, _ = model_trees(z, quantiles, normed=False,
    dbhfile=r"C:\Users\L1656\Documents\Git_repos\CCFPeat\parameters\runkolukusarjat\letto2014.txt",
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