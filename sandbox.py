# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:40:57 2018

@author: L1656
    """
import pandas as pd
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from plotting import plotresults, plotxarray, plotxarray2, plotresultsMLM, plot_columns
from timeseries_tools import fill_gaps
import datetime

#outputfile=driver(create_ncf=True)
#outputfile = 'results/201804241836_CCFPeat_results.nc'

def read_results(outputfiles):
    direc='C:/Users/L1656/Documents/Git_repos/CCFPeat/'
    if type(outputfiles) != list:
        outputfiles = [outputfiles]
    results = []
    for outputfile in outputfiles:
        fp =direc + outputfile
        result=xr.open_dataset(fp)
        result.coords['simulation']=result.simulation.values
        result.coords['soil']=result.soil_z.values
        result.coords['canopy']=result.canopy_z.values
        results.append(result)
    if len(results) == 1:
        return results[0]
    else:
        return results

plotresults(results.isel(simulation=0))

plotresultsMLM(results.isel(simulation=0))

output_control1 = driver(create_ncf=True, dbhfile="letto2009.txt")
output_control2 = driver(create_ncf=True, dbhfile="letto2014.txt")
output_partial = driver(create_ncf=True, dbhfile="letto2016_partial.txt")
output_clearcut = driver(create_ncf=True, dbhfile="letto2016_clearcut.txt")

results = read_results([output_control1, output_control2, output_partial, output_clearcut])
pal = sns.color_palette("hls", 5)
plotxarray2(results, 'soil_ground_water_level', colors=pal, xticks=True)

import seaborn as sns
pal = sns.color_palette("hls", 6)
from parameters.sensitivity_sampling import LAIcombinations
for i in range(len(LAIcombinations)):
    k=int(sum(LAIcombinations[i]))-1
    plotxarray(results.isel(simulation=i), ['soil_ground_water_level'], colors=pal[k:], xticks=True)

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


# readind data
# filepaths
forc_fp = ["H:/Lettosuo/Forcing_data/Annalea1/avohakkuu_EC.csv",
           "H:/Lettosuo/Forcing_data/Annalea2/Letto1_EC.csv",
           "H:/Lettosuo/Forcing_data/Annalea1/avohakkuu_meteo.csv",
           "H:/Lettosuo/Forcing_data/Annalea1/Letto1_meteo.csv",
           "H:/Lettosuo/Forcing_data/Annalea1/Letto1_metsanpohja.csv",
           "H:/Lettosuo/Forcing_data/Annalea1/Letto2_metsanpohja.csv",
           "H:/Lettosuo/Forcing_data/FMI/jokioinen_meteo.txt",
           "H:/Lettosuo/Forcing_data/FMI/jokioinen_rad.txt",
           "H:/Lettosuo/Forcing_data/Annalea2/Letto1_meteo_gapfilled.csv",
           "H:/Lettosuo/Forcing_data/FMI/jokioinen_prec1.txt",
           "H:/Lettosuo/Forcing_data/FMI/jokioinen_prec2.txt",
           "H:/Lettosuo/Forcing_data/FMI/somero_meteo.txt",
           "H:/Lettosuo/Forcing_data/FMI/hameenlinna_meteo.txt",
           "H:/Lettosuo/Forcing_data/FMI/salo_kiikala_meteo.txt"]

index=pd.date_range('01-01-2009','06-01-2018',freq='0.5H')
lettosuo_data=pd.DataFrame(index=index, columns=[])

for fp in forc_fp:
    dat = pd.read_csv(fp, sep=',', header='infer')
    if fp.split("/")[-2] == 'Annalea2':
        dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
        # period end
        dat.index = dat.index - pd.Timedelta(hours=0.5)
    else:
        dat.index = pd.to_datetime(dat.ix[:,0], yearfirst=True)
    if fp.split("/")[-2] == 'FMI':
        # UTC -> UTC + 2
        dat.index = dat.index + pd.Timedelta(hours=2)
    dat.index=dat.index.map(lambda x: x.replace(second=0))
    dat.ix[:,0]=dat.index
    dat = dat.drop_duplicates(subset=dat.columns[0])
    dat = dat.drop(dat.columns[0], axis=1)
    dat.columns = fp.split("/")[-1].split(".")[0] + ': ' + dat.columns
    lettosuo_data=lettosuo_data.merge(dat, how='outer', left_index=True, right_index=True)

# divide hourly precipitation to half hour
for i in [89, 107, 108, 109, 115, 121]:
    lettosuo_data.ix[:,i]=lettosuo_data.ix[:,i].replace(-1,0)
    lettosuo_data.ix[1:-1:2,i]=lettosuo_data.ix[0:-1:2,i].values
    lettosuo_data.ix[:,i]=lettosuo_data.ix[:,i].values/2.0

lettosuo_data = lettosuo_data[(lettosuo_data.index >= '01-01-2010') & 
                              (lettosuo_data.index <= '01-01-2018')]

gap_fill_lettosuo_meteo(lettosuo_data)
create_lettosuo_forcingfile()



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

idx = np.isfinite(lettosuo_data['salo_kiikala_meteo: Pressure (msl)']) & np.isfinite(lettosuo_data['Letto1_meteo: avg(Press (hPa))'])
diff = lettosuo_data['salo_kiikala_meteo: Pressure (msl)'][idx].values - lettosuo_data['Letto1_meteo: avg(Press (hPa))'][idx].values
plt.figure()
plt.plot(diff)
diff.mean()

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

""" --- inspecting precipitation data ---- """

def continuous_prec(Prec_data):

    """ --- continuous prec timeseries for jokioinen --- """
    # FMI open data jokioinen + FMI jokioinen gauge 1
    Prec_data['jokioinen_prec'] = Prec_data['jokioinen_meteo: Precipitation amount'].fillna(Prec_data['jokioinen_prec1: Prec [mm h-1]'])
    Prec_data['jokioinen_prec_flag']=np.where(Prec_data['jokioinen_prec'].isnull(), 10.0, 0.0)
    # + FMI jokioinen gauge 2
    Prec_data['jokioinen_prec'] = Prec_data['jokioinen_prec'].fillna(Prec_data['jokioinen_prec2: Prec [mm h-1]'])
    Prec_data['jokioinen_prec_flag']=np.where(Prec_data['jokioinen_prec'].isnull(), 20.0, Prec_data['jokioinen_prec_flag'])
    # + FMI open data somero
    Prec_data['jokioinen_prec'] = Prec_data['jokioinen_prec'].fillna(Prec_data['somero_prec: Precipitation amount'])
    Prec_data['jokioinen_prec_flag']=np.where(Prec_data['jokioinen_prec'].isnull(), 30.0, Prec_data['jokioinen_prec_flag'])
    # fill rest with zero
    Prec_data['jokioinen_prec'] = Prec_data['jokioinen_prec'].fillna(0)
    # plot
    #plot_columns(Prec_data[['jokioinen_prec_flag','jokioinen_prec']])
    
    """ --- continuous prec timeseries for lettosuo --- """
    Prec_data['lettosuo_prec_flag']=np.where(Prec_data['Letto1_metsanpohja: avg(Rain (mm))'].isnull(), 10.0, 0.0)
    # consider prec data unrealiable when Tair < 2C
    Prec_data['lettosuo_prec_flag']=np.where(Prec_data['Letto1_metsanpohja: avg(Temp (C))'] < 2.0,
             20.0, Prec_data['lettosuo_prec_flag'])
    Prec_data['lettosuo_prec']=np.where(Prec_data['Letto1_metsanpohja: avg(Temp (C))'] < 2.0,
             np.nan, Prec_data['Letto1_metsanpohja: avg(Rain (mm))'])
    # if rolling 3 day mean prec less than 10% of jokionen rolling mean, remove
    Prec_data_daily = pd.rolling_sum(Prec_data.fillna(0), 3 * 48, 1)
    Prec_data_daily['lettosuo_prec'] = np.where(Prec_data['lettosuo_prec'].isnull(), np.nan, Prec_data_daily['lettosuo_prec'])
    Prec_data['lettosuo_prec_flag']=np.where(Prec_data_daily['lettosuo_prec'] < 0.1 * Prec_data_daily['jokioinen_prec'],
             30.0, Prec_data['lettosuo_prec_flag'])
    Prec_data['lettosuo_prec']=np.where(Prec_data_daily['lettosuo_prec'] < 0.1 * Prec_data_daily['jokioinen_prec'],
             np.nan, Prec_data['lettosuo_prec'])
    # fill gaps with jokioinen data
    Prec_data['lettosuo_prec'] = Prec_data['lettosuo_prec'].fillna(Prec_data['jokioinen_prec'])
    # plot
    #plot_columns(Prec_data[['lettosuo_prec_flag','lettosuo_prec']])

    return Prec_data

# readind daily FMI prec
Prec_data = continuous_prec(lettosuo_data[['Letto1_metsanpohja: avg(Temp (C))',
                                           'Letto1_metsanpohja: avg(Rain (mm))',
                                           'jokioinen_meteo: Precipitation amount',
                                           'jokioinen_prec1: Prec [mm h-1]',
                                           'jokioinen_prec2: Prec [mm h-1]',
                                           'somero_prec: Precipitation amount']])

# filepaths
forc_fp = ["H:/Lettosuo/Forcing_data/FMI_jokioinen/jokioinen_daily_prec.txt",
           "H:/Lettosuo/Forcing_data/FMI_jokioinen/somero_daily_prec.txt"]

index=pd.date_range('01-01-2009','01-01-2018',freq='1.0D')
daily_FMI_prec=pd.DataFrame(index=index, columns=[])

for fp in forc_fp:
    dat = pd.read_csv(fp, sep=',', header='infer')
    dat = dat.replace(-1,0)
    dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
    dat = dat.drop(dat.columns[0], axis=1)
    dat.columns = fp.split("/")[-1].split(".")[0] + ': ' + dat.columns
    daily_FMI_prec=daily_FMI_prec.merge(dat, how='outer', left_index=True, right_index=True)

# append daily series from gapfilled filled hourly data[['jokioinen_prec','lettosuo_prec']]
Prec_daily = pd.rolling_sum(Prec_data, 48, 1)
Prec_daily = Prec_daily[Prec_daily.index.hour == 8.0]
Prec_daily = Prec_daily[Prec_daily.index.minute == 0.0]
Prec_daily.index = Prec_daily.index - pd.Timedelta(hours=24+8)

Prec_daily=Prec_daily.merge(daily_FMI_prec, how='outer', left_index=True, right_index=True)

plot_columns(Prec_daily,[0,2,3,4,8,9,10])

from timeseries_tools import yearly_cumulative

var = [#'jokioinen_prec',
       'lettosuo_prec',
       'jokioinen_daily_prec: Precipitation amount',
#       'somero_daily_prec: Precipitation amount',
       'jokioinen_prec1: Prec [mm h-1]',
       'jokioinen_prec2: Prec [mm h-1]']
Prec = Prec_daily.fillna(0)
Prec_cum = yearly_cumulative(Prec, var)

plt.figure()
for i in range(len(Prec_cum)):
    plt.plot(Prec_daily.index, Prec_cum[i,:], label=var[i])
plt.legend()

def plot_xy(x,y,lim):
    plt.scatter(x, y, marker='o', alpha=.2)
    idx = np.isfinite(x) & np.isfinite(y)
    p = np.polyfit(x[idx], y[idx], 1)
    corr = np.corrcoef(x[idx], y[idx])
    plt.annotate("y = %.2fx + %.2f \nR2 = %.2f" % (p[0], p[1], corr[1,0]**2), (0.3, 0.9), xycoords='axes fraction', ha='center', va='center')
    plt.plot(lim, [p[0]*lim[0] + p[1], p[0]*lim[1] + p[1]], 'r', linewidth=1)
    plt.plot(lim, lim, 'k--', linewidth=1)
    plt.ylim(lim)
    plt.xlim(lim)

ax1 = plt.subplot(131)
lim = [0, 45]
plot_xy(Prec_daily['jokioinen_prec2: Prec [mm h-1]'][Prec_daily.index.year < 2013],
        Prec_daily['jokioinen_daily_prec: Precipitation amount'][Prec_daily.index.year < 2013],
        lim)
ax2 = plt.subplot(132)
plot_xy(Prec_daily['jokioinen_prec2: Prec [mm h-1]'][Prec_daily.index.year > 2013],
        Prec_daily['jokioinen_daily_prec: Precipitation amount'][Prec_daily.index.year > 2013],
        lim)
ax3 = plt.subplot(133)
plot_xy(Prec_daily['jokioinen_prec1: Prec [mm h-1]'][Prec_daily.index.year > 2013],
        Prec_daily['jokioinen_daily_prec: Precipitation amount'][Prec_daily.index.year > 2013],
        lim)

ax2.set_yticklabels([])
ax3.set_yticklabels([])
ax1.set_ylabel('Daily open FMI Jokioinen (mm/d)')
ax1.set_xlabel('Jokioinen gauge 50 before 2013 (mm/d)')
ax2.set_xlabel('Jokioinen gauge 50 after 2013 (mm/d)')
ax3.set_xlabel('Jokioinen gauge 60 after 2013 (mm/d)')

labels=['Lettosuo filled','Daily open FMI Jokioinen','Jokioinen gauge 50','Jokioinen gauge 60']
plt.figure()
for i in reversed(range(len(Prec_cum))):
    if i == 1:
        plt.plot(Prec_daily.index, Prec_cum[i,:], 'k:', label=labels[i])
    else:
        plt.plot(Prec_daily.index, Prec_cum[i,:], label=labels[i])
plt.legend()
plt.plot(Prec_daily.index, Prec_cum[1,:], 'k:', label=labels[1])

plt.ylabel('Cumulative precipitation (mm)')

# lad profile test
from parameters.parameter_utils import model_trees

def plot_lad(filename="letto2016_partial.txt"):

    z = np.linspace(0, 30.0, 100)
    quantiles = [1.0]
    lad_p, lad_s, lad_d, _, _, _, _, _, _ = model_trees(z, quantiles, normed=False,
        dbhfile="C:/Users/L1656/Documents/Git_repos/CCFPeat/parameters/runkolukusarjat/" + filename,
        plot=True)

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

def gap_fill_lettosuo_meteo(lettosuo_data):
    """ gap filling meteo """
    frames = []
    readme = ""
    plot = False
    
    # --- Precipitation --- 
    # Jokioinen
    df, info = fill_gaps(lettosuo_data[['jokioinen_prec1: Prec [mm h-1]',
                                        'jokioinen_prec2: Prec [mm h-1]',
                                        'somero_meteo: Precipitation amount']],
                         'Prec_ref', 'Jokioinen gapfilled precipitaion [mm/30min]', 
                         fill_nan = 0.0, plot=plot)
    frames.append(df)
    readme += info
    
    # Lettosuo
    # consider prec data unrealiable when Tair < 2C
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!']=np.where(
            lettosuo_data['Letto1_metsanpohja: avg(Temp (C))'] < 2.0,
            np.nan, lettosuo_data['Letto1_metsanpohja: avg(Rain (mm))'])
    # if rolling 3 day mean prec less than 10% of jokionen rolling mean, remove
    Prec_daily_ref = pd.rolling_sum(df[['Prec_ref']], 3 * 48, 1)
    Prec_daily= pd.rolling_sum(lettosuo_data[['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!']].fillna(0), 3 * 48, 1)
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!']=np.where(
            Prec_daily['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'] < 0.1 * Prec_daily_ref['Prec_ref'],
            np.nan, lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'])
    
    df, info = fill_gaps(lettosuo_data[['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!',
                                        'jokioinen_prec1: Prec [mm h-1]',
                                        'jokioinen_prec2: Prec [mm h-1]',
                                        'somero_meteo: Precipitation amount']],
                         'Prec', 'Lettosuo gapfilled precipitaion [mm/30min]',
                         fill_nan = 0.0, plot=plot)
    frames.append(df)
    readme += info
    
    # --- Air temperature --- 
    df, info = fill_gaps(lettosuo_data[['Letto1_meteo: avg(Temp (C))',
                                        'somero_meteo: Air temperature',
                                        'jokioinen_meteo: Air temperature',
                                        'Letto1_meteo_gapfilled: PaikAirT T']],
                         'Tair', 'Air temperature [degC]', fill_nan = 'linear', plot=plot)
    frames.append(df)
    readme += info
    
    # --- Relative humidity --- 
    df, info = fill_gaps(lettosuo_data[['Letto1_meteo: avg(RH (%))',
                                        'somero_meteo: Relative humidity',
                                        'jokioinen_meteo: Relative humidity',
                                        'Letto1_meteo_gapfilled: Paik RH']],
                         'RH', 'Relative humidity [%]', fill_nan = 'linear', plot=plot)
    # Check for values greater than 100%
    df['RH'][df['RH'] > 100.0] = 100.0
    frames.append(df)
    readme += info
    
    # --- Global radiation --- 
    df, info = fill_gaps(lettosuo_data[['Letto1_meteo: avg(Glob (W/m2))',
                                        'Letto1_meteo_gapfilled: PaikGlob2',
                                        'jokioinen_rad: Global radiation']],
                         'Rg', 'Global radiation [W/m2]', fill_nan = 'linear', plot=plot)
    df['Rg'][df['Rg'] < 0.0] = 0.0
    frames.append(df)
    readme += info
    
    # --- Wind speed --- 
    lettosuo_data['Letto1_EC: wind speed (m/s) !u > 10 removed!']=np.where(
            lettosuo_data['Letto1_EC: wind speed (m/s)'] > 10.0,
            np.nan, lettosuo_data['Letto1_EC: wind speed (m/s)'])
    lettosuo_data['salo_kiikala_meteo: Wind speed'][lettosuo_data['salo_kiikala_meteo: Wind speed'] == 0.0] = np.nan
    lettosuo_data['Derived from salo_kiikala_meteo: wind speed (0.57*U_ref + 0.55)'] = 0.57 * lettosuo_data['salo_kiikala_meteo: Wind speed'] + 0.55
    df, info = fill_gaps(lettosuo_data[['Letto1_EC: wind speed (m/s) !u > 10 removed!',
                                        'somero_meteo: Wind speed',
                                        'Derived from salo_kiikala_meteo: wind speed (0.57*U_ref + 0.55)']],
                         'U', 'Wind speed [m/s]', fill_nan = 'linear', plot=plot)
    frames.append(df)
    readme += info
    
    # --- Friction velocity --- 
    lettosuo_data['Ustar = 0.2 * U'] = 0.2 * df['U']
    df, info = fill_gaps(lettosuo_data[['Letto1_EC: friction velocity (m/s)',
                                        'Ustar = 0.2 * U']],
                         'Ustar', 'Friction velocity [m/s]', fill_nan = 'linear', plot=plot)
    frames.append(df)
    readme += info
    
    # --- Ambient pressure --- 
    lettosuo_data['Derived from salo_kiikala_meteo: Pressure (msl) (P_ref - 16.9)'] = lettosuo_data['salo_kiikala_meteo: Pressure (msl)'] - 16.9
    df, info = fill_gaps(lettosuo_data[['Letto1_meteo: avg(Press (hPa))',
                                        'Derived from salo_kiikala_meteo: Pressure (msl) (P_ref - 16.9)']],
                         'P', 'Ambient pressure [hPa]', fill_nan = 'linear', plot=plot)
    frames.append(df)
    readme += info
    
    letto_data=pd.concat(frames, axis=1)
    letto_data[['Prec_ref', 'Prec', 'Tair', 'Rg', 'U', 'Ustar', 'RH', 'P']].plot(subplots=True,kind='line')
    
    letto_data.insert(0,'yyyy',letto_data.index.year.values)
    letto_data.insert(1,'mo',letto_data.index.month.values)
    letto_data.insert(2,'dd',letto_data.index.day.values)
    letto_data.insert(3,'hh',letto_data.index.hour.values)
    letto_data.insert(4,'mm',letto_data.index.minute.values)
    
    fp = "C:/Users/L1656/Documents/Git_repos/CCFPeat/forcing/"
    fn = "Lettosuo_meteo_2010_2018"
    letto_data.to_csv(path_or_buf=fp + fn + ".csv", sep=',', na_rep='NaN', index=False)
    Readme = "Readme for " + fn + ".csv"
    Readme += "\n\nKersti Haahti, Luke " + str(datetime.datetime.now().date())
    Readme += "\n\nyyyy, mo, dd, hh, mm: datetime [UTC + 2.0]" + readme
    outF = open(fp + fn + "_readme.txt", "w")
    print >>outF, Readme
    outF.close()

from canopy.radiation import solar_angles, compute_clouds_rad
from canopy.evapotranspiration import e_sat

def create_lettosuo_forcingfile():
    """ Lettosuo forcing file from meteo"""

    fpar = 0.45
    lat = 60.63
    lon = 23.95

    forc_fp ="C:/Users/L1656/Documents/Git_repos/CCFPeat/forcing/"
    fn_meteo = "Lettosuo_meteo_2010_2018.csv"
    dat = pd.read_csv(forc_fp + fn_meteo, sep=',', header='infer')

    # set to dataframe index
    dat.index = pd.to_datetime({'year': dat['yyyy'],
                                'month': dat['mo'],
                                'day': dat['dd'],
                                'hour': dat['hh'],
                                'minute': dat['mm']})

    readme = "\n\nyyyy, mo, dd, hh, mm: datetime [UTC + 2.0]"
    cols = ['yyyy','mo','dd','hh','mm']

    # day of year
    dat['doy'] = dat.index.dayofyear
    cols.append('doy')
    readme += "\ndoy: Day of year [days]"

    # precipitaion unit from [mm/dt] to [m/s]
    dt = (dat.index[1] - dat.index[0]).total_seconds()
    dat['Prec'] = dat['Prec'] * 1e-3 / dt
    cols.append('Prec')
    readme += "\nPrec: Precipitation [m/s]"

    # atm. pressure unit from [hPa] to [Pa]
    dat['P'] = dat['P'] * 1e2
    cols.append('P')
    readme += "\nP: Ambient pressure [Pa]"

    # air temperature: instant and daily [degC]
    cols.append('Tair')
    readme += "\nTair: Air temperature [degC]"
    dat['Tdaily'] = pd.rolling_mean(dat['Tair'], int((24*3600)/dt), 1)
    cols.append('Tdaily')
    readme += "\nTdaily: Daily air temperature [degC]"

    # wind speend and friction velocity
    cols.append('U')
    readme += "\nU: Wind speed [m/s]"
    cols.append('Ustar')
    readme += "\nUstar: Friction velocity [m/s]"

    # ambient H2O [mol/mol] from RH
    esat, _, _ = e_sat(dat['Tair'])
    dat['H2O'] = (dat['RH'] / 100.0) * esat / dat['P']
    cols.append('H2O')
    readme += "\nH2O: Ambient H2O [mol/mol]"

    # ambient CO2 [ppm]
    dat['CO2'] = 380.0
    cols.append('CO2')
    readme += "\nCO2: Ambient CO2 [ppm] - set constant!"

    # zenith angle
    jday = dat.index.dayofyear + dat.index.hour / 24.0 + dat.index.minute / 1440.0
    dat['Zen'], _, _, _, _, _ = solar_angles(lat, lon, jday, timezone=+2.0)
    cols.append('Zen')
    readme += "\nZen: Zenith angle [rad]"

    # radiation components

    # global radiation
    cols.append('Rg')
    readme += "\nRg: Global radiation [W/m2]"

    # estimated long wave budget
    f_cloud, f_diff, emi_sky = compute_clouds_rad(dat['doy'].values,
                                                  dat['Zen'].values,
                                                  dat['Rg'].values,
                                                  dat['H2O'].values * dat['P'].values)
    b = 5.6697e-8  # Stefan-Boltzman constant (W m-2 K-4)
    dat['LWin'] = 0.98 * emi_sky * b *(dat['Tair'] + 273.15)**4 # Wm-2 downwelling LW
    dat['LWout'] = 0.98 * b * (dat['Tair'] + 273.15)**4  # Wm-2 upwelling LW
    cols.extend(('LWin', 'LWout'))
    readme += "\nLWin: Downwelling long wave radiation [W/m2]"
    readme += "\nLWout: Upwelling long wave radiation [W/m2]"

    # Short wave radiation; separate direct and diffuse PAR & NIR
    dat['diffPar'] = f_diff * fpar * dat['Rg']
    dat['dirPar'] = (1 - f_diff) * fpar * dat['Rg']
    dat['diffNir'] = f_diff * (1 - fpar) * dat['Rg']
    dat['dirNir'] = (1 - f_diff) * (1 - fpar) * dat['Rg']
    cols.extend(('diffPar', 'dirPar', 'diffNir', 'dirNir'))
    readme += "\ndiffPar: Diffuse PAR [W/m2] \ndirPar: Direct PAR [W/m2]"
    readme += "\ndiffNir: Diffuse NIR [W/m2] \ndirNir: Direct NIR [W/m2]"

    # approximate net radiation [W/m2]
    dat['Rnet'] = (1.0 - 0.08) * (dat['dirPar'] + dat['diffPar'] \
                   + dat['dirNir'] + dat['diffNir']) \
                   + dat['LWin'] - dat['LWout']
    cols.append('Rnet')
    readme += "\nRnet: Net radiation [W/m2]"

    dat = dat[cols]
    dat[cols[5:]].plot(subplots=True,kind='line')

    fp = "C:/Users/L1656/Documents/Git_repos/CCFPeat/forcing/"
    fn = "Lettosuo_forcing_2010_2018"
    dat.to_csv(path_or_buf=fp + fn + ".csv", sep=',', na_rep='NaN', index=False)
    Readme = "Readme for " + fn + ".csv"
    Readme += "\n\nDerived from " + fn_meteo
    Readme += "\nKersti Haahti, Luke " + str(datetime.datetime.now().date()) + readme
    outF = open(fp + fn + "_readme.txt", "w")
    print >>outF, Readme
    outF.close()