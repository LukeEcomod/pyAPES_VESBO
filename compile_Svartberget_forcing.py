# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 10:10:39 2019

@author: L1656
"""
# %% collect short time data from directories to larger files

#from tools.dataprocessing_scripts import gather_data
#
#directory = 'H:/Muut projektit/Natalia/Svartberget_data/'
#
#foldernames = [
#        'SE-Svb_eco/',
#        'SE-Svb_fluxes/',
#        'SE-Svb_meteo/',
#        'SE-Svb_profile/',
#        'SE-Svb_T-profile/',
#        'SE-Deg_meteo/',
#        'SE-Deg_fluxes/'
#        ]
#
#for fn in foldernames:
#    print(fn)
#    gather_data(directory + fn, sort=1)
#    
# %% read data from files

from tools.dataprocessing_scripts import read_Svb_data

fp = ["H:/Muut projektit/Natalia/Svartberget_data/SE-Svb_eco/concat.csv",
      "H:/Muut projektit/Natalia/Svartberget_data/SE-Svb_fluxes/concat.csv",
      "H:/Muut projektit/Natalia/Svartberget_data/SE-Deg_fluxes/concat.csv",
      "H:/Muut projektit/Natalia/Svartberget_data/SE-Svb_meteo/concat.csv",
      "H:/Muut projektit/Natalia/Svartberget_data/SE-Deg_meteo/concat.csv",
      "H:/Muut projektit/Natalia/Svartberget_data/SE-Svb_profile/concat.csv",
      "H:/Muut projektit/Natalia/Svartberget_data/SE-Svb_T-profile/concat.csv"]

Svb_eco = read_Svb_data([fp[0]])
Svb_fluxes= read_Svb_data([fp[1],fp[2]])
Svb_meteo = read_Svb_data([fp[3],fp[4]])

# %% processing data

from matplotlib import pyplot as plt
import numpy as np
from tools.timeseries_tools import fill_gaps

# soil temperature
plt.figure()
titles = ['Tsoil at 5 cm', 'Tsoil at 10 cm', 'Tsoil at 15 cm', 'Tsoil at 30 cm', 'Tsoil at 50 cm']
ax=[]
for i in range(5):
    if i > 0:
        ax.append(plt.subplot(5,1,i+1, sharex=ax[0], sharey=ax[0]))
    else:
        ax.append(plt.subplot(5,1,i+1))
    Svb_eco.ix[:,i:20:5].plot(ax=ax[i])
    Svb_eco['SE-Svb_meteo: TS_avg_' + str(i) +'_1']=Svb_eco.ix[:,i:20:5].mean(axis=1)
    plt.title(titles[i])

# soil volumetric water content
plt.figure()
titles = ['Wliq at 5 cm (vert)', 'Wliq at 5 cm (hor)', 'Wliq at 10 cm', 'Wliq at 30 cm', 'Wliq at 50 cm']
ax=[]
for i in range(5):
    if i > 0:
        ax.append(plt.subplot(5,1,i+1, sharex=ax[0], sharey=ax[0]))
    else:
        ax.append(plt.subplot(5,1,i+1))
    Svb_eco.ix[:,24+i:24+20:5].plot(ax=ax[i])
    Svb_eco['SE-Svb_meteo: SWC_avg_' + str(i) +'_1']=Svb_eco.ix[:,24+i:24+20:5].mean(axis=1)
    plt.title(titles[i])

# top layer soil moisture from SE-Svb_eco: SWC_2_2_1, quality check
Svb_eco['Wliq']=np.where((Svb_eco['SE-Svb_eco: SWC_2_2_1'] < 5), np.nan,Svb_eco['SE-Svb_eco: SWC_2_2_1'])
Svb_eco['Wliq'][:-1]=np.where(np.isnan(Svb_eco['Wliq'][1:]), np.nan,Svb_eco['Wliq'][:-1])
Svb_eco['Wliq'][1:]=np.where(np.isnan(Svb_eco['Wliq'][:-1]), np.nan,Svb_eco['Wliq'][1:])
Svb_eco['Wliq']= Svb_eco['Wliq'].interpolate()
Svb_eco['Wliq_12h'] = Svb_eco['Wliq'].rolling(24, 1).mean()
Svb_eco['Wliq']=np.where((Svb_eco['SE-Svb_eco: SWC_2_2_1'] + 1 < Svb_eco['Wliq_12h']), np.nan,Svb_eco['SE-Svb_eco: SWC_2_2_1'])
Svb_eco['SE-Svb_eco: SWC_2_2_1_screened'] = Svb_eco['Wliq']

# meteo
from tools.plotting import plot_columns
import pandas as pd

# Variables screened using Degerö air pressure in comparison to Svb
variables=['SE-Svb_meteo: Pa_1_1_1',
           'SE-Svb_meteo: Ta_1_1_1',
           'SE-Svb_meteo: RH_1_1_1',
           'SE-Svb_meteo: Swin_1_1_1',
           'SE-Svb_meteo: Swin_1_2_1',
           'SE-Svb_meteo: Lwin_1_2_1',
           'SE-Svb_meteo: Lwout_1_2_1',
           'SE-Svb_meteo: PPFD_IN_1_2_1',
           'SE-Svb_meteo: PPFD_DIR_1_1_1',
           'SE-Svb_meteo: PPFD_DIFF_1_1_1']
# screen out unrealistic Pressure values
Svb_meteo['SE-Svb_meteo: Pa_1_1_1']=np.where(Svb_meteo['SE-Svb_meteo: Pa_1_1_1'] < 900,
                                             np.nan, Svb_meteo['SE-Svb_meteo: Pa_1_1_1'])
for var in variables[3:]:
    Svb_meteo[var]=np.where(Svb_meteo[var] < 0.0,
                            np.nan, Svb_meteo[var])
for var in variables:
    Svb_meteo[var + '_screened']=np.where(abs(
        Svb_meteo['SE-Svb_meteo: Pa_1_1_1']-Svb_meteo['SE-Deg_meteo: Pa_1_1_1']) > 5,
        np.nan, Svb_meteo[var])

# radiation
#plt.figure()
#ax1 = plt.subplot(7,1,1)
#Svb_meteo[[
#           'SE-Svb_meteo: Pa_1_1_1',
#           'SE-Deg_meteo: Pa_1_1_1']].plot(ax=ax1)
#ax = plt.subplot(7,1,2, sharex=ax1)
#Svb_meteo[[
#           'SE-Svb_meteo: Ta_1_1_1',
#           'SE-Deg_meteo: Ta_1_1_1']].plot(ax=ax)
#ax = plt.subplot(7,1,3, sharex=ax1)
#Svb_meteo[[
#           'SE-Svb_meteo: RH_1_1_1',
#           'SE-Deg_meteo: RH_1_1_1']].plot(ax=ax)
#ax = plt.subplot(7,1,4, sharex=ax1)
#Svb_meteo[[
#           'SE-Svb_meteo: Swin_1_1_1',
#           'SE-Svb_meteo: Swin_1_2_1',
#           'SE-Deg_meteo: Swin_1_1_1',
#           'SE-Deg_meteo: Swin_1_2_1']].plot(ax=ax)
#ax = plt.subplot(7,1,5, sharex=ax1)
#Svb_meteo[[
#           'SE-Svb_meteo: Lwin_1_2_1',
#           'SE-Svb_meteo: Lwout_1_2_1',
#           'SE-Deg_meteo: Lwin_1_2_1',
#           'SE-Deg_meteo: Lwout_1_2_1']].plot(ax=ax)
#ax = plt.subplot(7,1,6, sharex=ax1)
#Svb_meteo[[
#           'SE-Svb_meteo: PPFD_IN_1_2_1',
#           'SE-Svb_meteo: PPFD_DIR_1_1_1',
#           'SE-Svb_meteo: PPFD_DIFF_1_1_1',
#           'SE-Deg_meteo: PPFD_IN_1_2_1',
#           'SE-Deg_meteo: PPFD_DIR_1_1_1',
#           'SE-Deg_meteo: PPFD_DIFF_1_1_1']].plot(ax=ax)
#ax = plt.subplot(7,1,7, sharex=ax1)
#Svb_meteo[[
#           'SE-Svb_meteo: P_1_1_1',
#           'SE-Deg_meteo: P_1_1_1']].plot(ax=ax)
#
#plot_columns(Svb_meteo[['SE-Svb_meteo: Pa_1_1_1_screened',
#                        'SE-Deg_meteo: Pa_1_1_1']],plot_timeseries=False)
#plot_columns(Svb_meteo[['SE-Deg_meteo: Ta_1_1_1',
#                        'SE-Svb_meteo: Ta_1_1_1_screened']],plot_timeseries=False)
#plot_columns(Svb_meteo[['SE-Deg_meteo: RH_1_1_1',
#                        'SE-Svb_meteo: RH_1_1_1_screened']],plot_timeseries=False)
#plot_columns(Svb_meteo[['SE-Deg_meteo: Lwin_1_2_1',
#                        'SE-Deg_meteo: Lwout_1_2_1',
#                        'SE-Svb_meteo: Lwin_1_2_1_screened',
#                        'SE-Svb_meteo: Lwout_1_2_1_screened']],plot_timeseries=False)
#plot_columns(Svb_meteo[['SE-Svb_meteo: Swin_1_1_1_screened',
#                        'SE-Svb_meteo: Swin_1_2_1_screened',
#                        'SE-Deg_meteo: Swin_1_1_1',
#                        'SE-Deg_meteo: Swin_1_2_1']],plot_timeseries=False)
#plot_columns(Svb_meteo[['SE-Svb_meteo: PPFD_IN_1_2_1_screened',
#                        'SE-Svb_meteo: PPFD_DIR_1_1_1_screened',
#                        'SE-Svb_meteo: PPFD_DIFF_1_1_1_screened',
#                        'SE-Deg_meteo: PPFD_IN_1_2_1',
#                        'SE-Deg_meteo: PPFD_DIR_1_1_1',
#                        'SE-Deg_meteo: PPFD_DIFF_1_1_1']],plot_timeseries=False)
#plot_columns(Svb_meteo[['SE-Svb_meteo: Swin_1_1_1_screened',
#                        'SE-Svb_meteo: Swin_1_2_1_screened',
#                        'SE-Svb_meteo: PPFD_IN_1_2_1_screened',
#                        'SE-Svb_meteo: PPFD_DIR_1_1_1_screened']],plot_timeseries=False)
# precipitation
Svb_meteo[['SE-Svb_meteo: P_1_1_1','SE-Deg_meteo: P_1_1_1']] = Svb_meteo[
        ['SE-Svb_meteo: P_1_1_1','SE-Deg_meteo: P_1_1_1']].fillna(-999)

Prec_daily = Svb_meteo[['SE-Svb_meteo: P_1_1_1',
                        'SE-Deg_meteo: P_1_1_1']].resample('D').mean()
Prec_daily['SE-Svb_meteo: P_1_1_1'] = np.where(
        Prec_daily['SE-Svb_meteo: P_1_1_1'] < 0.0, np.nan, Prec_daily['SE-Svb_meteo: P_1_1_1'])
Prec_daily['SE-Deg_meteo: P_1_1_1'] = np.where(
        Prec_daily['SE-Deg_meteo: P_1_1_1'] < 0.0, np.nan, Prec_daily['SE-Deg_meteo: P_1_1_1'])

Prec_daily = Prec_daily.rename(columns={'SE-Svb_meteo: P_1_1_1':'SE-Svb_meteo: P_daily', 
                                        'SE-Deg_meteo: P_1_1_1':'SE-Deg_meteo: P_daily'})

Svb_meteo['SE-Svb_meteo: P_1_1_1'] = np.where(
        Svb_meteo['SE-Svb_meteo: P_1_1_1'] < 0.0, np.nan, Svb_meteo['SE-Svb_meteo: P_1_1_1'])
Svb_meteo['SE-Deg_meteo: P_1_1_1'] = np.where(
        Svb_meteo['SE-Deg_meteo: P_1_1_1'] < 0.0, np.nan, Svb_meteo['SE-Deg_meteo: P_1_1_1'])
    
df = pd.read_csv('H:/Muut projektit/Natalia/Svartberget_data/daily_prec.csv',
                 sep=',', header='infer')
df.index = pd.to_datetime(df.ix[:,0], dayfirst=True)
df.ix[:,0]=df.index
df = df.drop([df.columns[0]], axis=1)
df['Prec (mm)']=df['Prec (mm)'] / 48
Prec_daily = Prec_daily.merge(df, how='outer', left_index=True, right_index=True)

#plot_columns(Svb_meteo[['SE-Svb_meteo: P_1_1_1',
#                        'SE-Deg_meteo: P_1_1_1']])
#plot_columns(Prec_daily)

Prec_daily[['SE-Svb_meteo: P_daily','SE-Deg_meteo: P_daily','Prec (mm)']] = Prec_daily[
        ['SE-Svb_meteo: P_daily','SE-Deg_meteo: P_daily','Prec (mm)']].fillna(-999)
Svb_meteo = Svb_meteo.merge(Prec_daily, how='outer', left_index=True, right_index=True)
Svb_meteo['SE-Svb_meteo: P_daily'] = Svb_meteo['SE-Svb_meteo: P_daily'].fillna(method='ffill')
Svb_meteo['SE-Deg_meteo: P_daily'] = Svb_meteo['SE-Deg_meteo: P_daily'].fillna(method='ffill')
Svb_meteo['Prec (mm)'] = Svb_meteo['Prec (mm)'].fillna(method='ffill')
Svb_meteo['SE-Svb_meteo: P_daily'] = np.where(
        Svb_meteo['SE-Svb_meteo: P_daily'] < 0.0, np.nan, Svb_meteo['SE-Svb_meteo: P_daily'])
Svb_meteo['SE-Deg_meteo: P_daily'] = np.where(
        Svb_meteo['SE-Deg_meteo: P_daily'] < 0.0, np.nan, Svb_meteo['SE-Deg_meteo: P_daily'])
Svb_meteo['Prec (mm)'] = np.where(
        Svb_meteo['Prec (mm)'] < 0.0, np.nan, Svb_meteo['Prec (mm)'])
Svb_meteo['SE-Deg_meteo: P_corrected to match daily Prec'] = (Svb_meteo['SE-Deg_meteo: P_1_1_1'] *
         Svb_meteo['Prec (mm)'] / Svb_meteo['SE-Deg_meteo: P_daily'])

from tools.plotting import plot_timeseries_df

# Precipitation
df, _ = fill_gaps(Svb_meteo[['SE-Svb_meteo: P_1_1_1',
                             'SE-Deg_meteo: P_corrected to match daily Prec',
                             'Prec (mm)']],
                    'Prec', 'Precipitation  [mm/30min]', fill_nan=0.0, 
                     plot=True)

#plt.figure()
#plot_timeseries_df(Svb_meteo, ['SE-Svb_meteo: P_1_1_1',
#                               'SE-Deg_meteo: P_1_1_1',
#                               'SE-Svb_meteo: P_daily',
#                               'SE-Deg_meteo: P_daily',
#                               'Prec (mm)',
#                               'SE-Deg_meteo: P_corrected to match daily Prec'],limits=False,cum=True)
#plot_timeseries_df(df, ['Prec'],limits=False,cum=True)
#
#plt.figure()
#plot_timeseries_df(Svb_meteo, ['SE-Svb_meteo: P_1_1_1',
#                               'SE-Deg_meteo: P_1_1_1',
#                               'SE-Svb_meteo: P_daily',
#                               'SE-Deg_meteo: P_daily',
#                               'Prec (mm)',
#                               'SE-Deg_meteo: P_corrected to match daily Prec'],limits=False)
#plot_timeseries_df(df, ['Prec'],limits=False, linestyle=':')
#
#
#plot_columns(Svb_meteo[['SE-Svb_meteo: P_1_1_1',
#                        'SE-Deg_meteo: P_1_1_1',
#                        'SE-Svb_meteo: P_daily',
#                        'SE-Deg_meteo: P_daily',
#                        'Prec (mm)']])
# wind
variables = ['SE-Svb_fluxes: Ustar_1_1_1',
             'SE-Svb_fluxes: WS_1_1_1']
for var in variables:
    Svb_fluxes[var]=np.where(Svb_fluxes[var] < 0.0,
                            np.nan, Svb_fluxes[var])
Svb_fluxes[['SE-Svb_fluxes: Ustar_1_1_1',
            'SE-Svb_fluxes: WS_1_1_1',
            'SE-Deg_fluxes: Ustar_1_1_1',
            'SE-Deg_fluxes: WS_1_1_1']].plot()
#plot_columns(Svb_fluxes[['SE-Svb_fluxes: Ustar_1_1_1',
#                         'SE-Svb_fluxes: WS_1_1_1',
#                         'SE-Deg_fluxes: Ustar_1_1_1',
#                         'SE-Deg_fluxes: WS_1_1_1']],plot_timeseries=False)

# %% gapfill and collect data

from tools.iotools import save_df_to_csv
from tools.dataprocessing_scripts import create_forcingfile

frames = []
readme = ""

# top layer soil temperature from average of four pits
df, info = fill_gaps(Svb_eco[['SE-Svb_meteo: TS_avg_1_1']],
                     'Tsoil', 'Soil temperature at 5cm [degC]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# top layer soil moisture from SE-Svb_eco: SWC_2_2_1
df, info = fill_gaps(Svb_eco[['SE-Svb_eco: SWC_2_2_1_screened']],
                     'Wliq', 'Soil volumetric moisture content at 5cm [%]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# Air temperature
df, info = fill_gaps(Svb_meteo[['SE-Svb_meteo: Ta_1_1_1_screened',
                                'SE-Deg_meteo: Ta_1_1_1']],
                     'Tair', 'Air temperature [degC]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# Relative humidity
df, info = fill_gaps(Svb_meteo[['SE-Svb_meteo: RH_1_1_1_screened',
                                'SE-Deg_meteo: RH_1_1_1']],
                     'RH', 'Relative humidity [%]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# Pressure
df, info = fill_gaps(Svb_meteo[['SE-Svb_meteo: Pa_1_1_1_screened',
                                'SE-Deg_meteo: Pa_1_1_1']],
                     'P', 'Ambient pressure [hPa]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# Incoming shortwave radiation
df, info = fill_gaps(Svb_meteo[['SE-Svb_meteo: Swin_1_1_1_screened',
                                'SE-Svb_meteo: Swin_1_2_1_screened',
                                'SE-Deg_meteo: Swin_1_1_1',
                                'SE-Deg_meteo: Swin_1_2_1']],
                     'Rg', 'Global radiation i.e. incoming shortwave radiation [W/m2]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# Downwelling longwave radiation
df, info = fill_gaps(Svb_meteo[['SE-Svb_meteo: Lwin_1_2_1_screened',
                                'SE-Deg_meteo: Lwin_1_2_1']],
                     'LWin', 'Downwelling long wave radiation [W/m2]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# Upwelling longwave radiation
df, info = fill_gaps(Svb_meteo[['SE-Svb_meteo: Lwout_1_2_1_screened',
                                'SE-Deg_meteo: Lwout_1_2_1']],
                     'LWout', 'Upwelling long wave radiation [W/m2]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# Wind speed
df, info = fill_gaps(Svb_fluxes[['SE-Svb_fluxes: WS_1_1_1',
                                'SE-Deg_fluxes: WS_1_1_1']],
                     'U', 'Wind speed [m/s]', fill_nan='linear',
                     plot=True)
frames.append(df)
readme += info

# Friction velocity
Svb_fluxes['Ustar = 0.17 * U'] = 0.17 * df['U']
df, info = fill_gaps(Svb_fluxes[['SE-Svb_fluxes: Ustar_1_1_1',
                                'Ustar = 0.17 * U']],
                    'Ustar', 'Friction velocity [m/s]', fill_nan='linear', 
                     plot=True)
frames.append(df)
readme += info

# Precipitation
df, info = fill_gaps(Svb_meteo[['SE-Svb_meteo: P_1_1_1',
                             'SE-Deg_meteo: P_corrected to match daily Prec',
                             'Prec (mm)']],
                    'Prec', 'Precipitation  [mm/30min]', fill_nan=0.0, 
                     plot=True)
frames.append(df)
readme += info

Svb_final=pd.concat(frames, axis=1)
Svb_final[['Tair', 'Rg', 'LWin', 'U', 'RH', 'P', 'Prec', 'Tsoil', 'Wliq']].plot(subplots=True,kind='line',
         title=['Air temperature [degC]',
                'Incoming shortwave radiation [W/m2]',
                'Downwelling long wave radiation [W/m2]',
                'Wind speed [m/s]',
                'Relative humidity [%]',
                'Ambient pressure [hPa]',
                'Precipitation [mm/30min]',
                'Soil temperature at 5cm [degC]',
                'Soil volumetric moisture content at 5cm [%]'],legend=False)
plt.tight_layout()
plt.xlim('1.1.2014','1.1.2017')

save_df_to_csv(Svb_final, "Svarberget_meteo_2014_2016", readme=readme, fp="C:/Users/L1656/Documents/Git_repos/pyAPES/forcing/")

create_forcingfile("Svarberget_meteo_2014_2016", "Svarberget_forcing_2014_2016",
                   lat=64.26, lon=19.77, P_unit = 1e2) # [hPa]
