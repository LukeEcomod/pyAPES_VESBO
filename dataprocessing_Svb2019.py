# -*- coding: utf-8 -*-
"""
Created on Thu May 14 09:02:23 2020

@author: 03110850
"""


import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pyAPES_utilities.timeseries_tools import fill_gaps
from pyAPES_utilities.plotting import plot_columns
from pyAPES_utilities.dataprocessing_scripts import create_forcingfile, save_df_to_csv

dat = pd.read_csv(r'O:\Projects\Antoine_SLU\Environmental_data_2019.csv',
                 sep=';', header='infer')

dat.index = pd.to_datetime({'year': dat['year'],
                            'month': dat['month'],
                            'day': dat['day'],
                            'hour': dat['hour'],
                            'minute': dat['minute']})
dat.index = dat.index.round('30T')

# soil temperature
plt.figure()
titles = ['Tsoil at 5 cm', 'Tsoil at 10 cm', 'Tsoil at 15 cm', 'Tsoil at 30 cm', 'Tsoil at 50 cm']
ax=[]
for i in range(5):
    if i > 0:
        ax.append(plt.subplot(5,1,i+1, sharex=ax[0], sharey=ax[0]))
    else:
        ax.append(plt.subplot(5,1,i+1))
    dat.ix[:,8+i:28:5].plot(ax=ax[i])
    dat['TS_avg_' + str(i+1) +'_1']=dat.ix[:,8+i:28:5].mean(axis=1)
    dat[['TS_avg_' + str(i+1) +'_1']].plot(ax=ax[i],color='k')
    plt.title(titles[i])

# soil volumetric water content
plt.figure()
titles = ['Wliq at 5 cm (vert)', 'Wliq at 5 cm (hor)', 'Wliq at 10 cm', 'Wliq at 30 cm', 'Wliq at 50 cm']
depths=[5,]
ax=[]
for i in range(5):
    if i > 0:
        ax.append(plt.subplot(5,1,i+1, sharex=ax[0], sharey=ax[0]))
    else:
        ax.append(plt.subplot(5,1,i+1))
    dat.ix[:,32+i:32+20:5].plot(ax=ax[i])
    dat['SWC_avg_' + str(i+1) +'_1']=dat.ix[:,32+i:32+20:5].mean(axis=1)
    dat[['SWC_avg_' + str(i+1) +'_1']].plot(ax=ax[i],color='k')
    dat[['VWC_%' + str(i+1)]] = 100*dat[['VWC_' + str(i+1)]]
    dat[['VWC_%' + str(i+1)]].plot(ax=ax[i],color='grey')
    plt.title(titles[i])

dat2 = pd.read_csv(r'O:\Projects\Antoine_SLU\EddyCovarianceData_30min_FinalAnalysis_201901_201912_decouplfil.csv',
                 sep=';', header='infer')

dat2.index = pd.to_datetime({'year': dat2['year'],
                            'month': dat2['month'],
                            'day': dat2['day'],
                            'hour': dat2['hour'],
                            'minute': dat2['minute']})
dat2.index = dat2.index.round('30T')

dat = dat[['SWC_avg_1_1','TS_avg_1_1']].merge(dat2, how='outer', left_index=True, right_index=True)

dat3 = pd.read_csv(r'O:\Projects\Antoine_SLU\ICOS_Carbon_portal\Svb_meteo2019.csv',
                 sep=';', header='infer')

dat3.index = pd.to_datetime({'year': dat3['year'],
                            'month': dat3['month'],
                            'day': dat3['day'],
                            'hour': dat3['hour'],
                            'minute': dat3['minute']})
dat3.index = dat3.index.round('30T')

dat = dat.merge(dat3, how='outer', left_index=True, right_index=True)

# dat.columns

dat['PPFD_IN_1_2_2'] = dat['PPFD_IN_1_2_1'] / 4.56
dat['Pa_1_1_1'] = np.where(dat['Pa_1_1_1']<200.0, np.nan, dat['Pa_1_1_1'])

# plot_columns(dat[['Swin_1_1_1','Swin_1_2_1','PPFD_IN_1_2_1']])
# -> fraction of PAR is 0.49

# meteo
dat[['Pa_1_1_1','Ta_1_1_1','RH_1_1_1','P_1_1_1']].plot(subplots=True)
dat[['wndspd_f_ms-1', 'wndspd_ms-1']].plot()
dat[['Swin_1_1_1','Swin_1_2_1','PPFD_IN_1_2_2']].plot()
dat[['Lwin_1_2_1','Lwout_1_2_1']].plot()
dat[['CO2_MixingRatio_umolmol-1']].plot()
dat[['H2O_MixingRatio_mmolmol-1']].plot()  # unit?!?!

#%% meteo file & forcing

plot=False

frames = []
readme = ""

# top layer soil temperature from average of four pits
df, info = fill_gaps(dat[['TS_avg_1_1']],
                     'Tsoil', 'Soil temperature at 5cm [degC]', fill_nan='linear',
                     plot=plot)
frames.append(df)
readme += info

# top layer soil moisture from SE-Svb_eco: SWC_2_2_1
df, info = fill_gaps(dat[['SWC_avg_1_1']],
                     'Wliq', 'Soil volumetric moisture content at 5cm [%]', fill_nan='linear',
                     plot=plot)
frames.append(df)
readme += info

# Air temperature
df, info = fill_gaps(dat[['Ta_1_1_1']],
                     'Tair', 'Air temperature [degC]', fill_nan='linear',
                     plot=plot)
frames.append(df)
readme += info

# Relative humidity
df, info = fill_gaps(dat[['RH_1_1_1']],
                     'RH', 'Relative humidity [%]', fill_nan='linear',
                     plot=plot)
frames.append(df)
readme += info

# Pressure
df, info = fill_gaps(dat[['Pa_1_1_1']],
                     'P', 'Ambient pressure [hPa]', fill_nan='linear',
                     plot=plot)
frames.append(df)
readme += info

# Incoming shortwave radiation
df, info = fill_gaps(dat[['Swin_1_1_1','Swin_1_2_1']],
                     'Rg', 'Global radiation i.e. incoming shortwave radiation [W/m2]', fill_nan='linear',
                     plot=plot)
frames.append(df)
readme += info

# Downwelling longwave radiation
df, info = fill_gaps(dat[['Lwin_1_2_1']],
                     'LWin', 'Downwelling long wave radiation [W/m2]', fill_nan='linear',
                     plot=plot)
frames.append(df)
readme += info

# Wind speed
df, info = fill_gaps(dat[['wndspd_f_ms-1']],
                      'U', 'Wind speed [m/s]', fill_nan='linear',
                      plot=plot)
frames.append(df)
readme += info

# Friction velocity
# plot_columns(dat[['wndspd_f_ms-1','Ustar_ms-1']])
dat['Ustar = 0.19 * U'] = 0.19 * df['U']
df, info = fill_gaps(dat[['Ustar_ms-1',
                           'Ustar = 0.19 * U']],
                    'Ustar', 'Friction velocity [m/s]', fill_nan='linear',
                      plot=plot)
frames.append(df)
readme += info

# Precipitation
df, info = fill_gaps(dat[['P_1_1_1']],
                    'Prec', 'Precipitation [mm/30min]', fill_nan=0.0,
                     plot=plot)
frames.append(df)
readme += info

# Mixing ration of CO2
df, info = fill_gaps(dat[['CO2_MixingRatio_umolmol-1']],
                    'CO2', 'Mixing ration of CO2 [ppm]', fill_nan='linear',
                     plot=plot)
frames.append(df)
readme += info

Svb_final=pd.concat(frames, axis=1)
Svb_final[['Tair', 'Rg', 'LWin', 'RH', 'P', 'U', 'Ustar', 'Prec', 'CO2', 'Tsoil', 'Wliq']].plot(subplots=True,kind='line',
         title=['Air temperature [degC]',
                'Incoming shortwave radiation [W/m2]',
                'Downwelling long wave radiation [W/m2]',
                'Relative humidity [%]',
                'Ambient pressure [hPa]',
                'Wind speed [m/s]',
                'Friction velocity [m/s]',
                'Precipitation [mm/30min]',
                'Mixing ration of CO2 [ppm]',
                'Soil temperature at 5cm [degC]',
                'Soil volumetric moisture content at 5cm [%]'],legend=False)
plt.tight_layout()
plt.xlim('1.1.2019','1.1.2020')

Svb_final=Svb_final[(Svb_final.index >= '1.1.2019') & (Svb_final.index <= '1.1.2020')]

save_df_to_csv(Svb_final, 'Svartberget_meteo_2019', readme=readme,
                fp='forcing/Svartberget/', timezone=+1)

create_forcingfile('forcing/Svartberget/Svartberget_meteo_2019.csv', 'Svartberget_forcing_2019',
                    'forcing/Svartberget/', lat=64.26, lon=19.77, P_unit=1e2,timezone=+1.0, fpar=0.49) # [hPa]

# %% flux data

dat[['H_f_Wm-2','H_orig_Wm-2']].plot()
dat[['LE_f_Wm-2','LE_orig_Wm-2']].plot()
dat[['NEE_f_umolm-2s-1','NEE_orig_umolm-2s-1']].plot()
dat['GPP_f_umolm-2s-1'] = -dat['GPP_f_umolm-2s-1']
dat['GPP_umolm-2s-1'] = np.where(np.isfinite(dat['NEE_orig_umolm-2s-1']), dat['GPP_f_umolm-2s-1'], np.nan)
dat[['GPP_f_umolm-2s-1','GPP_umolm-2s-1', 'Reco_umolm-2s-1']].plot()

frames = []
readme = ""

# LWnet
df, info = fill_gaps(dat[['Lwnet_1_2_1']],
                    'LWnet', 'Net longwave radiation [W m-2]', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info

# SWnet
df, info = fill_gaps(dat[['Swnet_1_2_1']],
                    'SWnet', 'Net shortwave radiation [W m-2]', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info

# Rnet
df, info = fill_gaps(dat[['NetRad_1_2_1']],
                    'Rnet', 'Net radiation [W m-2]', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info

# sensible heat
df, info = fill_gaps(dat[['H_orig_Wm-2']],
                    'SH', 'Sensible heat flux [W m-2]', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info
df, info = fill_gaps(dat[['H_f_Wm-2']],
                    'SH_gapfilled', 'Sensible heat flux [W m-2], gapfilled', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info

# latent heat
df, info = fill_gaps(dat[['LE_orig_Wm-2']],
                    'LE', 'Latent heat flux [W m-2]', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info
df, info = fill_gaps(dat[['LE_f_Wm-2']],
                    'LE_gapfilled', 'Latent heat flux [W m-2], gapfilled', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info

# NEE
df, info = fill_gaps(dat[['NEE_orig_umolm-2s-1']],
                    'NEE', 'Net ecosystem exchange [umol m-2 s-1]', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info
df, info = fill_gaps(dat[['NEE_f_umolm-2s-1']],
                    'NEE_gapfilled', 'Net ecosystem exchange [umol m-2 s-1], gapfilled', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info

# GPP
df, info = fill_gaps(dat[['GPP_umolm-2s-1']],
                    'GPP', 'Gross primary production [umol m-2 s-1]', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info
df, info = fill_gaps(dat[['GPP_f_umolm-2s-1']],
                    'GPP_gapfilled', 'Gross primary production [umol m-2 s-1], gapfilled', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info

# Reco
df, info = fill_gaps(dat[['Reco_umolm-2s-1']],
                    'Reco', 'Ecosystem respiration [umol m-2 s-1]', fill_nan=np.nan,
                    plot=plot)
frames.append(df)
readme += info

Svb_EC=pd.concat(frames, axis=1)

Svb_EC = Svb_EC[(Svb_EC.index >= '1.1.2019') & (Svb_EC.index <= '1.1.2020')]

save_df_to_csv(Svb_EC, 'Svartberget_EC_2019', readme=readme,
               fp='forcing/Svartberget/', timezone=+1)