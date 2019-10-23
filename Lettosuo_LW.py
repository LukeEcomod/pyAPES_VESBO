# -*- coding: utf-8 -*-
"""
Created on Wed May 29 09:53:24 2019

@author: L1656
"""

import numpy as np
import pandas as pd
from os import listdir
import matplotlib.pyplot as plt

from pyAPES_utilities.timeseries_tools import fill_gaps
from pyAPES_utilities.plotting import plot_xy, plot_diurnal

#: machine epsilon
EPS = np.finfo(float).eps

prop_cycle = plt.rcParams['axes.prop_cycle']
default = prop_cycle.by_key()['color']

direc = "C:/Users/L1656/Documents/Git_repos/pyAPES_Kersti/"
meteo_file="Lettosuo_meteo_2010_2018"
lat=60.63
lon=23.95
timezone=+2.0

from canopy.radiation import solar_angles, compute_clouds_rad
from canopy.micromet import e_sat

fpar = 0.45

forc_fp = direc + "forcing/" + meteo_file +".csv"
dat = pd.read_csv(forc_fp, sep=',', header='infer', encoding = 'ISO-8859-1')

# set to dataframe index
dat.index = pd.to_datetime({'year': dat['yyyy'],
                            'month': dat['mo'],
                            'day': dat['dd'],
                            'hour': dat['hh'],
                            'minute': dat['mm']})

dt = (dat.index[1] - dat.index[0]).total_seconds()

# day of year
dat['doy'] = dat.index.dayofyear

# ambient H2O [mol/mol] from RH
esat, _ = e_sat(dat['Tair'])
dat['H2O'] = (dat['RH'] / 100.0) * esat / dat['P']

# zenith angle
jday = dat.index.dayofyear + dat.index.hour / 24.0 + dat.index.minute / 1440.0
# TEST (PERIOD START)
jday = dat.index.dayofyear + dat.index.hour / 24.0 + dat.index.minute / 1440.0 + dt / 2.0 / 86400.0
dat['Zen'], _, _, _, _, _ = solar_angles(lat, lon, jday, timezone=timezone)

b = 5.6697e-8  # Stefan-Boltzman constant (W m-2 K-4)

# Tervalammen suo
fp = "O:/Projects/Lettosuo/Forcing_data/TervisLW_2013-2018.csv"
LW_meas = pd.read_csv(fp, sep=',', header='infer', encoding = 'ISO-8859-1')
LW_meas.index = pd.to_datetime(LW_meas.ix[:,0], dayfirst=True)
dat=dat.merge(LW_meas, how='outer', left_index=True, right_index=True)
dat['emi_sky_meas'] = dat['LW sky (W m-2)'] / b /(dat['Tair'] + 273.15)**4
dat['emi0'] = 1.24 * (dat['H2O'].values * dat['P'].values / 100 /(dat['Tair'].values + 273.15))**(1./7.)
dat['fcloud_meas'] = (dat['emi_sky_meas'] - dat['emi0']) / (0.84 - 0.84 * dat['emi0'])

# %%
from canopy.radiation import compute_clouds_rad
# radiation components
f_cloud, f_diff, emi_sky = compute_clouds_rad(dat['doy'].values,
                                              dat['Zen'].values,
                                              dat['Rg'].values,
                                              dat['H2O'].values * dat['P'].values,
                                              dat['Tair'].values)

# solar constant at top of atm.
So = 1367
# clear sky Global radiation at surface
dat['Qclear'] = np.maximum(0.0,
                    (So * (1.0 + 0.033 * np.cos(2.0 * np.pi * (np.minimum(dat['doy'].values, 365) - 10) / 365)) * np.cos(dat['Zen'].values)))

tau_atm = dat['Rg'].rolling(4,1).sum() / (dat['Qclear'].rolling(4,1).sum() + EPS)
dat['tau_atm']=tau_atm
# cloud cover fraction
dat['f_cloud2'] = 1.0 - (tau_atm -0.2) / (0.7 - 0.2)
dat['f_cloud2'][dat['Qclear'] < 10] = np.nan

dat['Qclear_12h'] = dat['Qclear'].resample('12H').sum()
dat['Qclear_12h'] = dat['Qclear_12h'].fillna(method='ffill')
dat['Rg_12h'] = dat['Rg'].resample('12H').sum()
dat['Rg_12h'] = dat['Rg_12h'].fillna(method='ffill')

tau_atm = dat['Rg_12h'] / (dat['Qclear_12h'] + EPS)
dat['f_cloud3'] = 1.0 - (tau_atm -0.2) / (0.7 - 0.2)

dat['f_cloud4'] = np.where((dat.index.hour > 12) & (dat['f_cloud3'] < 0.2), 0.0, dat['f_cloud2'])
dat['f_cloud4'] = dat['f_cloud4'].interpolate()
dat['f_cloud4'][dat['f_cloud4'] < 0.0] = 0.0
dat['f_cloud4'][dat['f_cloud4'] > 1.0] = 1.0

dat['emi_sky2'] = (1 - 0.84 * dat['f_cloud4']) * dat['emi0'] + 0.84 * dat['f_cloud4']
dat['LWin2'] = dat['emi_sky2'] * b *(dat['Tair'] + 273.15)**4 # Wm-2 downwelling LW

dat['emi_sky'] = emi_sky
dat['f_cloud'] = f_cloud

# estimated long wave budget
dat['LWin'] = emi_sky * b *(dat['Tair'] + 273.15)**4 # Wm-2 downwelling LW

plt.figure()
ax = plt.subplot(3,1,1)
dat[['LWin', 'LWin2', 'LW sky (W m-2)']].plot(ax=ax)
ax = plt.subplot(3,1,2, sharex=ax)
dat[['emi_sky', 'emi_sky2', 'emi_sky_meas']].plot(ax=ax)
ax = plt.subplot(3,1,3, sharex=ax)
dat[['f_cloud', 'f_cloud4', 'fcloud_meas', 'f_cloud3']].plot(ax=ax)

months = dat.index.month
ix = np.where((months >= 4) & (months <= 10) & np.isfinite(dat['LW sky (W m-2)']))[0]
plt.figure()
plt.subplot(3,2,1)
plot_diurnal(dat['LWin'][ix], color=default[0], legend=False)
plot_diurnal(dat['LWin2'][ix], color=default[1], legend=False)
plot_diurnal(dat['LW sky (W m-2)'][ix], color=default[2], legend=False)
plt.subplot(3,2,2)
plot_diurnal(dat['emi_sky'][ix], color=default[0], legend=False)
plot_diurnal(dat['emi_sky2'][ix], color=default[1], legend=False)
plot_diurnal(dat['emi_sky_meas'][ix], color=default[2], legend=False)
plt.subplot(3,2,3)
plot_xy(dat['LW sky (W m-2)'][ix], dat['LWin'][ix])
plt.subplot(3,2,5)
plot_xy(dat['LW sky (W m-2)'][ix], dat['LWin2'][ix],color=default[1])
plt.subplot(3,2,4)
plot_xy(dat['emi_sky_meas'][ix], dat['emi_sky'][ix])
plt.subplot(3,2,6)
plot_xy(dat['emi_sky_meas'][ix], dat['emi_sky2'][ix],color=default[1])