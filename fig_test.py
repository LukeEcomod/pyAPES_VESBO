# -*- coding: utf-8 -*-
"""
Created on Thu Feb 01 10:49:48 2018

@author: L1656
"""
import numpy as np
from datetime import datetime
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
pal = sns.color_palette("hls", 5)

yearly_cum = np.concatenate(([np.cumsum(res['Efloor'])],
                            [np.cumsum(res['Transp'])],
                            [np.cumsum(res['CanEvap'])],
                            [np.cumsum(res['Drain'])],
                            [np.cumsum(res['Roff'])]), axis=0)
for k in range (Forc.index.year[0], Forc.index.year[-1]):
    ix = np.where(Forc.index.year > k)
    for i in range(0, len(yearly_cum)):
        yearly_cum[i,ix] = yearly_cum[i,ix] - yearly_cum[i,ix[0][0]]

# Read weir
fp = "C:/Users/L1656/Documents/Git_repos/CCFPeat/forcing/Lettosuo_weir.csv"
dat = pd.read_csv(fp, sep=',', header='infer', na_values=-999)
# parse first columns to datetime
N = len(dat)
t = []
for k in range(N):
    t.append(
        datetime(dat['yyyy'].iloc[k],
                 dat['mo'].iloc[k],
                 dat['dd'].iloc[k],
                 dat['hh'].iloc[k],
                 dat['mm'].iloc[k]))
dat.index = t
weir = dat[['runf']] * 3600 / 1000  # m/s

# Read gwl
fp = "C:/Users/L1656/Documents/Git_repos/CCFPeat/forcing/Lettosuo_gwl.csv"
dat = pd.read_csv(fp, sep=',', header='infer', na_values=-999)
# parse first columns to datetime
N = len(dat)
t = []
for k in range(N):
    t.append(
        datetime(dat['yyyy'].iloc[k],
                 dat['mo'].iloc[k],
                 dat['dd'].iloc[k],
                 dat['hh'].iloc[k],
                 dat['mm'].iloc[k]))
dat.index = t
gwl_meas = dat[['WT_E','WT_N','WT_W','WT_S']]

plt.clf()
plt.figure(1)
plt.subplot(5,1,1)
plt.plot(Forc.index, np.array(res['Prec']) - 40,color=pal[3], linewidth=0.5, label='Precipitation')
plt.plot(Forc.index, Forc['Tair'],color=pal[0], linewidth=0.5, label='Air temperature')
plt.ylabel('(mm) or (degC)')
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.legend(loc=2, fontsize=8)
plt.subplot(5,1,2)
plt.stackplot(Forc.index, yearly_cum,labels=['Efloor','Transp','IntEvap', 'Drainage', 'Surface runoff'], colors=pal)
plt.ylabel('(mm)')
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.legend(loc=2, fontsize=8)
plt.subplot(5,1,3)
plt.plot(Forc.index, res['Transp'], color=pal[1], linewidth=0.5, label='Transpiration')
plt.plot(Forc.index, res['Efloor'], color=pal[0], linewidth=0.5, label='Efloor')
plt.plot(Forc.index, np.add(res['Drain'], res['Roff']) / dt0 * 3600, 'k', linewidth=1, label='Total runoff')
plt.plot(weir.index, weir['runf'],':k', linewidth=1, label='Runoff (meas)')
#plt.ylim([0.0, 0.5])
plt.ylabel('(mm h$^{-1}$)')
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.subplot(5,1,4)
plt.plot(gwl_meas.index, gwl_meas['WT_E'],color=pal[2], linewidth=1)
plt.plot(gwl_meas.index, gwl_meas['WT_N'],color=pal[2], linewidth=1)
plt.plot(gwl_meas.index, gwl_meas['WT_W'],color=pal[2], linewidth=1)
plt.plot(gwl_meas.index, gwl_meas['WT_S'],color=pal[2], linewidth=1)
plt.plot(Forc.index, res['gwl'],color=pal[3],label='gwl')
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])
#plt.ylim([-1.5, 0.0])
plt.ylabel('(m)')
plt.subplot(5,1,5)
plt.plot(Forc.index, Forc['SnowD'] / 1000, 'gray', linewidth=1, label='Snow depth (meas)')
plt.plot(Forc.index, res['SWE'], 'k', label='SWE')
plt.ylabel('(cm) or (mm)')
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])
