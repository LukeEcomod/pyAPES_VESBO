# -*- coding: utf-8 -*-
"""
Created on Thu Feb 01 10:49:48 2018

@author: L1656
"""
import numpy as np
from matplotlib import pyplot as plt
from forcing.forc_utils import read_forcing
from parameters.general_parameters import gpara
import seaborn as sns

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
weir = read_forcing("Lettosuo_weir.csv",
                       gpara['start_time'],
                       gpara['end_time'],
                       cols=['runf'])
weir.loc[:,'runf'] = weir.loc[:,'runf'] * 3600  # mm/h

# Read gwl
gwl_meas = read_forcing("Lettosuo_gwl.csv",
                           gpara['start_time'],
                           gpara['end_time'],
                           cols=['WT_E','WT_N','WT_W','WT_S'])

# Read ET
ET_hyde = read_forcing(gpara['forc_filename'],
                            gpara['start_time'],
                            gpara['end_time'],
                            cols=['ET'])
ET_hyde.loc[:,'ET'] = ET_hyde.loc[:,'ET'] * 0.0324  # mm/30min
yearly_cumET = np.cumsum(ET_hyde['ET'])
for k in range (Forc.index.year[0], Forc.index.year[-1]):
    ix = np.where(Forc.index.year > k)
    yearly_cumET[ix[0]] = yearly_cumET[ix[0]] - yearly_cumET[ix[0][0]]
    yearly_cumET[ix[0][0]-1] = float('nan')

# Read snow
snow = read_forcing("FMI_jokioinen.csv",
                            gpara['start_time'],
                            gpara['end_time'],
                            cols=['SnowD'])

plt.clf()
plt.figure(1)
plt.subplot(5,1,1)
plt.plot(Forc.index, 1000 * np.array(res['Prec']) / dt0 * 3600 - 40, color=pal[3], linewidth=0.5, label='Precipitation')
plt.plot(Forc.index, Forc['Tair'],color=pal[0], linewidth=0.5, label='Air temperature')
plt.ylabel('(mm h$^{-1}$) or (degC)')
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.legend(loc=2, fontsize=8)
plt.subplot(5,1,2)
plt.stackplot(Forc.index, 1000 * yearly_cum,labels=['Efloor','Transp','IntEvap', 'Drainage', 'Surface runoff'], colors=pal)
plt.plot(Forc.index, yearly_cumET, 'k', linewidth=1, label='ET_hyde')
plt.ylabel('(mm)')
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.legend(loc=2, fontsize=8)
plt.subplot(5,1,3)
plt.plot(Forc.index, 1000*np.array(res['Transp']) / dt0 * 3600, color=pal[1], linewidth=0.5, label='Transpiration')
plt.plot(Forc.index, 1000*np.array(res['Efloor']) / dt0 * 3600, color=pal[0], linewidth=0.5, label='Efloor')
plt.plot(Forc.index, 1000 * np.add(res['Drain'], res['Roff']) / dt0 * 3600, 'k', linewidth=1, label='Total runoff')
plt.plot(weir.index, weir['runf'],':k', linewidth=1, label='Runoff (meas)')
plt.ylim([0.0, 0.5])
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
plt.ylim([-1, 0.0])
plt.ylabel('(m)')
plt.subplot(5,1,5)
plt.plot(snow.index, snow['SnowD'], 'gray', linewidth=1, label='Snow depth (meas)')
plt.plot(Forc.index, 1000*np.array(res['SWE']), 'k', label='SWE')
plt.ylabel('(cm) or (mm)')
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])

plt.figure()
plt.plot(Forc.index, res['LAI'], 'r', label='LAI')
plt.plot(Forc.index, res['Phenof'], 'k', label='Pheno state')
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])