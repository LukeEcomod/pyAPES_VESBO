# -*- coding: utf-8 -*-
"""
Created on Thu Feb 01 10:49:48 2018

@author: L1656
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from forcing.forc_utils import read_forcing, diurnal_cycle
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
ET_hyde = read_forcing("Hyde_data_2010_2016.csv",
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

if cmodel.Switch_MLM:
    Data = read_forcing("Hyde_data_2010_2016.csv",
                            gpara['start_time'],
                            gpara['end_time'],
                            cols=['NEE','GPP','LE','ET'],
                            na_values="NaN")

    # plot some results as well
    fmonth = 4
    lmonth = 5

    flux=pd.DataFrame(cm_flx, index=Forc.index)

    ix = pd.rolling_mean(Forc['Prec'], 48, 1)
    dryc = np.ones(len(Forc))
    f = np.where(ix > 0)[0]  # wet canopy indices
    dryc[f] = 0.0

    x = Data.GPP
    y = flux['GPP']
    p = np.polyfit(x, y, 1)

    plt.figure()
    plt.subplot(331)
    plt.plot(x, y, 'mo')
    plt.plot([0, 25], [0, 25], 'k-')
    plt.xlabel('GPP meas')
    plt.ylabel('GPP mod')
    plt.title('y = %.2f + %.2f' % (p[0], p[1]))

    ff = np.where(~np.isnan(Data.LE))[0]
    x = Data.LE
    y = flux.LE
    p = np.polyfit(x[ff], y[ff], 1)

    plt.subplot(334)
    plt.plot(x[ff], y[ff], 'co')
    plt.plot([0, 600], [0, 600], 'k-')
    plt.xlabel('LE meas')
    plt.ylabel('LE mod')
    plt.title('y = %.2f + %.2f' % (p[0], p[1]))
    
    # water fluxes
    del x, y
    f = np.where((Forc.index.month >= fmonth) & (Forc.index.month <= lmonth) & (dryc == 1))[0]  # & (Forc.Rnet.values > 20.0)
    
    flux['ETmod'] = (flux.Transp + flux.Efloor) # + flux.Evap)
    Data['ETmeas'] = Data.ET * 0.0324  # mm/30min
    x = Data.ETmeas
    y = flux.ETmod
    p = np.polyfit(x[f], y[f], 1)
    
    plt.subplot(337)
    plt.plot(x[f], y[f], 'ro')
    plt.plot([0, 0.25], [0, 0.25], 'k-')
    plt.xlabel('ET meas')
    plt.ylabel('ET mod')
    plt.title('y = %.2f + %.2f' % (p[0], p[1]))
    
    plt.subplot(332)
    plt.plot(Data.GPP, 'k.-')
    plt.plot(flux.GPP, 'm-')
    plt.ylabel('GPP')
    plt.title("Timeseries")
    
    plt.subplot(335)
    plt.plot(Data.LE, 'k.-')
    plt.plot(flux.LE, 'c-')
    plt.ylabel('LE')
    
    plt.subplot(338)
    plt.plot(Data.ETmeas, 'k.-')
    plt.plot(flux.ETmod, 'r-')
    plt.ylabel('ET')
    
    #plt.figure()
    #plt.plot(ETmod[f], 'r.-', dt*flux.Efloor[f], 'g-')
    
    # ensemble diurnal cycles
    meaa = diurnal_cycle(Data[['GPP', 'LE', 'ETmeas']].iloc[f], ap='hour')
    moda = diurnal_cycle(flux[['GPP', 'LE', 'ETmod']].iloc[f], ap='hour')
    
    plt.subplot(333); plt.ylabel('GPP'); plt.xlabel('time')
    plt.plot(meaa['GPP']['mean'], 'ko-', moda['GPP']['mean'], 'mo-')
    plt.xlim([0, 24])
    plt.title("Diurnal cycles")
    
    plt.subplot(336); plt.ylabel('LE'); plt.xlabel('time')
    plt.plot(meaa['LE']['mean'], 'ko-', moda['LE']['mean'], 'co-')
    plt.xlim([0, 24])
    
    plt.subplot(339); plt.ylabel('ET'); plt.xlabel('time')
    plt.plot(meaa['ETmeas']['mean'], 'ko-', moda['ETmod']['mean'], 'ro-')
    plt.xlim([0, 24])