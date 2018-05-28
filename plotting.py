# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:09:06 2018

@author: L1656
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from timeseries_tools import diurnal_cycle, yearly_cumulative
from forcing.forc_utils import read_forcing
from parameters.general_parameters import gpara
import seaborn as sns

pal = sns.color_palette("hls", 5)

def plotresults(results):
    start_time=results.date.values[0]
    end_time=results.date.values[-1]
    # Read ET
    ET_hyde = read_forcing("Hyde_data_1997_2016.csv",
                           start_time,
                           end_time,
                           cols=['ET'])
    ET_hyde.loc[:,'ET'] = ET_hyde.loc[:,'ET'] * 0.0324  # mm/30min
    # Read snow
    snow = read_forcing("FMI_jokioinen.csv",
                                start_time,
                                end_time,
                                cols=['SnowD'])
    # Read weir
    weir = read_forcing("Lettosuo_weir.csv",
                           start_time,
                           end_time,
                           cols=['runf'])
    weir.loc[:,'runf'] = weir.loc[:,'runf'] * 1e-3 *gpara['dt']
    # Read gwl
    gwl_meas = read_forcing("Lettosuo_gwl.csv",
                               start_time,
                               end_time,
                               cols=['WT_E','WT_N','WT_W','WT_S'])
    
    dates = results.date.values
    variables=['canopy_moss_evaporation',
               'canopy_transpiration',
               'canopy_evaporation',
               'soil_subsurface_drainage',
               'soil_surface_runoff']
    yearly_cum = yearly_cumulative(results, variables)
    yearly_cum_ET = yearly_cumulative(ET_hyde, ['ET'])
    
    plt.figure()
    plt.subplot(8,3,(1,2))
    plotxarray(results, ['forcing_air_temperature'], colors=pal, xticks=False)
    plt.subplot(8,3,(4,5))
    plotxarray(results, ['forcing_precipitation', 'canopy_throughfall'], colors=pal[3:], xticks=False, m_to='mmperh')
    plt.subplot(8,3,(7,8))
    plt.stackplot(dates, 1000 * yearly_cum, labels=variables, colors=pal)
    plt.xlim([results.date.values[0], results.date.values[-1]])
#    plt.ylim(0,700)
    plt.plot(ET_hyde.index, yearly_cum_ET[0], 'k', linewidth=1, label='ET_hyde')
    plt.ylabel('[mm]')
    plt.legend(bbox_to_anchor=(1.01,0.5), loc="center left", fontsize=8)
    plt.subplot(8,3,(10,11))
    plotxarray(results, ['canopy_transpiration', 'canopy_evaporation', 'canopy_moss_evaporation'],
               colors=[pal[1]] + [pal[2]] + [pal[0]], xticks=False, m_to='mmperh')
    plt.subplot(8,3,(13,14))
    plotxarray(results, ['soil_total_runoff'],
               colors='k', xticks=False, m_to='mmperh')
    plt.plot(weir.index, weir['runf'],':k', linewidth=1, label='meas_runoff')
    plt.legend(bbox_to_anchor=(1.01,0.5), loc="center left", fontsize=8)
    plt.subplot(8,3,(16,17))
    plotxarray(results, ['canopy_phenostate', 'canopy_LAI'], colors=pal, xticks=False)
    plt.subplot(8,3,(19,20))
    plt.plot(snow.index, 1e-3 * snow['SnowD'], 'gray', linewidth=1, label='meas_snow_depth [cm]')
    plotxarray(results, ['canopy_snow_water_equivalent'], colors='k', xticks=False, m_to='mm')
    plt.subplot(8,3,(22,23))
    plt.plot(gwl_meas.index, gwl_meas['WT_E'],color=pal[2], linewidth=1)
    plt.plot(gwl_meas.index, gwl_meas['WT_N'],color=pal[2], linewidth=1)
    plt.plot(gwl_meas.index, gwl_meas['WT_W'],color=pal[2], linewidth=1)
    plt.plot(gwl_meas.index, gwl_meas['WT_S'],color=pal[2], linewidth=1)
    plotxarray(results, ['soil_ground_water_level'], colors=pal[3:], xticks=True)

def plotxarray(results, variables, colors, xticks=True, m_to=False):
    ymin, ymax = 0.0, 0.0
    for k in range (0, len(variables)):
        results[variables[k]].plot(color=colors[k], label=results[variables[k]].name,linewidth=1)
        ymin = min(results[variables[k]].values.min(), ymin)
        ymax = max(results[variables[k]].values.max(), ymax)
    plt.ylim(ymin, ymax)
    plt.xlim([results.date.values[0], results.date.values[-1]])
    plt.title('')
    plt.ylabel('[' + results[variables[-1]].units.split('[')[-1])
    if xticks is False:
        frame1 = plt.gca()
        frame1.axes.xaxis.set_ticklabels([])
    if m_to=='mm' or m_to=='mmperh':
        if m_to=='mm':
            conversion = 1e3
            plt.ylabel('[mm]')
        if m_to=='mmperh':
            dt = (results.date.values[1] - results.date.values[0]) / np.timedelta64(1, 's')
            conversion = 1e3 / dt * 3600
            plt.ylabel('[mm h-1]')
        frame1 = plt.gca()
        _ = frame1.axes.yaxis.set_ticklabels(
                frame1.axes.yaxis.get_ticklocs()[:] * conversion)
    plt.xlabel('')
    plt.legend(bbox_to_anchor=(1.01,0.5), loc="center left", fontsize=8)

def plotresultsMLM(results):
    Data = read_forcing("Hyde_data_1997_2016.csv",
                            gpara['start_time'],
                            gpara['end_time'],
                            cols=['NEE','GPP','LE','ET'],
                            na_values="NaN")

    dates = results.date.values

    # plot some results as well
    fmonth = 4
    lmonth = 5

    ix = pd.rolling_mean(results.forcing_precipitation.values, 48, 1)
    dryc = np.ones(len(dates))
    f = np.where(ix > 0)[0]  # wet canopy indices
    dryc[f] = 0.0

    x = Data.GPP
    y = results.canopy_GPP.values
    p = np.polyfit(x, y, 1)

    plt.figure()
    plt.subplot(331)
    plt.plot(x, y, 'mo')
    plt.plot([0, 25], [0, 25], 'k-')
    plt.xlabel('GPP meas')
    plt.ylabel('GPP mod')
    plt.title('y = %.2f + %.2f' % (p[0], p[1]))
    plt.axis('equal')

    ff = np.where(~np.isnan(Data.LE))[0]
    x = Data.LE
    y = results.canopy_LE.values
    p = np.polyfit(x[ff], y[ff], 1)

    plt.subplot(334)
    plt.plot(x[ff], y[ff], 'co')
    plt.plot([0, 600], [0, 600], 'k-')
    plt.xlabel('LE meas')
    plt.ylabel('LE mod')
    plt.title('y = %.2f + %.2f' % (p[0], p[1]))
    plt.axis('equal')

    # water fluxes
    del x, y
    months = results.date.dt.month.values
    f = np.where((months >= fmonth) & (months <= lmonth) & (dryc == 1))[0]  # & (Forc.Rnet.values > 20.0)

    ETmod = (results.canopy_transpiration.values + results.canopy_moss_evaporation.values)*1e3
    Data['ETmeas'] = Data.ET * 0.0324  # mm/30min
    x = Data.ETmeas
    y = ETmod
    p = np.polyfit(x[f], y[f], 1)
    
    plt.subplot(337)
    plt.plot(x[f], y[f], 'ro')
    plt.plot([0, 0.25], [0, 0.25], 'k-')
    plt.xlabel('ET meas')
    plt.ylabel('ET mod')
    plt.title('y = %.2f + %.2f' % (p[0], p[1]))
    plt.axis('equal')
    
    plt.subplot(3,3,(2,3))
    plt.plot(dates, Data.GPP, 'k.-')
    plt.plot(dates,results.canopy_GPP.values, 'm-')
    plt.ylabel('GPP')
    plt.title("Timeseries")
    
    plt.subplot(3,3,(5,6))
    plt.plot(dates,Data.LE, 'k.-')
    plt.plot(dates,results.canopy_LE.values, 'c-')
    plt.ylabel('LE')
    
    plt.subplot(3,3,(8,9))
    plt.plot(dates, Data.ETmeas, 'k.-')
    plt.plot(dates, ETmod, 'r-')
    plt.ylabel('ET')
    
    #plt.figure()
    #plt.plot(ETmod[f], 'r.-', dt*flux.Efloor[f], 'g-')
    # ensemble diurnal cycles
#    meaa = diurnal_cycle(Data[['GPP', 'LE', 'ETmeas']].iloc[f], ap='hour')
#    moda = diurnal_cycle(results.iloc[f], ap='hour')
#    
#    plt.subplot(333); plt.ylabel('GPP'); plt.xlabel('time')
#    plt.plot(meaa['GPP']['mean'], 'ko-', moda['GPP']['mean'], 'mo-')
#    plt.xlim([0, 24])
#    plt.title("Diurnal cycles")
#    
#    plt.subplot(336); plt.ylabel('LE'); plt.xlabel('time')
#    plt.plot(meaa['LE']['mean'], 'ko-', moda['LE']['mean'], 'co-')
#    plt.xlim([0, 24])
#    
#    plt.subplot(339); plt.ylabel('ET'); plt.xlabel('time')
#    plt.plot(meaa['ETmeas']['mean'], 'ko-', moda['ETmod']['mean'], 'ro-')
#    plt.xlim([0, 24])

def plot_columns(data, col_index=None):
    col_names=[]
    if col_index == None:
        col_index = range(len(data.columns))
    for i in col_index:
        col_names.append(data.columns[i])
    data[col_names].plot(kind='line',marker='o',markersize=1)
    plt.legend()
    axes = pd.plotting.scatter_matrix(data[col_names], figsize=(10, 10), alpha=.2)
    corr = data[col_names].corr().as_matrix()
    lim = [data[col_names].min().min(), data[col_names].max().max()]
    for i in range(len(col_index)):
        for j in range(len(col_index)):
            if i != j:
                idx = np.isfinite(data[col_names[i]]) & np.isfinite(data[col_names[j]])
                if idx.sum() != 0.0:
                    p = np.polyfit(data[col_names[j]][idx], data[col_names[i]][idx], 1)
                    axes[i, j].annotate("y = %.2fx + %.2f \n R2 = %.2f" % (p[0], p[1], corr[i,j]**2), (0.3, 0.9), xycoords='axes fraction', ha='center', va='center')
                    axes[i, j].plot(lim, [p[0]*lim[0] + p[1], p[0]*lim[1] + p[1]], 'r', linewidth=1)
                axes[i, j].plot(lim, lim, 'k--', linewidth=1)
            axes[i, j].set_ylim(lim)
            axes[i, j].set_xlim(lim)