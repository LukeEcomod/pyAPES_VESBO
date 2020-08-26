# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 15:42:17 2020

@author: 03081268
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools

from pyAPES_utilities.bigleaf_functions import saturation_vapor_pressure, water_use_efficiency,\
    light_use_efficiency, canopy_conductance, eq_evap

from pyAPES_utilities.timeseries import running_mean

eps = np.finfo(float).eps  # machine epsilon

cfact = 0.021618 # umol CO2 m-2 s-1 to gC m-2 30min-1
etfact = 0.0324 # mmmol H2O m-2 s-1 to mm 30min-1
#: [J mol\ :sup:`-1`\ ], latent heat of vaporization at 20\ :math:`^{\circ}`\ C
LATENT_HEAT = 44100.0
#: [kg mol\ :sup:`-1`\ ], molar mass of H\ :sub:`2`\ O
MOLAR_MASS_H2O = 18.015e-3
#: [umol m\ :sup:`2` s\ :sup:`-1`\ ], conversion from watts to micromol
PAR_TO_UMOL = 4.56

def scen_differences(results, scens, rsim=0, nsim=None):
    """
    Differences between scenarios
    Args:
        res - xarray dataset
        nsim - list of int
        rsim - reference simulation index
    Returns:
        none
    """
    
    dt = 1800.0
    wet_dur = 24
    v = ['LAI', 'CO2', 'Vmax', 'GPP', 'Reco', 'ET', 'tr', 'WUE', 'LUE', 'Gs', 'alpha', 'CiCa']
    
    if nsim is None:
        nsim = len(results.simulation)
    
    cres = {k: [] for k in v }
    header = []
    for k in range(nsim):
        print(k, nsim)
        par = results['forcing_par'][:,k].values * PAR_TO_UMOL
        rg = results['forcing_par'][:,k].values + results['forcing_nir'][:,k].values
        prec = results['forcing_precipitation'][:,k].values
        pamb = 1e-3 * results['forcing_pressure'][:,k].values #kPa
        ta = results['forcing_air_temperature'][:,k].values
        h2o = results['forcing_h2o'][:,k].values
        co2 = results['forcing_co2'][:,k].values
        rnet = results['canopy_Rnet'][:,k].values
        gpp = results['canopy_GPP'][:,k].values
        reco = results['canopy_Reco'][:,k].values
        #nee = results['canopy_NEE'][:,k].values
        et = 1e3 * results['canopy_LE'][:,k].values / LATENT_HEAT #mmol s-1
        es,_ = saturation_vapor_pressure(ta) 
        vpd = 1e-3*es - h2o *pamb # kPa
        
        # dry canopy times:
        N = len(prec)
        ix = running_mean(prec, wet_dur, center=False)
        dryc = np.ones(N)
        f = np.where(ix > 0)[0]  # wet canopy indices
        dryc[f] = 0
        del ix
        ix = np.where(dryc == 1)
        
        
        header.append(scens[k][0] + ' ' + scens[k][1] + ' ' + scens[k][2])
        
        cres['LAI'].append(scens[k][0].split('_')[1])
        cres['CO2'].append(scens[k][1].split('_')[1])
        cres['Vmax'].append(scens[k][2].split('_')[1])
        
        # compute big-leaf parameters
        cres['GPP'].append(sum(gpp*cfact))
        cres['Reco'].append(sum(reco*cfact))
        cres['ET'].append(sum(et[et>0]*etfact))
        cres['tr'].append(sum(results['canopy_transpiration'][:,k].values * 1e3 / MOLAR_MASS_H2O *dt)) #molm-2
        
        cres['LUE'].append(np.nanmean(light_use_efficiency(gpp[ix], par[ix]))) # mol CO2 / mol Par

        cres['WUE'].append(np.nanmean(water_use_efficiency(gpp[ix], et[ix], par[ix])))  # mmolCO2 / mol H2O

        cres['Gs'].append(np.nanmean(canopy_conductance(et[ix], vpd[ix], pamb[ix], rg[ix])))
        ETeq = 1e3 * eq_evap(rnet, ta, P=1e3*pamb, units='mol')
        cres['alpha'].append(sum(et[et>0]) / sum(ETeq[ETeq>0]))
        
        ix = np.where((dryc == 1) & (par > 100))
        gsc = 1 / 1.6 * canopy_conductance(et[ix], vpd[ix], pamb[ix], rg[ix])
        
        cica = 1 - gpp[ix] / (gsc * co2[ix])
        
        
        cres['CiCa'].append(np.nanmean(cica))
        
    #print(len(cres['GPP']), len(cres['LAI']), cres['LAI'])
    resu = pd.DataFrame.from_dict(cres)
    resu.index = header
        
    # make figure
    msize = [8, 5, 10]
    symb = ['o', 's', '*']
    colo = ['k', 'r', 'g']

    spec = itertools.product(msize, colo, symb)
    spec = list(spec)

    fig, ax = plt.subplots(3,3, sharex=True)
    fig.set_size_inches(11.6, 9.0)
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
    
    k = 0
    m = 0

    for v in ['GPP', 'Reco', 'ET', 'tr', 'WUE', 'LUE', 'Gs', 'alpha', 'CiCa']:
        for n in range(nsim):
            ax[k,m].plot(n, resu[v].iloc[n] / resu[v].iloc[rsim], markersize=spec[n][0], marker=spec[n][2], color=spec[n][1])
        
        ax[k,m].set_ylabel(v)
        m +=1
        if m==3:
            k +=1
            m = 0
    
    for m in range(3):
        ax[2,m].set_xticks(range(nsim))
        ax[2,m].set_xticklabels(resu.index, rotation=90.0, fontsize=8)
        ax[2,m].tick_params(direction='in')
    
    return cres, resu


def summarize_scenarios(results, para, nsim=None):
    
    dt = 1800.0
    wet_dur = 24
    v = ['LAI', 'Ca', 'Vmax', 'g1', 'NEE', 'GPP', 'Reco', 'ET', 'tr', 'WUE',
         'LUE', 'WUEs', 'LUEs', 'Gs', 'alpha','EF', 'CiCa']
    
    if nsim is None:
        nsim = len(results.simulation)
    

    cres = {i: [] for i in v }
    
    for k in range(nsim):
        print(k, nsim)
        X = results.isel(simulation=k)
        par = X['forcing_par'].values * PAR_TO_UMOL
        rg = X['forcing_par'].values + X['forcing_nir'].values
        prec = X['forcing_precipitation'].values
        pamb = 1e-3 * X['forcing_pressure'].values  #kPa
        ta = X['forcing_air_temperature'].values
        h2o = X['forcing_h2o'].values
        co2 = X['forcing_co2'].values
        rnet = X['canopy_Rnet'].values
        gpp = X['canopy_GPP'].values
        reco = X['canopy_Reco'].values
        nee = X['canopy_NEE'].values
        et = 1e3 * X['canopy_LE'].values / LATENT_HEAT  #mmol s-1
        es,_ = saturation_vapor_pressure(ta) 
        vpd = 1e-3*es - h2o *pamb #  kPa
        
        # dry canopy times:
        N = len(prec)
        ix = running_mean(prec, wet_dur, center=False)
        dryc = np.ones(N)
        f = np.where(ix > 0)[0]  # wet canopy indices
        dryc[f] = 0
        del ix
        ix = np.where((dryc == 1) & (par >100))
        
        # compute big-leaf parameters
        cres['NEE'].append(sum(nee*cfact))
        cres['GPP'].append(sum(gpp*cfact))
        cres['Reco'].append(sum(reco*cfact))
        cres['ET'].append(sum(et[et>0]*etfact))
        cres['tr'].append(sum(X['canopy_transpiration'].values * 1e3 / MOLAR_MASS_H2O *dt))  #molm-2
        
        cres['LUE'].append(np.nanmean(light_use_efficiency(gpp[ix], par[ix]))) # mol CO2 / mol Par
        cres['WUE'].append(np.nanmean(water_use_efficiency(gpp[ix], et[ix], par[ix])))  # mmolCO2 / mol H2O
        cres['Gs'].append(np.nanmean(canopy_conductance(et[ix], vpd[ix], pamb[ix], rg[ix])))
    
        cres['LUEs'].append(sum(gpp[ix]) / sum(par[ix]))
        cres['WUEs'].append(sum(gpp[ix]) / sum(et[ix]))
        
        ETeq = 1e3 * eq_evap(rnet, ta, P=1e3*pamb, units='mol')
        cres['alpha'].append(sum(et[et>0]) / sum(ETeq[et>0]))
        cres['EF'].append(sum(X['canopy_LE'][ix] / sum(rnet[ix])))
        
        gsc = 1 / 1.6 * canopy_conductance(et[ix], vpd[ix], pamb[ix], rg[ix])
        
        cica = 1 - gpp[ix] / (gsc * co2[ix])
        
        cres['CiCa'].append(np.nanmean(cica))
        
        # scenario params
        cres['LAI'].append(para['LAI'][k])
        cres['Ca'].append(para['Ca'][k])
        cres['Vmax'].append(para['Vmax'][k])
        cres['g1'].append(para['g1'][k])
        
    resu = pd.DataFrame.from_dict(cres)
    resu.index = np.arange(0, nsim)

    return resu