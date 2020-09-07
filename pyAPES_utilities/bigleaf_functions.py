# -*- coding: utf-8 -*-
"""
Micromet- and flux-data analysis functions

Created on Wed Apr  8 22:22:08 2020

@author: slauniai
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import linregress
 
eps = np.finfo(float).eps  # machine epsilon

cfact = 0.021618 # umol CO2 m-2 s-1 to gC m-2 30min-1
etfact = 0.0324 # mmmol H2O m-2 s-1 to mm 30min-1

def vpd_from_rh(T, RH):
    """
    Computes vapor pressure deficit from temperature and relative humidity
    INPUT:
        T - temperature (degC), array or scalar
        RH - relative humidity (%), array or scalar
    OUTPUT:
        VPD - vapor pressure deficit (Pa), array
    USES:
        saturation_vapor_pressure
    """
    RH = np.array(RH)
    T = np.array(T)
    VPD, _ = (1.0-RH / 100.0)*saturation_vapor_pressure(T)  # Pa
    return VPD

def saturation_vapor_pressure(T):
    """
    Computes saturation vapor pressure with respect to free and flat water surface for given temperature T
    INPUT:
        T - temperature (degC), scalar or array
    OUTPUT: 
        esat - saturation vapor pressure (Pa)
        delta - slope of saturation vapor pressure curve (Pa degC-1)
    SOURCE:
        Campbell & Norman, 1998. Introduction to Environmental Biophysics.
    """
    # constants
    a = 611.0  # Pa
    b = 17.502  # (-)
    c = 240.97  # degC

    esat = a*np.exp(b*T / (T+c))  # Pa
    delta = b*c*esat / ((c + T)**2)  # Pa degC-1
    return esat, delta

def light_use_efficiency(gpp, par, par_lim=100.0):
    """
    Args:
        gpp [umol m-2 s-1]
        par [umol m-2 s-1]
        par_lim is threshold
    Returns:
        lue = gpp / par [umol/mmol]
    """
    x = 1e3*np.maximum(0.0, gpp / (par + eps))
    x[par < par_lim] = np.NaN
    return x

def water_use_efficiency(gpp, et, par, par_lim=100.0, et_lim=0.01):
    """
    Args:
        gpp [umol m-2 s-1]
        et [mmol m-2 s-1]
        par [umol m-2 s-1]
        par_lim is threshold
    Returns:
        wue = gpp/et [umol / mmol]
    """
    et[et < et_lim] = np.NaN
    x = np.maximum(0.0, gpp / (et + eps))
    x[par < par_lim] = np.NaN
    return x

def intrinsic_water_use_efficiency(gpp, et, vpd, par, vpd_lim=0.1, par_lim=100.0, gpp_lim=0.1, et_lim=0.1):
    """
    Args:
        gpp [umolm-2s-1]
        et [mmolm-2s-1]
        vpd [kPa]
        par [umolm-2s-1]
        par_lim, vpd_lim, gpp_lim: thresholds for acceptable periods
    Returns:
        iwue = gpp / gv = gpp / et * vpd [gC-1 kgH2O-1 hPa-1]
    """
    # iwue = gpp / gv = gpp / et * vpd
    x = cfact / etfact * gpp / et * 10 * vpd # gC-1 kgH2O-1 hPa-1
    x = np.maximum(0.0, x)
    x[(par < par_lim) | (vpd < vpd_lim) | (gpp < gpp_lim)] = np.NaN
    return x

def canopy_conductance(et, vpd, p, rg, vpd_lim=0.1, et_lim=0.1, rg_lim=100.0):
    """
    canopy conductance assuming well-coupled conditions
    Args:
        et [mmolm-2s-1]
        vpd and p [kPa]
        rg [Wm-2]
        vpd_lim [kPa]: threshold-vpd
        et_lim [mmolm-2s-1]: threshold-et
        rg_lim [Wm-2]: threshold-global rad
    Returns:
        conductance [mol m-2 s-1]
    """
    x = np.maximum(0.0, 1e-3*et / (vpd / p + eps))
    x[(rg < rg_lim) | (vpd < vpd_lim) | (et < et_lim)] = np.NaN
    return x

def priestley_taylor(et, AE, T, P=101.3):
    """ 
    priestley-taylor alpha: ratio of actual to equilibrium et
    Args:
        et [mmolm-2s-1]
        AE [Wm-2]
        T [degC]
        P [kPa]
    Returns:
        et/et_eq
    """    
    et0 = eq_evap(AE, T, P=1e3*P, units='mol')
    
    x = 1e-3*et / et0
    return x

def eq_evap(AE, T, P=101300.0, units='W'):
    """
    Calculates the equilibrium evaporation according to McNaughton & Spriggs, 1986. \n
    INPUT: 
        AE - Available energy (Wm-2)  
        T - air temperature (C)
        P - pressure (Pa)
        units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    OUTPUT: 
        equilibrium evaporation rate (Wm-2)
        constants
        NT=273.15; %0 degC in K
    """
    CP_AIR_MASS = 1004.67  # J kg-1 K-1 heat capasity of the air at constant pressure
    NT = 273.15
    Mw = 18e-3  # kg mol-1
    L = 1e3*(3147.5 - 2.37*(T + NT))  # latent heat of vaporization of water [J/kg]
    _, s = saturation_vapor_pressure(T)  # des / dT, Pa
  
    g = P*CP_AIR_MASS / (0.622*L)  # psychrom. const (Pa)

    x = np.divide((AE*s), (s+g))  # Wm-2 = Js-1m-2
    if units == 'mm':
        x = x / L  # kg m-2 s-1 = mm s-1
    elif units == 'mol':
        x = x / L / Mw  # mol m-2 s-1

    return x    

def water_stress_index(rew):
    """returns water stress index """
    x_min = 0.02
    rew_crit = 0.25 # lagergren & lindroth
    x = np.maximum(x_min, np.minimum(1.0, rew / rew_crit))
    return x

def cup_length(flx, v='GPP', crit=0.15):
    """
    estimates carbon uptake period length from flux data
    Args:
        doy, flx (gpp or nee)
        v = variable (str)
        crit - threshold
    returns: list [cup_start, cup_end, cup_length]
    """
    flx = flx.resample('1D').sum()*cfact # daily gpp
    doy = np.array(flx.index.dayofyear)

    if v == 'GPP':
        ym = np.percentile(flx['GPP'], 98.0)
        #print(ym)
        ix0 = np.where((flx['GPP'] >= crit*ym) & (doy > 50.0))[0][0]
        ix1 = np.where(flx['GPP'] >= crit*ym)[0][-1]

    elif v == 'NEE':
        #crit = 0.5 # gC m-2 d-1
        flx = -1.0*flx
        ym = np.percentile(flx['NEE'], 98.0)
        #print(ym)
        ix0 = np.where(flx['NEE'] >= crit*ym)[0][0]
        ix1 = np.where(flx['NEE'] >= crit*ym)[0][-1]

    return [doy[ix0], doy[ix1], doy[ix1] - doy[ix0]]

def thermal_growingseason(ta, tlim=5.0):
    """
    length of thermal growing season following Linderholm et al. 2008 Climatic Change
    returns gs_start, gs_end and gs_length in days
    """

    ta = ta.resample('1D').mean()
    doy = list(ta.index.dayofyear)
    ta = np.array(ta.values)

    gs_start = np.NaN
    gs_end = np.NaN
    # begining of gs is last day of first 6-day period with ta>5degC
    for d in range(6, 185):
        if all(ta[d-6:d] >= tlim):
            gs_start = doy[d]
            break
    # end of growing season is the last day of first 6-day period with ta <5degC
    for d in range(186, len(ta)):
        if all(ta[d-6:d] <= tlim):
            gs_end = doy[d]
            break
    return [gs_start, gs_end, gs_end - gs_start]

def cup_season(dat, crit=5):
    """
    Extracts CUP period start and end dates from gap-filled NEE following
    Zhu et al. 2013 Plos One:
        1) 15-day running mean filter to daily NEE
        2) 10 day period where NEE >0 for 5 first days and NEE <0 for last 5 days (CUP start)
        3) fit linear model, obtain CUP start as zero-crossing date
    Args:
        dat - pd.Dataframe with datetimeindex and column 'NEE'
        crit - number of consequtive days flux sign must be consistent to determine CUP start/end
    Returns:
        dataframe; 'SCU' = doy of CUP start, 'ECU' = doy of CUP end
    """
    dat = dat['NEE']
    D = dat.resample('1D').sum() * cfact # daily NEE gC m-2 d-1
    yrs = np.unique(D.index.year.values)
    
    res = pd.DataFrame(index=yrs, columns=['SCU','ECU'], data=None)
    for yr in yrs:
        tmp = D.loc[D.index.year == yr]
        
        # compute 15-day running mean
        y = tmp.rolling(15, win_type=None).mean()
        x = np.array(y.index.dayofyear) # doy
        y = np.array(y.values)
        
        # find SCU
        N = len(x)
        n = crit - 1
        ix = 0
        while ix < N - n: 
            a = np.all(y[ix:ix+n] > 0)
            b = np.all(y[ix+n+1:ix+2*n-1] < 0)
            if a and b:
                k = np.arange(ix, ix+2*n-1)
                s, i0, _, _, _ = linregress(x[k], y[k])
                scu = -i0 / s

                res.loc[yr, 'SCU'] = scu
                break

            ix += 1

        # compute ECU
        while ix < N-n: 
            a = np.all(y[ix:ix+n] < 0)
            b = np.all(y[ix+n+1:ix+2*n-1] > 0)
            if a and b and ix > 180:
                k = np.arange(ix, ix+2*n-1)
                s, i0, _, _, _ = linregress(x[k], y[k])
                ecu = - i0 / s

                res.loc[yr, 'ECU'] = ecu  
                break

            ix += 1      
    return res

    

def energy_balance_closure(data, fmonth=1, lmonth=12, figs=True, ust_lim=0.3):
    """
    Energy balance closure evaluated as linear regression between Rn-G and H + LE
    on annual or growing-season basis
    """
    from scipy.stats import linregress
    import matplotlib.pyplot as plt
    
    yrs = np.unique(data.year.values)
    ebc = pd.DataFrame(index=yrs, columns=['ebc', 'n'])
    for yr in yrs:
        ix = np.where((data['year'] == yr) & (data['month'] >= fmonth) & (data['month'] <= lmonth))[0]
        d = data.iloc[ix] [['Rnet', 'H', 'LE', 'Gflux', 'Qc_H', 'Qc_ET', 'ust']]
        del ix
        
        ix = np.where((d['Qc_H'] == 0) & (d['Qc_ET'] == 0) & (d['ust'] >= ust_lim))[0]
        d = d.iloc[ix]
        
        x = d['Rnet']# - d['Gflux']
        y = d['H'] + d['LE']
        
        ix = np.where((np.isnan(x) == False) & (np.isnan(y) == False))[0]
        x = x[ix]
        y = y[ix]
        n = len(y)
        if n > 2:
            slope, b0, r_value, p_value, serr = linregress(x[ix], y[ix])  
            if figs:
                plt.figure(yr)
                plt.plot(x, y, 'o'); plt.axis('square')
                xx = np.array([min(x), max(x)])
                fit = slope * xx + b0
                plt.plot(xx, fit, 'k-')
                plt.title('%d: %.2f + %.2f, r2=%.2f' %(yr, slope, b0, r_value**2))
            
            ebc.loc[yr]['ebc'] = slope
            ebc.loc[yr]['n'] = n
        
    return ebc
        
        
    
    
    