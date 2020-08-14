# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 13:11:23 2019

@author: slauniai
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
sys.path.append(r'c:\MLM')

from pathlib import Path

# statistical tools
import scipy as sp

def linear_regression(x, y, alpha=0.05):
    """
    Univariate linear least-squares regression with confidence intervals using
    scipy.stats.linregress
    adapted from: http://bagrow.info/dsv/LEC10_notes_2014-02-13.html
    Args: x, y, alpha
    Returns: res (dict)
    """
    
    n = len(x)    
    
    slope, b0, r_value, p_value, serr = sp.stats.linregress(x, y)
    
    # mean squared deviation
    s2 = 1./n * sum([(y[i] - b0 - slope * x[i])**2 for i in range(n)])
    
    #confidence intervals of slope and intercept
    xx = x * x
    c = -1 * sp.stats.t.ppf(alpha/2.,n-2)
    bb1 = c * (s2 / ((n-2) * (xx.mean() - (x.mean())**2)))**.5
    ci_slope = [slope - bb1, slope +bb1]
    
    bb0 = c * ((s2 / (n-2)) * (1 + (x.mean())**2 / (xx.mean() - (x.mean())**2)))**.5
    ci_interc = [b0 - bb0, b0 + bb0]
    
    res = {'model':[slope, b0], 'r': r_value, 'p': p_value, 'ci': [ci_slope, ci_interc]}
    return res

def mann_kendall_test(x, alpha=0.05):
    """
    This function is derived from code originally posted by Sat Kumar Tomer
    (satkumartomer@gmail.com)
    See also: http://vsp.pnnl.gov/help/Vsample/Design_Trend_Mann_Kendall.htm

    The purpose of the Mann-Kendall (MK) test (Mann 1945, Kendall 1975, Gilbert
    1987) is to statistically assess if there is a monotonic upward or downward
    trend of the variable of interest over time. A monotonic upward (downward)
    trend means that the variable consistently increases (decreases) through
    time, but the trend may or may not be linear. The MK test can be used in
    place of a parametric linear regression analysis, which can be used to test
    if the slope of the estimated linear regression line is different from
    zero. The regression analysis requires that the residuals from the fitted
    regression line be normally distributed; an assumption not required by the
    MK test, that is, the MK test is a non-parametric (distribution-free) test.
    Hirsch, Slack and Smith (1982, page 107) indicate that the MK test is best
    viewed as an exploratory analysis and is most appropriately used to
    identify stations where changes are significant or of large magnitude and
    to quantify these findings.

    Input:
        x:   a vector of data
        alpha: significance level (0.05 default)

    Output:
        trend: tells the trend (increasing, decreasing or no trend)
        h: True (if trend is present) or False (if trend is absence)
        p: p value of the significance test
        z: normalized test statistics

    Examples
    --------
      >>> x = np.random.rand(100)
      >>> trend,h,p,z = mk_test(x,0.05)

    """
    n = len(x)

    # calculate S
    s = 0
    for k in range(n-1):
        for j in range(k+1, n):
            s += np.sign(x[j] - x[k])

    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g:  # there is no tie
        var_s = (n*(n-1) *(2*n+5))/18
    else:  # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = sum(x == unique_x[i])
        var_s = (n*(n-1)*(2*n+5) - np.sum(tp*(tp-1)*(2*tp+5)))/18

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else: # s == 0:
        z = 0

    # calculate the p_value
    p = 2*(1 - sp.stats.norm.cdf(abs(z)))  # two tail test
    h = abs(z) > sp.stats.norm.ppf(1-alpha/2)

    if (z < 0) and h:
        trend = 'decreasing'
    elif (z > 0) and h:
        trend = 'increasing'
    else:
        trend = 'no trend'

    return trend, h, p, z
    
def trend_breakpoint(x, y, min_obs=4, figs=True):
    """ 
    Seeks for trend breakpoints using Chow test. Trends are computed using linear
    least-squares estimate
    Args:
        x - independent variable
        y - dependent variable
        min_obs - minimum values for computing fit
    Returns:
        res (dict)
    """
    def find_rss(x, y, model):
        """ sum of squared residuals
        """
        ymod = model[1] + model[0] * x
        return sum((y - ymod)**2), len(x)
    
    def chow_test(x1, y1, x2, y2, conf=0.9):
        
        # linear least-square regressions
        f = linear_regression(np.append(x1, x2), np.append(y1, y2), alpha=0.05)
        f1 = linear_regression(x1, y1, alpha=0.05)
        f2 = linear_regression(x2, y2, alpha=0.05)
        
        # sum of squared residuals
        rss_total, n_total = find_rss(np.append(x1, x2), np.append(y1, y2), f['model'])
        rss_1, n_1 = find_rss(x1, y1, f1['model'])
        rss_2, n_2 = find_rss(x2, y2, f2['model'])
        
        # F-test
        chow_nom = (rss_total - (rss_1 + rss_2)) / 2
        chow_denom = (rss_1 + rss_2) / (n_1 + n_2 - 4)
        
        F = chow_nom / chow_denom
        
        if not F:
            F = 1.0
        
        # Chow-test p-value
        df1 = 2
        df2 = len(x1) + len(x2) - 4
    
        # The survival function (1-cdf) is more precise than using 1-cdf,
        # this helps when p-values are very close to zero.
        # -f.logsf would be another alternative to directly get -log(pval) instead.
        p_val = sp.stats.f.sf(F, df1, df2)
        
        return p_val, f1, f2          
        
    # compute
    
    x0 = x - x[0]
    N = len(x)

    results = {'p_value': [], 'mod_1': [], 'mod_2': []}
    
    k = 0 # - min_obs
    while k < N:
        x1 = x0[0:k]
        y1 = y[0:k]
        x2 = x0[k:]
        y2 = y[k:]
        
        if len(x2) >= min_obs and len(x1) >= min_obs:
            print(k, x1, x2)
            p, f1, f2 = chow_test(x1, y1, x2, y2)
            f1['x'] = x[0:k]
            f2['x'] = x[k:]
        elif len(x2) < min_obs:
            p, f2 = np.NaN, None
            f1 = linear_regression(x1, y1, alpha=0.05)
            f1['x'] = x[0:k]
        elif len(x1) < min_obs:
            p, f1 = np.NaN, None
            f2 = linear_regression(x2, y2, alpha=0.05)
            f2['x'] = x[k:]            
        #print(t[0:k], t[k:N], p)
        
        results['p_value'].append(p)
        results['mod_1'].append(f1)
        results['mod_2'].append(f2)
        k += 1
    # breakpoint is defined as point where p_value == min(p_value) and p_value < 0.05
    min_p = np.nanmin(results['p_value'])
    if min_p <= 0.05:
        #ix = int(np.where(results['p_value'] == min_p) [0])
        ix = int(max(np.where(np.array(results['p_value']) <= 0.05)[0]))

        results['changepoint'] =  float(x[ix])
    else:
        results['changepoint'] = None
    
    if figs:
        plt.figure()
        plt.plot(x, y, 'ro')
        f = linear_regression(x0, y, alpha=0.05)        
        fit = f['model'][0] * x0 + f['model'][1]
        txt = '%.2f (%.2f / %.2f) x, p=%.2f, r2=%.2f' %(f['model'][0], f['ci'][0][0], f['ci'][0][1],
                                                                             f['p'], f['r']**2)
        plt.plot(x, fit, 'k-', label=txt)

        if min_p <= 0.05:
            c1 = results['mod_1'][ix]
            c2 = results['mod_2'][ix]
            
            fit1 = c1['model'][0] * (c1['x'] - x[0]) + c1['model'][1]
            txt1 = '%.2f (%.2f / %.2f) x, p=%.2f, r2=%.2f' %(c1['model'][0], c1['ci'][0][0], c1['ci'][0][1],
                                                                     c1['p'], c1['r']**2)
            fit2 = c2['model'][0] * (c2['x'] - x[0]) + c2['model'][1]
            txt2 = '%.2f (%.2f / %.2f) x, p=%.2f, r2=%.2f' %(c2['model'][0], c2['ci'][0][0], c2['ci'][0][1],
                                                                     c2['p'], c2['r']**2)            
            plt.plot(c1['x'], fit1, 'k--', label=txt1)
            plt.plot(c2['x'], fit2, 'k:', label=txt2)
        
        plt.legend(fontsize=8)
    
    return results

def compute_trends(D, yrs=None, crit=None, fname=None, figs=False,
                   folder=r'c:\repositories\Hyde_trends\Results\figs'):
    """
    compute trends in timeseries.
    Tests:
    Theil-Sen (theilslopes)
    Mann-Kendall test (mk_test)
    linear regression (linregress)
    
    Args:
        D - pd.dataframe
        yrs - list of years to be used
        figs - True plots
    Returns:
        
    """
        
    crit0 = dict()
    crit0['sign_lev_1'] = 0.01
    crit0['sign_lev_2'] = 0.05
    crit0['sign_lev_3'] = 0.1
    
    if crit is None:
        crit = crit0.copy()

    D = D.replace([np.inf, -np.inf], np.nan)
    
    # select only years desired
    if yrs:
        tv = list(D.index.year)
        ix = [j for j, e in enumerate(tv) if e in set(yrs)]
        D = D.iloc[ix]
    
    D.replace([np.inf, -np.inf], np.nan)
    
    t = D.index.values.astype(float)
    x = t - t[0]
    
    # create dataframe for storing trend results     
    cols = ['ts', 'ts0', 'tsl', 'tsu', 'ts_p', 'ls', 'ls0', 'ls_r2', 'ls_p', 'ls_err', 'mk_tr', 'mk_p']
    T = pd.DataFrame(data=None, columns=cols, index=D.columns) 
    
    for key in D.columns:
        m = D[key].dropna(axis=0)
        t = m.index.values.astype(float)
        x = t - t[0]
        y = m.values

        # theil-sen estimator (median slope between pair of points)
        for a in [crit['sign_lev_1'], crit['sign_lev_2'], crit['sign_lev_3']]:
            ts, ts0, tsl, tsu = sp.stats.theilslopes(y, x, alpha=a)
            ts_a = a
            if np.sign(tsl) == np.sign(tsu): 
                print(key,a); break
        
        # mann-kendall test
        for a in [crit['sign_lev_2'], crit['sign_lev_1']]:
            mk_tr, _, mk_p, _, = mann_kendall_test(y, alpha=a)
            if mk_p <= a: break
    
        # linear least-squares regression
        for a in [crit['sign_lev_2'], crit['sign_lev_1']]:
            ls, ls0, r, ls_p, ls_err = sp.stats.linregress(x, y)

            if ls_p <= a: break
        
        T.loc[key][cols] = [ts, ts0, tsl, tsu,ts_a, ls, ls0, r**2, ls_p, ls_err, mk_tr, mk_p]
            
    # drop rows where no significant trends
    T = T.dropna(axis=0, how='all')
    
    # save to file
    if fname:
        T.to_csv(Path(folder) / (fname + '.csv'), float_format='%.3f', sep='\t')
        T.to_html(Path(folder) / (fname + '.html'))

    if figs:
        figlist = []

        # plot trends
        for key in list(T.index):
            # plot data and Theil-Sen slopes and linregress slopes
            m = D[key].dropna(axis=0)
            t = m.index.values.astype(float)
            x = t - t[0]
            x0 = np.mean(x)
            xx = x - x0
            m = m.values
            y = T.loc[key]
            
            figname = key
            ttxt = 'TS: [%.3f (%.3f/%.3f)] LR: [%.3f (+/-%.3f); r2=%.2f, p=%.3f]' %(
                    y.ts, y.tsl, y.tsu, y.ls, y.ls_err, y.ls_r2, y.ls_p)
            fig = plt.figure('A_' + key)
            
            f = y.ls0 + y.ls*x
            f0 = y.ls0 + y.ls*x0
    
            fu = f0 + xx*(y.ls + y.ls_err) # * N**0.5)
            fl =  f0 + xx*(y.ls - y.ls_err) # *  N**0.5)

            plt.fill_between(t, fl, fu, color='b', alpha=0.2)
            plt.plot(t, m, 'ko', t, f, 'b-')
            
            # theil-sen intercepts for boundary lines
            ts0l = np.median(m) - y.tsl*np.median(x)
            ts0u = np.median(m) - y.tsu*np.median(x)
            plt.plot(t, y.ts0 + y.ts*x, 'r-', t, ts0l + y.tsl*x, 'r--', t, ts0u + y.tsu*x, 'r--')
            plt.title(ttxt, fontsize=8); plt.ylabel(figname)
            plt.xticks(t, rotation=45.0)
            
            #fn = Path(folder) / (fname + '_' + figname + '.png')
            #fig.savefig(fn, dpi=600)

            figlist.append(fig)

        # create multipage pdf
        from matplotlib.backends.backend_pdf import PdfPages
        fn = Path(folder) / (fname + '.pdf')

        with PdfPages(fn) as pdf:
            # As many times as you like, create a figure fig and save it:
            for f in figlist:
                pdf.savefig(f)
    return T

def seek_critical_periods(data, cols, window_size=30, step=1, remove_trend=False, normalize=False):
    """
    Find critical periods creating annual anomalies. 
    Follows Le Maire et al. 2010 JGR.
    
    Args:
        data - 30min flux dataframe
        meteo  - 30min meteo dataframe
        windowsize - window days
        remove_trend - True removes linear trends from variables
        normalize - True normalizes all variables to zero mean unit variance before computations
    """
    from scipy.signal import detrend
    
#    cfact = 0.021618 # umol m-2 s-1 to gC m-2 30min-1
#    etfact = 0.0324 # mmmol m-2 s-1 to mm 30min-1
#    
#
#    dcols0 = ['NEE', 'NEE_1', 'GPP', 'GPP_1', 'GPP_2', 'GPP_4','Reco', 'Reco_1', 'Reco_2',
#              'Reco_4', 'ET']
#    mcols0 = ['U', 'ust', 'Ta', 'RH', 'CO2', 'H2O', 'O3', 'Prec', 'P', 'dirPar', 'diffPar',
#              'dirNir', 'diffNir', 'Rnet', 'LWin', 'LWout', 'LWnet', 'Tsh', 'Tsa', 'Tsc',
#              'Wh', 'Wa', 'Wc', 'emiatm', 'cloudfract', 'Rew', 'diff_fr']
    
    # resample to daily 
    data = data[cols].resample('D').mean()
    
    # remove leap years
    data = data[~((data.index.month == 2) & (data.index.day == 29))]
    n = len(np.unique(data.index.year))
    doy = np.tile(np.arange(1, 366), n)
    data['doy'] = doy
    data['year'] = data.index.year
    #print(max(data.doy))
    
    # annual data
    A = data.resample('A').mean()

    # make data array
    dcols = ['year', 'doy']
    dcols.extend(cols)
    data = data[dcols].values
    N = len(data)
    
    # compute means for 'windows': the value is for past 'window_size'
    
    dd = np.ones(np.shape(data)) * np.NaN
    k = 0
    while k < N:
        #print(k)
        dd[k,0:2] = data[k,0:2]

        if k >= window_size:
            dd[k,2:] = np.mean(data[k-window_size:k, 2:], axis=0)

        k += step
    
    #convert to dataframe
    dd = pd.DataFrame(dd, columns=dcols)
    
    ix = dd.year >= 2001
    dd = dd.loc[ix,:]
    
    ix = A.index.year >= 2001
    A = A.loc[ix,:]
    
    # seek period of maximum correlation
    res = dict.fromkeys(cols, np.ones((366, 3))*np.NaN)

    for key in res.keys():
        x = A[key].values
        #x = x[1:]
        if remove_trend:
            x = detrend(x) 
        x = x - np.mean(x)

        for k in range(1, 366): #doy
            # seek index
            f = (dd.doy == k)
            y = dd[key].loc[f].values
            #y = y[1:] # omit 1st year
            del f
            if remove_trend:
                y = detrend(y)          
            y = y - np.mean(y)

            # pearson corr. coeff and period to annual variance ratio
            r = np.corrcoef(x, y)
            v = np.var(y) / np.var(x)
            
            res[key][k,0] = r[0,1]
            res[key][k,1] = v
     
        plt.figure()
        plt.subplot(211); plt.plot(res[key][:,0], 'g-'); plt.ylim(-0.5, 1.0)
        plt.title(key); plt.ylabel('r [-]')
        plt.subplot(212); plt.plot(res[key][:,1], 'g-'); plt.ylabel(r'\sigma_v/\sigma_a')

    return dd, res

#def seek_critical_periods(data, meteo, window_size=30, step=1, detrend=False, normalize=False):
#    """
#    Find critical periods creating annual anomalies. 
#    Follows Le Maire et al. 2010 JGR.
#    
#    Args:
#        data - 30min flux dataframe
#        meteo  - 30min meteo dataframe
#        windowsize - window days
#        detrend - True removes linear trends from variables
#        normalize - True normalizes all variables to zero mean unit variance before computations
#    """
#    from scipy.signal import detrend
#    
#    cfact = 0.021618 # umol m-2 s-1 to gC m-2 30min-1
#    etfact = 0.0324 # mmmol m-2 s-1 to mm 30min-1
##    def is_leap_and_29Feb(s):
##        # returns mask for leap year days = True
##        return ((s.index.year % 4 == 0) & ((s.index.year % 100 != 0) | 
##                (s.index.year % 400 == 0)) & (s.index.month == 2) & 
##                (s.index.day == 29))
#    
#
#    dcols0 = ['NEE', 'NEE_1', 'GPP', 'GPP_1', 'GPP_2', 'GPP_4','Reco', 'Reco_1', 'Reco_2',
#              'Reco_4', 'ET']
#    mcols0 = ['U', 'ust', 'Ta', 'RH', 'CO2', 'H2O', 'O3', 'Prec', 'P', 'dirPar', 'diffPar',
#              'dirNir', 'diffNir', 'Rnet', 'LWin', 'LWout', 'LWnet', 'Tsh', 'Tsa', 'Tsc',
#              'Wh', 'Wa', 'Wc', 'emiatm', 'cloudfract', 'Rew', 'diff_fr']
#    
#    # daily 
#    data = data[dcols0].resample('D').sum() * cfact
#    meteo = meteo[mcols0].resample('D').mean()
#    
#    # remove leap years
#    data = data[~((data.index.month == 2) & (data.index.day == 29))]
#    n = len(np.unique(data.index.year))
#    doy = np.tile(np.arange(1, 366), n)
#    data['doy'] = doy
#    data['year'] = data.index.year
#    #print(max(data.doy))
#    
#    # annual data
#    A = data.resample('A').sum()
#    #A.index = A.index.year.values
#    print(A.index, A.head())
#    plt.figure()
#    plt.subplot(211); plt.plot(A.NEE)
#    plt.subplot(212); plt.plot(data.NEE)
#    
#    # make data array
#    dcols = ['year', 'doy']
#    dcols.extend(dcols0)
#    data = data[dcols].values
#    N = len(data)
#    
#    meteo = meteo[~((meteo.index.month == 2) & (meteo.index.day == 29))]  
#    meteo['doy'] = doy
#    meteo['year'] = meteo.index.year
#    
#    M = meteo.resample('A').mean()
#    # make meteo array
#    mcols = ['year', 'doy']
#    mcols.extend(mcols0)
#    #print(mcols)
#    #print(dcols)
#    meteo = meteo[mcols].values
#    
#    # compute cumulative sums and means for 'windows': the value is for past 'window_size'
#    
#    dd = np.ones(np.shape(data)) * np.NaN
#    dm = np.ones(np.shape(meteo)) * np.NaN
#    k = window_size
#    k = 0
#    while k < N:
#        #print(k)
#        dd[k,0:2] = data[k,0:2]
#        dm[k,0:2] = data[k,0:2]
#
#        if k >= window_size:
#            dd[k,2:] = np.sum(data[k-window_size:k, 2:], axis=0)
#            dm[k,2:] = np.mean(meteo[k-window_size:k, 2:], axis=0)
#
#        k += step
#    
#    #convert to dataframes
#    dd = pd.DataFrame(dd, columns=dcols)
#    dm = pd.DataFrame(dm, columns=mcols)
#    
#    ix = dd.year >= 2001
#    dd = dd.loc[ix,:]
#    dm = dm.loc[ix,:]
#    
#    ix = A.index.year >= 2001
#    A = A.loc[ix,:]
#    #print(A)
#    
#    # seek maximum period of maximum correlation
#    res = dict.fromkeys(dcols0, np.ones((366, 3))*np.NaN)
#
#    for key in res.keys():
#        x = A[key].values
#        #x = x[1:]
#        x = detrend(x) 
#        x = x - np.mean(x)
#
#        for k in range(1, 366): #doy
#            # seek index
#            f = (dd.doy == k)
#            y = dd[key].loc[f].values
#            #y = y[1:] # omit 1st year
#            del f
#
#            y = detrend(y)          
#            y = y - np.mean(y)
#            #print(key, k,  np.shape(x), np.shape(y))
#
#            # pearson corr. coeff and variance ratio
#            r = np.corrcoef(x, y)
#            v = np.var(y) / np.var(x)
#            
#            res[key][k,0] = r[0,1]
#            res[key][k,1] = v
#                #            if len(ix) > 2:
##                # deviations around mean
##                y = y[ix]
##                # y = detrend(y)
##                y = y - np.mean(y)
##                
##                x = A[key].values[ix] 
##                #x = detrend(x)
##                x = x - np.mean(x)
#                
#                
#            #res[key][k,2] = np.var(x)             
#        plt.figure()
#        plt.subplot(211); plt.plot(res[key][:,0], 'g-'); plt.ylim(-0.5, 1.0)
#        plt.title(key); plt.ylabel('r')
#        plt.subplot(212); plt.plot(res[key][:,1], 'g-'); plt.ylabel('vr')
#        #plt.subplot(313);  plt.plot(res[key][:,1]*res[key][:,2], 'r-'); plt.plot(res[key][:,2], 'r-');
#
#    # seek maximum period of maximum correlation
#    mres = dict.fromkeys(mcols0, np.ones((366, 3))*np.NaN)
#
#    for key in mres.keys():
#        for k in range(0, 366): #doy
#            # seek index
#            f = (dm.doy == k)
#            y = dm[key].loc[f].values
#            del f
#            ix = np.where(np.isnan(y) == 0)[0] # good obs
#            
#            if len(ix) > 2:
#                # deviations around mean
#                y = y[ix]
#                # y = detrend(y)
#                y = y - np.mean(y)
#                
#                x = M[key].values[ix] 
#                #x = detrend(x)
#                x = x - np.mean(x)
#                
#                # pearson corr. coeff and variance ratio
#                r = np.corrcoef(x, y)
#                v = np.var(y) / np.var(x)
#                
#                mres[key][k,0] = r[0,1]
#                mres[key][k,1] = v
#                
#            #res[key][k,2] = np.var(x)
#          
#        plt.figure()
#        plt.subplot(211); plt.plot(mres[key][:,0], 'k-'); plt.ylim(-0.5, 1.0)
#        plt.title(key); plt.ylabel('r')
#        plt.subplot(212); plt.plot(mres[key][:,1], 'k-'); plt.ylabel('vr')
#
#
#    return dd, dm
    
    
    
    
    
    