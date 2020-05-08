# -*- coding: utf-8 -*-
"""
Created on Fri Jun 08 15:00:35 2018

Note:
    migrated to python3
    - absolute imports
    - print()

@author: L1656
"""

import numpy as np
import pandas as pd
from os import listdir

from pyAPES_utilities.timeseries_tools import fill_gaps
import datetime

#: machine epsilon
EPS = np.finfo(float).eps

def create_forcingfile(meteo_fp, output_file, dir_save, lat, lon, P_unit, timezone=+2.0):
    """
    Create forcing file from meteo.
    Args:
        meteo_file (str): name of file with meteo (.csv not included)
        output_file (str): name of output file (.csv not included)
        lat (float): latitude
        lon (float): longitude
        P_unit (float): unit conversion needed to get to [Pa]
    """

    from canopy.radiation import solar_angles, compute_clouds_rad
    from canopy.micromet import e_sat

    fpar = 0.45

    dat = pd.read_csv(meteo_fp, sep=',', header='infer', encoding = 'ISO-8859-1')

    # set to dataframe index
    dat.index = pd.to_datetime({'year': dat['yyyy'],
                                'month': dat['mo'],
                                'day': dat['dd'],
                                'hour': dat['hh'],
                                'minute': dat['mm']})

    readme = ''
    cols = []

    # day of year
    dat['doy'] = dat.index.dayofyear
    cols.append('doy')
    readme += "\ndoy: Day of year [days]"

    # precipitaion unit from [mm/dt] to [m/s]
    dt = (dat.index[1] - dat.index[0]).total_seconds()
    dat['Prec'] = dat['Prec'] * 1e-3 / dt
    cols.append('Prec')
    readme += "\nPrec: Precipitation [m/s]"

    # atm. pressure unit from [XPa] to [Pa]
    dat['P'] = dat['P'] * P_unit
    cols.append('P')
    readme += "\nP: Ambient pressure [Pa]"

    # air temperature: instant and daily [degC]
    cols.append('Tair')
    readme += "\nTair: Air temperature [degC]"

#    dat['Tdaily'] = dat['Tair'].rolling(int((24*3600)/dt), 1).mean()
    dat['Tdaily'] = dat['Tair'].resample('D').mean()
    dat['Tdaily'] = dat['Tdaily'].fillna(method='ffill')

    cols.append('Tdaily')
    readme += "\nTdaily: Daily air temperature [degC]"

    # wind speend and friction velocity
    cols.append('U')
    readme += "\nU: Wind speed [m/s]"
    cols.append('Ustar')
    readme += "\nUstar: Friction velocity [m/s]"

    # ambient H2O [mol/mol] from RH
    esat, _ = e_sat(dat['Tair'])
    dat['H2O'] = (dat['RH'] / 100.0) * esat / dat['P']
    cols.append('H2O')
    readme += "\nH2O: Ambient H2O [mol/mol]"

    # ambient CO2 [ppm]
    readme += "\nCO2: Ambient CO2 [ppm]"
    if 'CO2' not in dat:
        dat['CO2'] = 400.0
        readme += " - set constant!"
    cols.append('CO2')

    # zenith angle
    jday = dat.index.dayofyear + dat.index.hour / 24.0 + dat.index.minute / 1440.0
# TEST (PERIOD START)
    jday = dat.index.dayofyear + dat.index.hour / 24.0 + dat.index.minute / 1440.0 + dt / 2.0 / 86400.0
    dat['Zen'], _, _, _, _, _ = solar_angles(lat, lon, jday, timezone=timezone)
    cols.append('Zen')
    readme += "\nZen: Zenith angle [rad], (lat = %.2f, lon = %.2f)" % (lat, lon)

    # radiation components

    if {'LWin','diffPar', 'dirPar', 'diffNir', 'dirNir'}.issubset(dat.columns) == False:
        f_cloud, f_diff, emi_sky = compute_clouds_rad(dat['doy'].values,
                                                      dat['Zen'].values,
                                                      dat['Rg'].values,
                                                      dat['H2O'].values * dat['P'].values,
                                                      dat['Tair'].values)

    if 'LWin' not in dat or dat['LWin'].isnull().any():
        if 'LWin' not in dat:
            dat['LWin']=np.nan
            print('Longwave radiation estimated')
        else:
            print('Longwave radiation partly estimated')
        # Downwelling longwve radiation
        # solar constant at top of atm.
        So = 1367
        # clear sky Global radiation at surface
        dat['Qclear'] = np.maximum(0.0,
                        (So * (1.0 + 0.033 * np.cos(2.0 * np.pi * (np.minimum(dat['doy'].values, 365) - 10) / 365)) * np.cos(dat['Zen'].values)))
        tau_atm = tau_atm = dat['Rg'].rolling(4,1).sum() / (dat['Qclear'].rolling(4,1).sum() + EPS)
        # cloud cover fraction
        dat['f_cloud'] = 1.0 - (tau_atm - 0.2) / (0.7 - 0.2)
        dat['f_cloud'][dat['Qclear'] < 10] = np.nan

        dat['Qclear_12h'] = dat['Qclear'].resample('12H').sum()
        dat['Qclear_12h'] = dat['Qclear_12h'].fillna(method='ffill')
        dat['Rg_12h'] = dat['Rg'].resample('12H').sum()
        dat['Rg_12h'] = dat['Rg_12h'].fillna(method='ffill')

        tau_atm = dat['Rg_12h'] / (dat['Qclear_12h'] + EPS)
        dat['f_cloud_12h'] = 1.0 - (tau_atm -0.2) / (0.7 - 0.2)

        dat['f_cloud'] = np.where((dat.index.hour > 12) & (dat['f_cloud_12h'] < 0.2), 0.0, dat['f_cloud'])
        dat['f_cloud'] = dat['f_cloud'].fillna(method='ffill')
        dat['f_cloud'] = dat['f_cloud'].fillna(method='bfill')
        dat['f_cloud'][dat['f_cloud'] < 0.0] = 0.0
        dat['f_cloud'][dat['f_cloud'] > 1.0] = 1.0

        emi0 = 1.24 * (dat['H2O'].values * dat['P'].values / 100 /(dat['Tair'].values + 273.15))**(1./7.)
        emi_sky = (1 - 0.84 * dat['f_cloud']) * emi0 + 0.84 * dat['f_cloud']

        # estimated long wave budget
        b = 5.6697e-8  # Stefan-Boltzman constant (W m-2 K-4)
        dat['LWin_estimated'] = emi_sky * b *(dat['Tair'] + 273.15)**4 # Wm-2 downwelling LW

        dat[['LWin','LWin_estimated']].plot(kind='line')

        dat['LWin'] = np.where(np.isfinite(dat['LWin']),dat['LWin'],dat['LWin_estimated'])

    cols.append('LWin')
    readme += "\nLWin: Downwelling long wave radiation [W/m2]"

    # Short wave radiation; separate direct and diffuse PAR & NIR
    if {'diffPar', 'dirPar', 'diffNir', 'dirNir'}.issubset(dat.columns) == False:
        print('Shortwave radiation components estimated')
        dat['diffPar'] = f_diff * fpar * dat['Rg']
        dat['dirPar'] = (1 - f_diff) * fpar * dat['Rg']
        dat['diffNir'] = f_diff * (1 - fpar) * dat['Rg']
        dat['dirNir'] = (1 - f_diff) * (1 - fpar) * dat['Rg']
    cols.extend(('diffPar', 'dirPar', 'diffNir', 'dirNir'))
    readme += "\ndiffPar: Diffuse PAR [W/m2] \ndirPar: Direct PAR [W/m2]"
    readme += "\ndiffNir: Diffuse NIR [W/m2] \ndirNir: Direct NIR [W/m2]"

    if {'Tsoil', 'Wliq'}.issubset(dat.columns):
        cols.extend(('Tsoil', 'Wliq'))
        dat['Wliq'] = dat['Wliq'] / 100.0
        readme += "\nTsoil: Soil surface layer temperature [degC]]"
        readme += "\nWliq: Soil surface layer moisture content [m3 m-3]"

    X = np.zeros(len(dat))
    DDsum = np.zeros(len(dat))
    for k in range(1,len(dat)):
        if dat['doy'][k] != dat['doy'][k-1]:
            X[k] = X[k - 1] + 1.0 / 8.33 * (dat['Tdaily'][k-1] - X[k - 1])
            if dat['doy'][k] == 1:  # reset in the beginning of the year
                DDsum[k] = 0.
            else:
                DDsum[k] = DDsum[k - 1] + max(0.0, dat['Tdaily'][k-1] - 5.0)
        else:
            X[k] = X[k - 1]
            DDsum[k] = DDsum[k - 1]
    dat['X'] = X
    cols.append('X')
    readme += "\nX: phenomodel delayed temperature [degC]"
    dat['DDsum'] = DDsum
    cols.append('DDsum')
    readme += "\nDDsum: degreedays [days]"

    dat = dat[cols]
    dat[cols].plot(subplots=True, kind='line')

    print("NaN values in forcing data:")
    print(dat.isnull().any())

    save_df_to_csv(dat, output_file, readme=readme,fp=dir_save)

def save_df_to_csv(df, fn, readme='', fp="forcing/", timezone = +2):
    """
    Save dataframe with datetime index to csv file with corresponding readme.txt.
    Args:
        df (DataFrame): data to save
        fn (str): filename of saved file (.csv not included)
        readme (str): readme corresponding to df to save as txt file
        fp (str): filepath, forcing folder used as default
        timezone (float): time zone in refernce to UTC, default UTC + 2
    """

    # add datetime as columns
    df.insert(0, 'yyyy', df.index.year.values)
    df.insert(1, 'mo', df.index.month.values)
    df.insert(2, 'dd', df.index.day.values)
    df.insert(3, 'hh', df.index.hour.values)
    df.insert(4, 'mm', df.index.minute.values)

    df.to_csv(path_or_buf=fp + fn + ".csv", sep=',', na_rep='NaN', index=False)
    Readme = "Readme for " + fn + ".csv"
    Readme += "\n\nKersti Haahti, Luke " + str(datetime.datetime.now().date())
    Readme += "\n\nyyyy, mo, dd, hh, mm: datetime [UTC + %.1f]" % timezone
    Readme += readme
    outF = open(fp + fn + "_readme.txt", "w")
    print(Readme, file=outF)
    outF.close()

