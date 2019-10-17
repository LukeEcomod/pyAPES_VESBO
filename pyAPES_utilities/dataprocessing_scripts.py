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

direc = "C:/Users/L1656/Documents/Git_repos/pyAPES_Kersti/"

"""

### LETTOSUO ###
lettosuo_data = read_lettosuo_data()
gap_fill_lettosuo_meteo(lettosuo_data)
create_forcingfile("Lettosuo_meteo_2010_2018", "Lettosuo_forcing_2010_2018",
                   lat=60.63, lon=23.95, P_unit = 1e2) # [hPa]

### HYYTIALA ###
gather_hyde_data()
create_forcingfile("Hyde_data_1997_2016", "Hyde_forcing_1997_2016",
                   lat=61.51, lon=24.0, P_unit = 1e3) # [kPa]

"""

def create_forcingfile(meteo_file, output_file, lat, lon, P_unit, timezone=+2.0):
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

    forc_fp = direc + "forcing/" + meteo_file +".csv"
    dat = pd.read_csv(forc_fp, sep=',', header='infer', encoding = 'ISO-8859-1')

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

    if 'LWin' not in dat:
        print('Longwave radiation estimated')
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
        dat['LWin'] = emi_sky * b *(dat['Tair'] + 273.15)**4 # Wm-2 downwelling LW

        dat[['LWin','f_cloud']].plot(subplots=True, kind='line')

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

    save_df_to_csv(dat, output_file, readme=readme,fp=direc + "forcing/")

def read_lettosuo_data(starttime='09-01-2009',endtime='01-01-2019'):

    """
    Reads data related to lettosuo case to dataframe.
    """

    # filepaths
    forc_fp = ["O:/Projects/Lettosuo/Forcing_data/Annalea2/Letto1_EC.csv",  # period end
               "O:/Projects/Lettosuo/Forcing_data/Annalea1/Letto1_meteo.csv",
               "O:/Projects/Lettosuo/Forcing_data/Annalea1/Letto1_metsanpohja.csv",
               "O:/Projects/Lettosuo/Forcing_data/FMI/jokioinen_meteo.txt",
               "O:/Projects/Lettosuo/Forcing_data/FMI/jokioinen_rad.txt",
               "O:/Projects/Lettosuo/Forcing_data/Annalea2/Letto1_meteo_gapfilled.csv",
               "O:/Projects/Lettosuo/Forcing_data/FMI/jokioinen_prec1.txt",
               "O:/Projects/Lettosuo/Forcing_data/FMI/jokioinen_prec2.txt",
               "O:/Projects/Lettosuo/Forcing_data/FMI/somero_meteo.txt",
               "O:/Projects/Lettosuo/Forcing_data/FMI/salo_kiikala_meteo.txt",
               "O:/Projects/Lettosuo/Forcing_data/MikaK/osittaishakkuu_NEE_GPP_RECO.csv",  # period end
               "O:/Projects/Lettosuo/Forcing_data/MikaK/meteo_concat.csv",
               "O:/Projects/Lettosuo/Forcing_data/MikaK/metsanpohja_concat.csv",
               "O:/Projects/Lettosuo/Forcing_data/FMI/jokioinen_daily_prec.txt"]

    index=pd.date_range(starttime, endtime, freq='0.5H')
    lettosuo_data=pd.DataFrame(index=index, columns=[])

    for fp in forc_fp:
        dat = pd.read_csv(fp, sep=',', header='infer', encoding = 'ISO-8859-1')
        if fp.split("/")[-2] == 'Annalea2':
            dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
            # period end
            dat.index = dat.index - pd.Timedelta(hours=0.5)
        else:
            dat.index = pd.to_datetime(dat.ix[:,0], yearfirst=True)
        if fp=="O:/Projects/Lettosuo/Forcing_data/MikaK/osittaishakkuu_NEE_GPP_RECO.csv":
            # period end
            dat.index = dat.index - pd.Timedelta(hours=0.5)
        if fp.split("/")[-2] == 'FMI':
            if fp=="O:/Projects/Lettosuo/Forcing_data/FMI/jokioinen_daily_prec.txt":
                # Vuorokauden sademäärä (mm) kuvaa vuorokauden aikana (aamun 06 UTC
                # havainnosta seuraavan vuorokauden 06 UTC havaintoon) kertynyttä sateen määrää.
                dat.index = dat.index + pd.Timedelta(hours=6.0)
            # UTC -> UTC + 2
            dat.index = dat.index + pd.Timedelta(hours=2)
        dat.index=dat.index.map(lambda x: x.replace(second=0))
        dat.ix[:,0]=dat.index
        dat = dat.drop_duplicates(subset=dat.columns[0])
        dat = dat.drop(dat.columns[0], axis=1)
        dat.columns = fp.split("/")[-1].split(".")[0] + ': ' + dat.columns
        dat = dat[(dat.index >= starttime) & (dat.index < endtime)]
        if len(np.setdiff1d(dat.index, index)) > 0:
            print(fp, np.setdiff1d(dat.index, index))
            raise ValueError("Error")
        lettosuo_data=lettosuo_data.merge(dat, how='outer', left_index=True, right_index=True)


    # divide hourly precipitation to half hour
    cols=['somero_meteo: Precipitation amount',
          'jokioinen_prec1: Prec [mm h-1]',
          'jokioinen_meteo: Precipitation amount',
          'jokioinen_prec2: Prec [mm h-1]']
    for col in cols:
        lettosuo_data[col]=lettosuo_data[col].replace(-1,0)
        lettosuo_data[col][1:-1:2]=lettosuo_data[col][0:-2:2].values
        lettosuo_data[col]=lettosuo_data[col].values/2.0

    lettosuo_data['jokioinen_daily_prec: Precipitation amount'] = lettosuo_data['jokioinen_daily_prec: Precipitation amount'].replace(-1,0)

    lettosuo_data = lettosuo_data[(lettosuo_data.index >= starttime) &
                                  (lettosuo_data.index < endtime)]

    return lettosuo_data

def read_Svb_data(forc_fp=None):

    """
    Reads data related to Svartberget to dataframe.
    """

    if forc_fp == None: # filepaths
        forc_fp = ["O:/Projects/Muut projektit/Natalia/Svartberget_data/SE-Svb_eco/concat.csv",
                   "O:/Projects/Muut projektit/Natalia/Svartberget_data/SE-Svb_fluxes/concat.csv",
                   "O:/Projects/Muut projektit/Natalia/Svartberget_data/SE-Svb_meteo/concat.csv",
                   "O:/Projects/Muut projektit/Natalia/Svartberget_data/SE-Svb_profile/concat.csv",
                   "O:/Projects/Muut projektit/Natalia/Svartberget_data/SE-Svb_T-profile/concat.csv"]

    index=pd.date_range('01-01-2014','01-01-2017',freq='0.5H')
    data=pd.DataFrame(index=index, columns=[])

    for fp in forc_fp:
        dat = pd.read_csv(fp, sep=',', header='infer', encoding = 'ISO-8859-1')
        dat.index = pd.to_datetime(dat.ix[:,0] + ' ' + dat.ix[:,1], dayfirst=True)
        dat.index=dat.index.map(lambda x: x.replace(second=0))
        dat.ix[:,0]=dat.index
        dat = dat.drop_duplicates(subset=dat.columns[0])
        dat = dat.drop([dat.columns[0], dat.columns[1]], axis=1)
        dat.columns = fp.split("/")[-2] + ': ' + dat.columns
        if len(np.setdiff1d(dat.index, index)) > 0:
            print(fp, np.setdiff1d(dat.index, index))
            raise ValueError("Error")
        data=data.merge(dat, how='outer', left_index=True, right_index=True)

    return data

def gap_fill_lettosuo_meteo(lettosuo_data, plot=False,starttime='09-01-2009',endtime='01-01-2019'):
    """
    Gap fills Lettosuo meteo (and collesponding flags)
    and save to file (with readme)
    """

    frames = []
    readme = ""

    # --- Precipitation ---

    # Jokioinen
    df, info = fill_gaps(lettosuo_data[['jokioinen_prec1: Prec [mm h-1]',
                                        'jokioinen_meteo: Precipitation amount',
                                        'jokioinen_prec2: Prec [mm h-1]',
                                        'somero_meteo: Precipitation amount']],
                         'Prec_ref', 'Jokioinen gapfilled precipitaion [mm/30min]',
                         fill_nan=0.0, plot=plot)

    Prec_daily = lettosuo_data[['jokioinen_daily_prec: Precipitation amount']].resample('24h', base=8, convention='end').mean()
    Prec_daily = Prec_daily.merge(df[['Prec_ref']].resample('24h', base=8).sum(), how='outer', left_index=True, right_index=True)
    Prec_daily['correction'] = Prec_daily['jokioinen_daily_prec: Precipitation amount']/Prec_daily['Prec_ref']
    Prec_daily['correction'] = Prec_daily['correction'].fillna(1.0)

    lettosuo_data[['jokioinen_daily_prec: Prec']] = lettosuo_data[['jokioinen_daily_prec: Precipitation amount']].fillna(method='ffill')/48.0
    lettosuo_data = lettosuo_data.merge(Prec_daily[['correction']], how='outer', left_index=True, right_index=True)
    lettosuo_data['correction'] = lettosuo_data['correction'].fillna(method='ffill')
    lettosuo_data = lettosuo_data[(lettosuo_data.index >= starttime) &
                                  (lettosuo_data.index < endtime)]
    df['Prec_ref'] = np.where(np.isinf(lettosuo_data['correction']),
                              lettosuo_data['jokioinen_daily_prec: Prec'],
                              lettosuo_data['correction']*df['Prec_ref'])

    lettosuo_data['jokioinen_prec: Prec'] = df['Prec_ref']
    frames.append(df)
    readme += info

    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm))'] = np.where(
            lettosuo_data.index < '01-01-2018',
            lettosuo_data['Letto1_metsanpohja: avg(Rain (mm))'],
            30.0 * lettosuo_data['metsanpohja_concat: avg(Rain (?))'])

    # Lettosuo, old
    # consider prec data unrealiable when Tair < 2C
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!']=np.where(
            lettosuo_data['Letto1_metsanpohja: avg(Temp (C))'] < 2.0,
            np.nan, lettosuo_data['Letto1_metsanpohja: avg(Rain (mm))'])

    # if rolling 3 day mean prec less than 10% of jokionen rolling mean, remove
    Prec_daily_ref = df['Prec_ref'].rolling(3*48, 1).sum()
    Prec_daily_let = lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'].fillna(0).rolling(3*48, 1).sum()
#    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!']=np.where(
#            Prec_daily_let < 0.1 * Prec_daily_ref,
#            np.nan, lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'])

    # clearly no rain event at lettosuo (wtd data)
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'][
            (lettosuo_data.index > '07-03-2011') & (lettosuo_data.index < '07-11-2011')]=0.0

    # problematic sections
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'][
            (lettosuo_data.index > '01-01-2016') & (lettosuo_data.index < '07-17-2016')]=np.nan
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'][
            (lettosuo_data.index > '01-01-2018') & (lettosuo_data.index < '07-10-2018')]=np.nan

    df, info = fill_gaps(lettosuo_data[['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!',
                                        'somero_meteo: Precipitation amount',
                                        'jokioinen_prec: Prec']],
                         'Prec', 'Lettosuo gapfilled precipitaion [mm/30min]',
                         fill_nan=0.0, plot=plot)

    frames.append(df)
    readme += info

    # Lettosuo
    # consider prec data unrealiable when Tair < 2C
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!']=np.where(
            lettosuo_data['Letto1_metsanpohja: avg(Temp (C))'] < 2.0,
            np.nan, lettosuo_data['Letto1_metsanpohja: avg(Rain (mm))'])

    # if rolling 3 day mean prec less than 10% of jokionen rolling mean, remove
    Prec_daily_let = lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'].fillna(0).rolling(3*48, 1).sum()
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!']=np.where(
            Prec_daily_let < 0.1 * Prec_daily_ref,
            np.nan, lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'])

    # clearly no rain event at lettosuo (wtd data)
    lettosuo_data['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!'][
            (lettosuo_data.index > '07-03-2011') & (lettosuo_data.index < '07-11-2011')]=0.0

    df, info = fill_gaps(lettosuo_data[['Letto1_metsanpohja: avg(Rain (mm)) !sections removed!',
                                        'somero_meteo: Precipitation amount',
                                        'jokioinen_prec: Prec']],
                         'Prec_old', 'Lettosuo gapfilled precipitaion [mm/30min]',
                         fill_nan=0.0, plot=plot)
    frames.append(df)
    readme += info

    # Somero
    df, info = fill_gaps(lettosuo_data[['somero_meteo: Precipitation amount',
                                        'jokioinen_prec: Prec']],
                         'Prec_ref2', 'Somero gapfilled precipitaion [mm/30min]',
                         fill_nan=0.0, plot=plot)
    frames.append(df)
    readme += info

    # --- Air temperature ---

    lettosuo_data['Letto1_meteo: avg(Temp (C))'] = np.where(
            lettosuo_data.index < '01-01-2018',
            lettosuo_data['Letto1_meteo: avg(Temp (C))'],
            lettosuo_data['meteo_concat: avg(Temp (C))'])

    df, info = fill_gaps(lettosuo_data[['Letto1_meteo: avg(Temp (C))',
                                        'somero_meteo: Air temperature',
                                        'jokioinen_meteo: Air temperature',
                                        'Letto1_meteo_gapfilled: PaikAirT T']],
                         'Tair', 'Air temperature [degC]', fill_nan='linear', plot=plot)
    frames.append(df)
    readme += info

    # --- Relative humidity ---
    lettosuo_data['Letto1_meteo: avg(RH (%))'] = np.where(
            lettosuo_data.index < '01-01-2018',
            lettosuo_data['Letto1_meteo: avg(RH (%))'],
            lettosuo_data['meteo_concat: avg(RH (%))'])

    df, info = fill_gaps(lettosuo_data[['Letto1_meteo: avg(RH (%))',
                                        'somero_meteo: Relative humidity',
                                        'jokioinen_meteo: Relative humidity',
                                        'Letto1_meteo_gapfilled: Paik RH']],
                         'RH', 'Relative humidity [%]', fill_nan='linear', plot=plot)
    # Check for values greater than 100%
    df['RH'][df['RH'] > 100.0] = 100.0
    frames.append(df)
    readme += info

    # --- Global radiation ---
    lettosuo_data['Letto1_meteo: avg(Glob (W/m2))'] = np.where(
            lettosuo_data.index < '01-01-2018',
            lettosuo_data['Letto1_meteo: avg(Glob (W/m2))'],
            lettosuo_data['meteo_concat: avg(Glob (W/m2))'])

    df, info = fill_gaps(lettosuo_data[['Letto1_meteo: avg(Glob (W/m2))',
                                        'Letto1_meteo_gapfilled: PaikGlob2',
                                        'jokioinen_rad: Global radiation']],
                         'Rg', 'Global radiation [W/m2]', fill_nan='linear', plot=plot)
    df['Rg'][df['Rg'] < 0.0] = 0.0
    frames.append(df)
    readme += info

    # --- Wind speed ---
    lettosuo_data['Letto1_EC: wind speed (m/s)'] = np.where(
            lettosuo_data.index < '01-01-2016',
            lettosuo_data['Letto1_EC: wind speed (m/s)'],
            lettosuo_data['osittaishakkuu_NEE_GPP_RECO: wind_speed [m/s]'])

    lettosuo_data['Letto1_EC: wind speed (m/s) !u > 10 removed!']=np.where(
            lettosuo_data['Letto1_EC: wind speed (m/s)'] > 10.0,
            np.nan, lettosuo_data['Letto1_EC: wind speed (m/s)'])
    lettosuo_data['salo_kiikala_meteo: Wind speed'][lettosuo_data['salo_kiikala_meteo: Wind speed'] == 0.0] = np.nan
    lettosuo_data['Derived from salo_kiikala_meteo: wind speed (0.57*U_ref + 0.55)'] = 0.57 * lettosuo_data['salo_kiikala_meteo: Wind speed'] + 0.55
    df, info = fill_gaps(lettosuo_data[['Letto1_EC: wind speed (m/s) !u > 10 removed!',
                                        'somero_meteo: Wind speed',
                                        'Derived from salo_kiikala_meteo: wind speed (0.57*U_ref + 0.55)']],
                         'U', 'Wind speed [m/s]', fill_nan='linear', plot=plot)
    frames.append(df)
    readme += info

    # --- Friction velocity ---
    lettosuo_data['Letto1_EC: friction velocity (m/s)'] = np.where(
            lettosuo_data.index < '01-01-2016',
            lettosuo_data['Letto1_EC: friction velocity (m/s)'],
            lettosuo_data['osittaishakkuu_NEE_GPP_RECO: ustar [m s-1]'])

    lettosuo_data['Ustar = 0.2 * U'] = 0.2 * df['U']
    df, info = fill_gaps(lettosuo_data[['Letto1_EC: friction velocity (m/s)',
                                        'Ustar = 0.2 * U']],
                         'Ustar', 'Friction velocity [m/s]', fill_nan='linear', plot=plot)
    frames.append(df)
    readme += info

    # --- Ambient pressure ---
    lettosuo_data['Letto1_meteo: avg(Press (hPa))'] = np.where(
            lettosuo_data.index < '01-01-2018',
            lettosuo_data['Letto1_meteo: avg(Press (hPa))'],
            lettosuo_data['meteo_concat: avg(Press (hPa))'])

    lettosuo_data['Derived from salo_kiikala_meteo: Pressure (msl) (P_ref - 16.9)'] = lettosuo_data['salo_kiikala_meteo: Pressure (msl)'] - 16.9
    df, info = fill_gaps(lettosuo_data[['Letto1_meteo: avg(Press (hPa))',
                                        'Derived from salo_kiikala_meteo: Pressure (msl) (P_ref - 16.9)']],
                         'P', 'Ambient pressure [hPa]', fill_nan='linear', plot=plot)
    frames.append(df)
    readme += info

    # -- Snow ---
    # Jokioinen
    df, info = fill_gaps(lettosuo_data[['jokioinen_meteo: Snow depth']],
                         'Snow_depth1', 'Jokioinen snow depth [cm]',
                         fill_nan=np.nan, plot=plot)
    frames.append(df)
    readme += info
    # Somero
    df, info = fill_gaps(lettosuo_data[['somero_meteo: Snow depth']],
                         'Snow_depth2', 'Somero snow depth [cm]',
                         fill_nan=np.nan, plot=plot)
    frames.append(df)
    readme += info
    # Salo kiikala
    df, info = fill_gaps(lettosuo_data[['salo_kiikala_meteo: Snow depth']],
                         'Snow_depth3', 'Salo-kiikala snow depth [cm]',
                         fill_nan=np.nan, plot=plot)
    frames.append(df)
    readme += info

    letto_data=pd.concat(frames, axis=1)
    letto_data[['Prec_ref', 'Prec', 'Tair', 'Rg', 'U', 'Ustar', 'RH', 'P','Snow_depth1','Snow_depth2','Snow_depth3']].plot(subplots=True,kind='line')

    cumPrec=letto_data[['Prec_ref', 'Prec_ref2', 'Prec_old', 'Prec']].groupby(df.index.year).cumsum()
    cumPrec.plot()

    save_df_to_csv(letto_data, "Lettosuo_meteo_2010_2018", readme=readme, fp=direc + "forcing/")

def gather_hyde_data():
    """
    Collects yearly hyde data to one file that is saved.
    """

    # read files in format Forcing_YYYY.dat
    directory = "O:/Projects/Samulilta/Hyde/Forcing_corrected/"
    frames = []
    columns =['U','ust','Ta','RH','CO2','H2O','Prec','P',
              'dirPar','diffPar','dirNir','diffNir','Rnet',
              'LWin','LWout','LWnet','Tsh','Tsa','Tsc','Wh','Wa']
    for year in range(1997,2017):
        forc_fp = directory + "Forcing_" + str(year) + ".dat"
        dat = pd.read_csv(forc_fp, sep=',', header='infer', encoding = 'ISO-8859-1')
        index = pd.date_range('01-01-' + str(year),'31-12-' + str(year) + ' 23:59:59',freq='0.5H')
        if len(index) != len(dat):
            print("Year " + str(year) + ": Missing values!")
        else:
            dat.index = index
            frames.append(dat[columns])

    hyde_data=pd.concat(frames)

    # read files in format FIHy_YYYY_pd.csv
    directory = "O:/Projects/Samulilta/Hyde/FIHy_1997_2016_new/"
    frames = []
    columns =['NEE','GPP','LE','ET','fRg']
    for year in range(1997,2017):
        forc_fp = directory + "FIHy_" + str(year) + "_pd.csv"
        dat = pd.read_csv(forc_fp, sep=',', header='infer', encoding = 'ISO-8859-1')
        index = pd.date_range('01-01-' + str(year),'31-12-' + str(year) + ' 23:59:59',freq='0.5H')
        if len(index) != len(dat):
            print("Year " + str(year) + ": Missing values!")
        else:
            dat.index = index
            frames.append(dat[columns])

    df = pd.concat(frames)
    hyde_data=hyde_data.merge(df, how='outer', left_index=True, right_index=True)
    hyde_data=hyde_data.rename(columns={'Ta':'Tair', 'ust':'Ustar', 'fRg':'Rg'})

    save_df_to_csv(hyde_data, "Hyde_data_1997_2016", fp=direc + "forcing/")

def read_lettosuo_EC():

    """
    Reads data related to lettosuo case to dataframe.
    """

    # filepaths
    forc_fp = ["O:/Projects/Lettosuo/Forcing_data/Annalea1/avohakkuu_EC.csv",
               "O:/Projects/Lettosuo/Forcing_data/Annalea2/Letto1_EC.csv",
               "O:/Projects/Lettosuo/Forcing_data/MikaK/Partial_EC_gapfilled_fluxes.csv",
               "O:/Projects/Lettosuo/Forcing_data/Annalea2/energyfluxes_lettosuo.csv"]

    index=pd.date_range('01-01-2009','06-01-2018',freq='0.5H')
    lettosuo_data=pd.DataFrame(index=index, columns=[])

    for fp in forc_fp:
        dat = pd.read_csv(fp, sep=',', header='infer', encoding = 'ISO-8859-1')
        if fp.split("/")[-2] == 'Annalea2':
            dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
            # period end
            dat.index = dat.index - pd.Timedelta(hours=0.5)
        else:
            dat.index = pd.to_datetime(dat.ix[:,0], yearfirst=True)
        if fp.split("/")[-2] == 'FMI':
            # UTC -> UTC + 2
            dat.index = dat.index + pd.Timedelta(hours=2)
        dat.index=dat.index.map(lambda x: x.replace(second=0))
        dat.ix[:,0]=dat.index
        dat = dat.drop_duplicates(subset=dat.columns[0])
        dat = dat.drop(dat.columns[0], axis=1)
        dat.columns = fp.split("/")[-1].split(".")[0] + ': ' + dat.columns
        lettosuo_data=lettosuo_data.merge(dat, how='outer', left_index=True, right_index=True)

    lettosuo_data = lettosuo_data[(lettosuo_data.index >= '01-01-2010') &
                                  (lettosuo_data.index <= '01-01-2018')]

    lettosuo_EC = lettosuo_data[
            ['Letto1_EC: PaikNEE_S',
             'Letto1_EC: PaikNEE_N',
             'Letto1_EC: GPP_S',
             'Letto1_EC: GPP_N',
             'Letto1_EC: Reco_S',
             'Letto1_EC: Reco_N',
             'Letto1_EC: LE (Wm-2)',
             'Letto1_EC: SH (W m-2)',
             'Partial_EC_gapfilled_fluxes: NEE [mg CO2 m-2 s-1]',
             'Partial_EC_gapfilled_fluxes: GPP [mg CO2 m-2 s-1]',
             'Partial_EC_gapfilled_fluxes: Reco [mg CO2 m-2 s-1]',
             'Partial_EC_gapfilled_fluxes: gapped',
             'energyfluxes_lettosuo: LE [W m-2]',
             'energyfluxes_lettosuo: SH [W m-2]',
             'avohakkuu_EC: NEE [mg CO2 m-2 s-1]',
             'avohakkuu_EC: GPP [mg CO2 m-2 s-1]',
             'avohakkuu_EC: Reco [mg CO2 m-2 s-1]',
             'avohakkuu_EC: gapped',
             'avohakkuu_EC: latent heat flux: FL [W m-2]',
             'avohakkuu_EC: sensible heat flux: FH [W m-2]']].copy()

    lettosuo_EC = lettosuo_EC.rename(columns={
             'Letto1_EC: PaikNEE_S': 'control-S_NEE',
             'Letto1_EC: GPP_S': 'control-S_GPP',
             'Letto1_EC: Reco_S': 'control-S_Reco',
             'Letto1_EC: PaikNEE_N': 'control-N_NEE',
             'Letto1_EC: GPP_N': 'control-N_GPP',
             'Letto1_EC: Reco_N': 'control-N_Reco',
             'Letto1_EC: LE (Wm-2)': 'control-N_LE',
             'Letto1_EC: SH (W m-2)': 'control-N_SH',
             'Partial_EC_gapfilled_fluxes: NEE [mg CO2 m-2 s-1]': 'partial_NEE',
             'Partial_EC_gapfilled_fluxes: GPP [mg CO2 m-2 s-1]': 'partial_GPP',
             'Partial_EC_gapfilled_fluxes: Reco [mg CO2 m-2 s-1]': 'partial_Reco',
             'Partial_EC_gapfilled_fluxes: gapped': 'partial_gapped',
             'energyfluxes_lettosuo: LE [W m-2]': 'partial_LE',
             'energyfluxes_lettosuo: SH [W m-2]': 'partial_SH',
             'avohakkuu_EC: NEE [mg CO2 m-2 s-1]': 'clearcut_NEE',
             'avohakkuu_EC: GPP [mg CO2 m-2 s-1]': 'clearcut_GPP',
             'avohakkuu_EC: Reco [mg CO2 m-2 s-1]': 'clearcut_Reco',
             'avohakkuu_EC: gapped': 'clearcut_gapped',
             'avohakkuu_EC: latent heat flux: FL [W m-2]': 'clearcut_LE',
             'avohakkuu_EC: sensible heat flux: FH [W m-2]': 'clearcut_SH'})

    lettosuo_EC['control-S_LE'] = lettosuo_EC['control-N_LE'].values
    lettosuo_EC['control-S_SH'] = lettosuo_EC['control-N_SH'].values
    lettosuo_EC['control-S_gapped'] = np.where(
            np.isfinite(lettosuo_data['Letto1_EC: NEE_South'].values), 0, 1)
    lettosuo_EC['control-N_gapped'] = np.where(
            np.isfinite(lettosuo_data['Letto1_EC: NEE_North'].values), 0, 1)

    for column in lettosuo_EC.columns:
        if column.split('_')[-1] == 'LE':
            lettosuo_EC[column][lettosuo_EC[column] < -100] = lettosuo_EC[column][lettosuo_EC[column] < -100] * np.nan
            lettosuo_EC[column][lettosuo_EC[column] > 1000] = lettosuo_EC[column][lettosuo_EC[column] > 1000] * np.nan

    save_df_to_csv(lettosuo_EC, "Lettosuo_EC", fp=direc + "forcing/")

def gather_data(dir_path="O:/Projects/Lettosuo/Forcing_data/datat/meteo/", cols=None, sort=0):
    """
    Collect files in one directory to one file.
    """

    filenames = listdir(dir_path)

    frames = []

    for fn in filenames:
        if fn != 'concat.csv':
            print(fn)
            dat = pd.read_csv(dir_path + fn, sep=',', header='infer', encoding = 'ISO-8859-1')
            if cols is not None:
                frames.append(dat[cols])
            else:
                frames.append(dat)

    data = pd.concat(frames, ignore_index=True, sort=False)
    data = data[data[data.columns[0]] != data.columns[0]]
    data = data[data[data.columns[0]] != 'dd/mm/yyyy']
    if sort == 0: # datetime in first column
        data = data.sort_values(by=data.columns[0])
    elif sort == 1: # date in first column, time in second
        data = data.sort_values(by=[data.columns[0], data.columns[1]])

    data.to_csv(path_or_buf=dir_path + 'concat.csv', sep=',', na_rep='NaN', index=False)

def rad_to_30min():
    fd="O:/Projects/Lettosuo/Forcing_data/FMI/"
    fn="jokioinen_rad5min.txt"
    dat = pd.read_csv(fd+fn, sep=',', header='infer', encoding = 'ISO-8859-1')
    dat.index=pd.to_datetime(dat.ix[:,0], yearfirst=True)
    dat2=dat.resample('30T').mean()
    dat2.to_csv(path_or_buf=fd + "jokioinen_rad.txt", sep=',', na_rep='NaN', index=True)

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

