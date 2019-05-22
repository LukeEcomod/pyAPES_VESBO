# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 09:23:39 2019

@author: L1656
"""

import pandas as pd
import numpy as np
import datetime
from pyAPES_utilities.plotting import plot_columns
import matplotlib.pyplot as plt
from pyAPES_utilities.dataprocessing_scripts import save_df_to_csv

def read_WTD_data():

    fpaths = [r"H:\Lettosuo\WTD_paavolta\Lettosuo_WTD_EC.csv",
               r"H:\Lettosuo\WTD_paavolta\Lettosuo_WTD_trans.csv"]

    index=pd.date_range('01-01-2010','06-01-2019',freq='1H')
    data=pd.DataFrame(index=index, columns=[])

    for fp in fpaths:
        dat = pd.read_csv(fp, sep=';', header='infer')
        dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
        dat.index = dat.index + datetime.timedelta(hours=0.5)
        dat.index = dat.index.map(lambda x: x.replace(minute=0, second=0))
        dat.ix[:,0]=dat.index
        dat = dat.drop_duplicates(subset=dat.columns[0])
        dat = dat.drop([dat.columns[0]], axis=1)
        if len(np.setdiff1d(dat.index, index)) > 0:
            print(fp, np.setdiff1d(dat.index, index))
            raise ValueError("Error")
        data=data.merge(dat, how='outer', left_index=True, right_index=True)

    return data

def process_WTD():
    # read data
    WTD=read_WTD_data()

    # putken pohja tulee vastaan
    WTD['ctrl8m'][WTD['ctrl8m'] < -83] = np.nan
    WTD['ctrl22.5m'][WTD['ctrl22.5m'] < -81] = np.nan

    # pitkäaikaiset loggerit
    plot_columns(WTD[['WT_E','WT_N','WT_S','WT_W']])  # -> korrelaatio 'WT_E','WT_S','WT_W' välillä > 0.95
    # keskiarvo ja keskihajonta pitkäaikaisista
    WTD['WT_ESW'] = np.nanmean([WTD['WT_E'].values,
                                WTD['WT_S'].values,
                                WTD['WT_W'].values], axis=0)
    WTD['WT_ESW_std'] = np.nanstd([WTD['WT_E'].values,
                                   WTD['WT_S'].values,
                                   WTD['WT_W'].values], axis=0)

    # ennen hakkuuta control, hakkuun jälkeen partial
    WTD['WT_ESW_part'] = np.where(WTD.index <= '03-15-2016', np.nan, WTD['WT_ESW'])
    WTD['WT_ESW_ctrl'] = np.where(WTD.index <= '03-15-2016', WTD['WT_ESW'], np.nan)

    # epämääräinen alkujakso avohakkuulla
    WTD['clear4m'][WTD.index < '12-01-2015'] = np.nan
    WTD['clear8m'][WTD.index < '12-01-2015'] = np.nan
    WTD['clear12m'][WTD.index < '12-01-2015'] = np.nan
    WTD['clear22.5m'][WTD.index < '12-01-2015'] = np.nan

    # kalibrointikausi hakkuun ajankohtaan asti
    WTD_calib = WTD[(WTD.index <= '03-15-2016')]

    # Tasokorjaus (slope=1) kalibrointi jakson perusteella, ehtona R2 > 0.7

    #plot_columns(WTD_calib[['ctrl4m','ctrl8m','ctrl12m','ctrl22.5m','WT_ESW']])
    #WTD['pred_ctrl8m'] = 0.99*WTD['ctrl8m'] + 11.68
    #WTD['pred_ctrl12m'] = 0.98*WTD['ctrl12m'] + 2.55
    #WTD['pred_ctrl22.5m'] = 0.80*WTD['ctrl22.5m'] + 1.53
    plot_columns(WTD_calib[['ctrl4m','ctrl8m','ctrl12m','ctrl22.5m','WT_ESW']],slope=1.0)
    WTD['pred_ctrl8m'] = WTD['ctrl8m'] + 12.33
    WTD['pred_ctrl12m'] = WTD['ctrl12m'] + 3.2
    WTD['pred_ctrl22.5m'] = WTD['ctrl22.5m'] + 11.27
    #WTD[['pred_ctrl8m','pred_ctrl12m','pred_ctrl22.5m','WT_ESW']].plot()

    #plot_columns(WTD_calib[['part4m','part8m','part12m','part22.5m','WT_ESW']],slope=1.0)
    #WTD['pred_part12m'] = 1.0*WTD['part12m'] + 10.84
    #WTD['pred_part22.5m'] = 0.85*WTD['part22.5m'] - 3.46
    plot_columns(WTD_calib[['part4m','part8m','part12m','part22.5m','WT_ESW']],slope=1.0)
    WTD['pred_part12m'] = WTD['part12m'] + 11.02
    WTD['pred_part22.5m'] = WTD['part22.5m'] + 2.51
    #WTD[['pred_part12m','pred_part22.5m','WT_ESW']].plot()

    #plot_columns(WTD_calib[['clear4m','clear8m','clear12m','clear22.5m','WT_ESW']])
    #WTD['pred_clear4m'] = 1.5*WTD['clear4m'] + 29.89
    #WTD['pred_clear12m'] = 1.32*WTD['clear12m'] + 15.85
    plot_columns(WTD_calib[['clear4m','clear8m','clear12m','clear22.5m','WT_ESW']],slope=1.0)
    WTD['pred_clear4m'] = WTD['clear4m'] + 8.69
    WTD['pred_clear12m'] = WTD['clear12m'] + 3.92
    #WTD[['pred_clear4m','pred_clear12m','WT_ESW']].plot()

    # keskiarvot ja vaihteluväli (keskihajonnan ja aikasarjojen vaihtelusta) käsittelyille

    WTD['control'] = 0.01 * np.nanmean([WTD['WT_ESW_ctrl'].values,
                                 WTD['pred_ctrl8m'].values,
                                 WTD['pred_ctrl12m'].values,
                                 WTD['pred_ctrl22.5m'].values], axis=0)
    WTD['control_min'] = 0.01 * np.nanmin([WTD['WT_ESW_ctrl'].values - WTD['WT_ESW_std'].values,
                                    WTD['pred_ctrl8m'].values - WTD['WT_ESW_std'].values,
                                    WTD['pred_ctrl12m'].values - WTD['WT_ESW_std'].values,
                                    WTD['pred_ctrl22.5m'].values - WTD['WT_ESW_std'].values], axis=0)
    WTD['control_max'] = 0.01 * np.nanmax([WTD['WT_ESW_ctrl'].values + WTD['WT_ESW_std'].values,
                                    WTD['pred_ctrl8m'].values + WTD['WT_ESW_std'].values,
                                    WTD['pred_ctrl12m'].values + WTD['WT_ESW_std'].values,
                                    WTD['pred_ctrl22.5m'].values + WTD['WT_ESW_std'].values], axis=0)

    WTD['partial'] = 0.01 * np.nanmean([WTD['WT_ESW_part'].values,
                                 WTD['pred_part12m'].values,
                                 WTD['pred_part22.5m'].values], axis=0)
    WTD['partial_min'] = 0.01 * np.nanmin([WTD['WT_ESW_part'].values - WTD['WT_ESW_std'].values,
                                    WTD['pred_part12m'].values - WTD['WT_ESW_std'].values,
                                    WTD['pred_part22.5m'].values - WTD['WT_ESW_std'].values], axis=0)
    WTD['partial_max'] = 0.01 * np.nanmax([WTD['WT_ESW_part'].values + WTD['WT_ESW_std'].values,
                                    WTD['pred_part12m'].values + WTD['WT_ESW_std'].values,
                                    WTD['pred_part22.5m'].values + WTD['WT_ESW_std'].values], axis=0)

    WTD['clearcut'] = 0.01 * np.nanmean([WTD['pred_clear4m'].values,
                                 WTD['pred_clear12m'].values], axis=0)
    WTD['clearcut_min'] = 0.01 * np.nanmin([WTD['pred_clear4m'].values - WTD['WT_ESW_std'].values,
                                    WTD['pred_clear12m'].values - WTD['WT_ESW_std'].values], axis=0)
    WTD['clearcut_max'] = 0.01 * np.nanmax([WTD['pred_clear4m'].values + WTD['WT_ESW_std'].values,
                                    WTD['pred_clear12m'].values + WTD['WT_ESW_std'].values], axis=0)

    # lopulliset kuvaan
    plt.figure()
    plt.fill_between(WTD.index, WTD['control_max'].values, WTD['control_min'].values,
                     facecolor='k', alpha=0.3)
    plt.plot(WTD.index, WTD['control'].values,'-k', linewidth=1.0)

    plt.fill_between(WTD.index, WTD['partial_max'].values, WTD['partial_min'].values,
                     facecolor='b', alpha=0.3)
    plt.plot(WTD.index, WTD['partial'].values,'-b', linewidth=1.0)

    plt.fill_between(WTD.index, WTD['clearcut_max'].values, WTD['clearcut_min'].values,
                     facecolor='r', alpha=0.3)
    plt.plot(WTD.index, WTD['clearcut'].values,'-r', linewidth=1.0)

    # ja tiedostoon
    save_df_to_csv(WTD[['control','control_max','control_min','partial','partial_max','partial_min','clearcut','clearcut_max','clearcut_min']],
                       'lettosuo_WTD_pred', readme=' - Check timezone!! \nSee Lettosuo_dataprocessing.process_WTD()')

def read_weir_data():

    fpaths = [r"H:\Lettosuo\WTD_paavolta\pato\logger_data.csv",
               r"H:\Lettosuo\WTD_paavolta\pato\tarkistusmittaukset.csv"]

    index=pd.date_range('01-01-2012','06-01-2019',freq='1H')
    data=pd.DataFrame(index=index, columns=[])

    for fp in fpaths:
        dat = pd.read_csv(fp, sep=';', header='infer', encoding = 'ISO-8859-1')
        dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
        dat.index = dat.index + datetime.timedelta(hours=0.5)
        dat.index = dat.index.map(lambda x: x.replace(minute=0, second=0))
        dat.ix[:,0]=dat.index
        dat = dat.drop_duplicates(subset=dat.columns[0])
        dat = dat.drop([dat.columns[0]], axis=1)
        if len(np.setdiff1d(dat.index, index)) > 0:
            print(fp, np.setdiff1d(dat.index, index))
            raise ValueError("Error")
        data=data.merge(dat, how='outer', left_index=True, right_index=True)

    return data

def process_runoff_data():

    weir_data = read_weir_data()

    plt.figure()
    ax=plt.subplot(4,1,1)
    plt.title('Loggeri- ja tarkistusmittausten erotus ja interpoloitu virhe (mm)')
    weir_data['error'] = weir_data['Water Height Point mm'] - weir_data['Tarkistusmittaus']
    weir_data['interpolated error'] = weir_data['error'].interpolate()
    weir_data['interpolated error'] = weir_data['interpolated error'].fillna(method='bfill')
    plt.plot(weir_data.index,weir_data['interpolated error'],'--r')
    plt.plot(weir_data.index,weir_data['Water Height Point mm'] - weir_data['Tarkistusmittaus'],'or')
    ax=plt.subplot(4,1,2,sharex=ax)
    plt.title('Loggerin vedenpinta, korjattu vedepinta ja tarkistusmittaukset (mm)')
    weir_data['Water Height Corrected mm'] = weir_data['Water Height Point mm'] -  weir_data['interpolated error']
    plt.plot(weir_data.index,weir_data['Water Height Point mm'])
    plt.plot(weir_data.index,weir_data['Water Height Corrected mm'],'--r')
    plt.plot(weir_data.index,weir_data['Tarkistusmittaus'],'or')
    ax=plt.subplot(4,1,3,sharex=ax)
    plt.title('Virtaama (l/s) - kaavalla Q=min(1381*(h[m])^2.5, 60)')
    weir_data['Discharge l/s uncorrected'] = 1381.0*(weir_data['Water Height Point mm']/1000.0)**2.5
    weir_data['Discharge l/s uncorrected'] = np.minimum(weir_data['Discharge l/s uncorrected'],60.0)
    weir_data['Discharge l/s'] = 1381.0*(weir_data['Water Height Corrected mm']/1000.0)**2.5
    weir_data['Discharge l/s'] = np.minimum(weir_data['Discharge l/s'],60.0)
    plt.plot(weir_data.index,weir_data['Discharge l/s uncorrected'])
    plt.plot(weir_data.index,weir_data['Discharge l/s'],'r')
    ax=plt.subplot(4,1,4,sharex=ax)
    plt.title('Valunta (mm/h) - valuma-alueen pinta-alana käytetty 17.1 ha')
    weir_data['Runoff mm/h'] = weir_data['Discharge l/s'] * 0.001 * 3600. / 171000. * 1000.
    weir_data['Runoff mm/h uncorrected'] = weir_data['Discharge l/s uncorrected'] * 0.001 * 3600. / 171000. * 1000.
    plt.plot(weir_data.index,weir_data['Runoff mm/h uncorrected'])
    plt.plot(weir_data.index,weir_data['Runoff mm/h'],'r')
    plt.tight_layout()

    plt.figure()
    ax=plt.subplot(2,1,1)
    plt.title('Erotus vedenpinta (mm)')
    plt.plot(weir_data.index,weir_data['Water Height Point mm'] - weir_data['Water Height Corrected mm'])
    print(np.nanmean(np.where(abs((weir_data['Water Height Point mm'] - weir_data['Water Height Corrected mm']).values) < 100,
                          abs((weir_data['Water Height Point mm'] - weir_data['Water Height Corrected mm']).values), np.nan)))
    print(np.nanmax(np.where(abs((weir_data['Water Height Point mm'] - weir_data['Water Height Corrected mm']).values) < 100,
                          abs((weir_data['Water Height Point mm'] - weir_data['Water Height Corrected mm']).values), np.nan)))

    plt.subplot(2,1,2,sharex=ax)
    plt.title('Erotus valunta (mm/h)')
    plt.plot(weir_data.index,weir_data['Runoff mm/h'] - weir_data['Runoff mm/h uncorrected'])
    print(np.nanmean(np.where(abs((weir_data['Runoff mm/h uncorrected'] - weir_data['Runoff mm/h']).values) < 0.5,
                          abs((weir_data['Runoff mm/h uncorrected'] - weir_data['Runoff mm/h']).values), np.nan)))
    print(np.nanmax(np.where(abs((weir_data['Runoff mm/h uncorrected'] - weir_data['Runoff mm/h']).values) < 0.5,
                          abs((weir_data['Runoff mm/h uncorrected'] - weir_data['Runoff mm/h']).values), np.nan)))

    # ja tiedostoon
    save_df_to_csv(weir_data[['Water Height Point mm','Tarkistusmittaus','Water Height Corrected mm','Discharge l/s','Runoff mm/h']],
                       'lettosuo_weir_data', readme=' - Check timezone!! \nSee Lettosuo_dataprocessing.process_runoff_data()')

def fit_pf_Laiho():
    from pyAPES_utilities.parameter_utilities import fit_pF
    # heads [kPa]
    head = [0.01, 0.3, 0.981, 4.905, 9.81, 33.0, 98.1]

    # volumetric water content measurements corresponding to heads for different peat types [%]
    watcont = [[94.69, 49.42, 29.61, 21.56, 20.05, 17.83, 16.54],
               [91.41, 66.26, 56.98, 45.58, 41.44, 39.32, 37.89],
               [89.12, -999, 72.83, 63.97, 54.40, 50.15, 48.80],
               [89.46, -999, 82.46, 76.79, 66.93, 63.61, 62.53],
               [92.22, -999, 87.06, 78.02, 74.76, 72.77, 71.70],
               [91.98, 66.70, 57.49, 39.95, 34.41, 29.83, 28.39],
               [88.75, 78.98, 78.77, 75.83, 72.37, 61.35, 45.66],
               [91.93, -999, 83.65, 78.39, 75.56, 74.08, 73.10],
               [93.45, 87.44, 87.33, 77.95, 76.46, 75.01, 73.01],
               [93.32, 87.15, 86.73, 82.90, 81.84, 80.51, 79.55],
               [93.05, 54.55, 42.88, 33.98, 29.80, 26.90, 25.87],
               [92.90, 78.15, 72.19, 54.47, 51.05, 49.66, 48.35],
               [90.11, -999, 80.86, 77.98, 69.66, 60.70, 53.40],
               [93.14, -999, 83.81, 78.69, 74.86, 73.26, 72.25],
               [93.17, -999, 89.76, 80.05, 76.67, 74.21, 72.84]]

    fit_pF(head, watcont[0:5], fig=True,percentage=True, kPa=True)
    fit_pF(head, watcont[5:10], fig=True,percentage=True, kPa=True)
    fit_pF(head, watcont[10:15], fig=True,percentage=True, kPa=True)

def EC_data():

    # PERIOD START / PERIOD END???
    # Energy fluxes on hourly time scale??

    starttime = '01-01-2010'
    endtime = '01-01-2019'

    fpaths = [r"H:\Lettosuo\Forcing_data\datat\lettosuo_28022019\pythonlukee\osittaishakkuu_energy.csv",  # period start/end??
              r"H:\Lettosuo\Forcing_data\datat\lettosuo_28022019\pythonlukee\osittaishakkuu_CO2.csv",
              r"H:\Lettosuo\Forcing_data\datat\lettosuo_28022019\pythonlukee\avohakkuu_energy.csv",  # period start/end??
              r"H:\Lettosuo\Forcing_data\datat\lettosuo_28022019\pythonlukee\avohakkuu_CO2.csv",
              r"H:\Lettosuo\Forcing_data\Annalea2\Letto1_EC.csv",
              r"H:\Lettosuo\Forcing_data\Annalea1\Letto1_meteo.csv"]  # period start

    fp_yearfirst = [r"H:\Forcing_data\datat\lettosuo_28.2.2019\osittaishakkuu_CO2.csv",
                    r"H:\Lettosuo\Forcing_data\Annalea1\Letto1_meteo.csv"]

    fp_hourly = [r"H:\Lettosuo\Forcing_data\datat\lettosuo_28022019\pythonlukee\osittaishakkuu_energy.csv",
                r"H:\Lettosuo\Forcing_data\datat\lettosuo_28022019\pythonlukee\avohakkuu_energy.csv"]

    fp_periodstart = [r"H:\Lettosuo\Forcing_data\Annalea1\Letto1_meteo.csv"]

    index = pd.date_range(starttime, endtime, freq='0.5H')
    data = pd.DataFrame(index=index, columns=[])

    for fp in fpaths:
        dat = pd.read_csv(fp, sep=',', header='infer', encoding = 'ISO-8859-1')
        if fp in fp_yearfirst:
            dat.index = pd.to_datetime(dat.ix[:,0], yearfirst=True)
        else:
            dat.index = pd.to_datetime(dat.ix[:,0], dayfirst=True)
        dat.index=dat.index.map(lambda x: x.replace(second=0))
        dat.ix[:,0]=dat.index
        if fp in fp_periodstart:  # period start to period end
            dat.index = dat.index + pd.Timedelta(hours=0.5)
        dat = dat.drop_duplicates(subset=dat.columns[0])
        dat = dat.drop([dat.columns[0]], axis=1)
        dat.columns = fp.split("\\")[-1].split(".")[0] + ': ' + dat.columns
        dat = dat[(dat.index >= starttime) & (dat.index < endtime)]
        if len(np.setdiff1d(dat.index, index)) > 0:
            print(fp, np.setdiff1d(dat.index, index))
            raise ValueError("Error")
        data=data.merge(dat, how='outer', left_index=True, right_index=True)
        if fp in fp_hourly:
            for col in dat.columns:
                data[col][1:-1:2] = data[col][0:-2:2].values

    variables = ['osittaishakkuu_energy: LE [W m-2]',
                 'osittaishakkuu_energy: SH [W m-2]',
                 'osittaishakkuu_CO2: NEE [mg CO2 m-2 s-1]',
                 'avohakkuu_energy: LE [W m-2]',
                 'avohakkuu_energy: SH [W m-2]',
                 'avohakkuu_CO2: NEE [mg CO2 m-2 s-1]',
                 ]
    flags = ['osittaishakkuu_energy: Flag_LE',
             'osittaishakkuu_energy: Flag_SH',
             'osittaishakkuu_CO2: gapped',
             'avohakkuu_energy: Flag_LE',
             'avohakkuu_energy: Flag_SH',
             'avohakkuu_CO2: Flag_NEE'
             ]

    for i in range(len(variables)):
        data[variables[i] + ', not gapfilled'] = np.where(
                data[flags[i]] == 0, data[variables[i]], np.nan)

    data['Letto1_meteo: avg(GlobRefl (W/m2))'] = np.where(data['Letto1_meteo: avg(GlobRefl (W/m2))'] < data['Letto1_meteo: avg(Glob (W/m2))'],
                            data['Letto1_meteo: avg(GlobRefl (W/m2))'], np.nan)
    data['Letto1_meteo: avg(GlobRefl (W/m2))'] = np.where(data['Letto1_meteo: avg(GlobRefl (W/m2))'] > 0.0,
                            data['Letto1_meteo: avg(GlobRefl (W/m2))'], np.nan)
    data['control_NSWRAD'] = data['Letto1_meteo: avg(Glob (W/m2))'] - data['Letto1_meteo: avg(GlobRefl (W/m2))']

    lettosuo_EC = data[[
           'osittaishakkuu_energy: LE [W m-2], not gapfilled',
           'osittaishakkuu_energy: SH [W m-2], not gapfilled',
           'osittaishakkuu_energy: SHF [W m-2]',
           'osittaishakkuu_energy: NRAD [W m-2]',
           'osittaishakkuu_CO2: NEE [mg CO2 m-2 s-1], not gapfilled',
           'osittaishakkuu_CO2: Reco [mg CO2 m-2 s-1]',
           'osittaishakkuu_CO2: GPP [mg CO2 m-2 s-1]',
           'avohakkuu_energy: LE [W m-2], not gapfilled',
           'avohakkuu_energy: SH [W m-2], not gapfilled',
           'avohakkuu_energy: SHF [W m-2]',
           'avohakkuu_energy: NRAD [W m-2]',
           'avohakkuu_CO2: NEE [mg CO2 m-2 s-1], not gapfilled',
           'avohakkuu_CO2: Reco [mg CO2 m-2 s-1]',
           'avohakkuu_CO2: GPP [mg CO2 m-2 s-1]',
           'Letto1_EC: LE (Wm-2)',
           'Letto1_EC: SH (W m-2)',
           'Letto1_meteo: avg(NetRad (W/m2))',
           'control_NSWRAD',
           'Letto1_EC: NEE_South',
           'Letto1_EC: NEE_North',
           'Letto1_EC: GPP_S',
           'Letto1_EC: GPP_N',
           'Letto1_EC: Reco_S',
           'Letto1_EC: Reco_N'
            ]].copy()

    lettosuo_EC['Letto1_EC: NEE'] = np.where(lettosuo_EC['Letto1_EC: NEE_South'].notnull(),
               lettosuo_EC['Letto1_EC: NEE_South'], lettosuo_EC['Letto1_EC: NEE_North'])

    lettosuo_EC = lettosuo_EC.rename(columns={
           'osittaishakkuu_energy: LE [W m-2], not gapfilled': 'partial_LE',
           'osittaishakkuu_energy: SH [W m-2], not gapfilled': 'partial_SH',
           'osittaishakkuu_energy: NRAD [W m-2]': 'partial_GHF',  # huom nimet!
           'osittaishakkuu_energy: SHF [W m-2]': 'partial_NRAD',  # huom nimet!
           'osittaishakkuu_CO2: NEE [mg CO2 m-2 s-1], not gapfilled': 'partial_NEE',
           'osittaishakkuu_CO2: GPP [mg CO2 m-2 s-1]': 'partial_GPP',
           'osittaishakkuu_CO2: Reco [mg CO2 m-2 s-1]': 'partial_Reco',
           'avohakkuu_energy: LE [W m-2], not gapfilled': 'clearcut_LE',
           'avohakkuu_energy: SH [W m-2], not gapfilled': 'clearcut_SH',
           'avohakkuu_energy: SHF [W m-2]': 'clearcut_GHF',
           'avohakkuu_energy: NRAD [W m-2]': 'clearcut_NRAD',
           'avohakkuu_CO2: NEE [mg CO2 m-2 s-1], not gapfilled': 'clearcut_NEE',
           'avohakkuu_CO2: GPP [mg CO2 m-2 s-1]': 'clearcut_GPP',
           'avohakkuu_CO2: Reco [mg CO2 m-2 s-1]': 'clearcut_Reco',
           'Letto1_EC: LE (Wm-2)': 'control_LE',
           'Letto1_EC: SH (W m-2)': 'control_SH',
           'Letto1_meteo: avg(NetRad (W/m2))': 'control_NRAD',
           'Letto1_EC: NEE_North': 'control_NEE',
           'Letto1_EC: GPP_N': 'control_GPP',
           'Letto1_EC: Reco_N': 'control_Reco',
           'Letto1_EC: NEE_South': 'control-S_NEE',
           'Letto1_EC: GPP_S': 'control-S_GPP',
           'Letto1_EC: Reco_S': 'control-S_Reco'
           })

    lettosuo_EC['partial_NRAD'] = lettosuo_EC['control_NRAD'].copy()
    lettosuo_EC['partial_NSWRAD'] = lettosuo_EC['control_NSWRAD'].copy()

    lettosuo_EC['partial_NRAD'][lettosuo_EC.index < '01-01-2016'] = np.nan
    lettosuo_EC['partial_NSWRAD'][lettosuo_EC.index < '01-01-2016'] = np.nan

    lettosuo_EC['control_NRAD'][lettosuo_EC.index >= '01-01-2016'] = np.nan
    lettosuo_EC['control_NSWRAD'][lettosuo_EC.index >= '01-01-2016'] = np.nan

    for column in lettosuo_EC.columns:
        if column.split('_')[-1] == 'LE':
            lettosuo_EC[column][lettosuo_EC[column] < -100] = lettosuo_EC[column][lettosuo_EC[column] < -100] * np.nan
            lettosuo_EC[column][lettosuo_EC[column] > 1000] = lettosuo_EC[column][lettosuo_EC[column] > 1000] * np.nan

    treatments=['control','partial','clearcut']
    variables=['_NRAD', '_LE', '_SH', '_NEE', '_GPP', '_Reco']
    plt.figure()
    for i in range(len(variables)):
        if i == 0:
            ax1 = plt.subplot(len(variables),1,i+1)
            lettosuo_EC[[treatment + variables[i] for treatment in treatments]].plot(ax=ax1)
        else:
            ax = plt.subplot(len(variables),1,i+1,sharex=ax1)
            lettosuo_EC[[treatment + variables[i] for treatment in treatments]].plot(ax=ax)

    readme = "\nTreatment energy balance from hourly data\nPeriod end/start unclear"
    save_df_to_csv(lettosuo_EC, "Lettosuo_EC_2010_2018", readme=readme)

    return data

def Tsoil_wsoil_data(starttime='01-01-2010',endtime='01-01-2019'):

    # filepaths
    forc_fp = ["H:/Lettosuo/Forcing_data/Annalea1/Letto1_metsanpohja.csv",
               "H:/Lettosuo/Forcing_data/Annalea1/Letto2_metsanpohja.csv",
               "H:/Lettosuo/Forcing_data/MikaK/metsanpohja_concat.csv",
               "H:/Lettosuo/Forcing_data/datat/kontrolli_meteo/concat.csv",
               "H:/Lettosuo/Forcing_data/datat/Avohakkuu_meteo.csv"]

    index=pd.date_range(starttime, endtime, freq='0.5H')
    lettosuo_data=pd.DataFrame(index=index, columns=[])

    for fp in forc_fp:
        dat = pd.read_csv(fp, sep=',', header='infer', encoding = 'ISO-8859-1')
        dat.index = pd.to_datetime(dat.ix[:,0], yearfirst=True)
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

    lettosuo_data = lettosuo_data[(lettosuo_data.index >= starttime) &
                                  (lettosuo_data.index < endtime)]

    columns1=[
           'Letto1_metsanpohja: avg(Moist1 7cm (m3/m3))',
           'Letto1_metsanpohja: avg(Moist2 20cm (m3/m3))',
           'Letto1_metsanpohja: avg(T1 5cm (C))',
           'Letto1_metsanpohja: avg(T2 15cm (C))',
           'Letto1_metsanpohja: avg(T3 30cm (C))',
           'Letto1_metsanpohja: avg(T4 40cm (C))',
           'Letto2_metsanpohja: avg(SoilMoist-10cm (m3/m3))',
           'Letto2_metsanpohja: avg(SoilMoist-20cm (m3/m3))',
           'Letto2_metsanpohja: avg(T-5cm (C))',
           'Letto2_metsanpohja: avg(T-10cm (C))',
           'Letto2_metsanpohja: avg(T-20cm (C))',
           'Letto2_metsanpohja: avg(T-30cm (C))'
           ]
    columns2=[
           'metsanpohja_concat: avg(Moist1 7cm (?))',
           'metsanpohja_concat: avg(Moist2 20cm (?))',
           'metsanpohja_concat: avg(T1 5cm (C))',
           'metsanpohja_concat: avg(T2 15cm (C))',
           'metsanpohja_concat: avg(T3 30cm (C))',
           'metsanpohja_concat: avg(T4 40cm (C))',
           'concat: avg(SoilMoist-10cm (m3/m3))',
           'concat: avg(SoilMoist-20cm (m3/m3))',
           'concat: avg(T-5cm (C))',
           'concat: avg(T-10cm (C))',
           'concat: avg(T-20cm (C))',
           'concat: avg(T-30cm (C))',
           ]

    for i in range(len(columns1)):
        lettosuo_data[columns1[i]] = np.where(
            lettosuo_data.index < '01-01-2018',
            lettosuo_data[columns1[i]],
            lettosuo_data[columns2[i]])

    T_W_data = lettosuo_data[[
           'Letto1_metsanpohja: avg(Moist1 7cm (m3/m3))',
           'Letto1_metsanpohja: avg(Moist2 20cm (m3/m3))',
           'Letto1_metsanpohja: avg(T1 5cm (C))',
           'Letto1_metsanpohja: avg(T2 15cm (C))',
           'Letto1_metsanpohja: avg(T3 30cm (C))',
           'Letto1_metsanpohja: avg(T4 40cm (C))',
           'Letto2_metsanpohja: avg(SoilMoist-10cm (m3/m3))',
           'Letto2_metsanpohja: avg(SoilMoist-20cm (m3/m3))',
           'Letto2_metsanpohja: avg(T-5cm (C))',
           'Letto2_metsanpohja: avg(T-10cm (C))',
           'Letto2_metsanpohja: avg(T-20cm (C))',
           'Letto2_metsanpohja: avg(T-30cm (C))',
           'Avohakkuu_meteo: STCH4avg',
           'Avohakkuu_meteo: STCH5avg',
           'Avohakkuu_meteo: STCH6avg',
           'Avohakkuu_meteo: STCH7avg',
           'Avohakkuu_meteo: SMCHA',
           'Avohakkuu_meteo: SMCHB'
           ]].copy()

    T_W_data = T_W_data.rename(columns={
           'Letto1_metsanpohja: avg(Moist1 7cm (m3/m3))':'Letto1_SMA',
           'Letto1_metsanpohja: avg(Moist2 20cm (m3/m3))':'Letto1_SMB',
           'Letto1_metsanpohja: avg(T1 5cm (C))':'Letto1_T5',
           'Letto1_metsanpohja: avg(T2 15cm (C))':'Letto1_T15',
           'Letto1_metsanpohja: avg(T3 30cm (C))':'Letto1_T30',
           'Letto1_metsanpohja: avg(T4 40cm (C))':'Letto1_T40',
           'Letto2_metsanpohja: avg(SoilMoist-10cm (m3/m3))':'Letto2_SMA',
           'Letto2_metsanpohja: avg(SoilMoist-20cm (m3/m3))':'Letto2_SMB',
           'Letto2_metsanpohja: avg(T-5cm (C))':'Letto2_T5',
           'Letto2_metsanpohja: avg(T-10cm (C))':'Letto2_T10',
           'Letto2_metsanpohja: avg(T-20cm (C))':'Letto2_T20',
           'Letto2_metsanpohja: avg(T-30cm (C))':'Letto2_T30',
           'Avohakkuu_meteo: SMCHA':'Clear_SMA',
           'Avohakkuu_meteo: SMCHB':'Clear_SMB',
           'Avohakkuu_meteo: STCH4avg':'Clear_T5',
           'Avohakkuu_meteo: STCH5avg':'Clear_T10',
           'Avohakkuu_meteo: STCH6avg':'Clear_T20',
           'Avohakkuu_meteo: STCH7avg':'Clear_T30'
           })

    sites=['Letto1','Letto2','Clear']
    variables=['_SMA', '_SMB', '_T5', '_T30']
    plt.figure()
    for i in range(len(variables)):
        if i == 0:
            ax1 = plt.subplot(len(variables),1,i+1)
            T_W_data[[site + variables[i] for site in sites]].plot(ax=ax1)
        else:
            ax = plt.subplot(len(variables),1,i+1,sharex=ax1)
            T_W_data[[site + variables[i] for site in sites]].plot(ax=ax)

    readme = "\nSoil moisture data not ok?"
    save_df_to_csv(T_W_data, "Lettosuo_Tsoil_2010_2018", readme=readme)

    return lettosuo_data

def fit_pf_paivanen():
    from pyAPES_utilities.parameter_utilities import fit_pF

    # heads [kPa]
    head = [0.0001, 1, 3.2, 10, 20, 60, 100, 200, 500, 1000, 1500]
    watcont_sedge = [[94.3, 68, 47.9, 35.5, 22, 18.4, 16.7, 12.6, 7.9, 6.3, -999],
               [91.7, 82.9, 61.8, 35.9, 31.7, 25.2, 23.4, 19.2, 17.3, 14.4, -999],
               [90.6, 86.2, 56.4, 36.4, 33.4, 29.8, 26.8, 23.5, 20.2, 16.4, -999],
               [89.7, 85, 74.5, 53.4, 36, 29, 24.7, 22.1, 17.6, 14.7, -999],
               [87.3, 85.4, 77.8, 64, 41.4, 28.8, 23.4, 22.7, 21.9, 17, -999],
               [89.3, 86.5, 80.7, 52.5, 45.6, 35.4, 32, 25.1, 20.6, 18.4, -999],
               [91, 89.9, 84.7, 60, -999, 33.8, 27.2, 29.3, -999, 17.2, 12.9],
               [89.3, 87.2, 79.2, 59.8, 53.8, 46.3, 41.5, 36.1, 32, 28.6, -999],
               [87.2, 77.7, 76.2, 56.4, -999, 33.3, 31.6, 28.8, -999, 17.3, 15.7],
               [84.2, 83.7, 76, 56.6, 44.2, 41.2, 37.4, -999, 36.3, 32.6, -999],
               [83.9, 80.7, 78.8, 67.1, 54.7, 41, 37.4, 35.2, 29.3, 27, -999],
               [87.2, 84.3, 81.6, 71.7, 55.1, 43.4, 36.8, 33.4, 29, 27.3, -999],
               [82.4, 81.8, 70.5, 57.5, 50.1, 48.2, 44.2, 42.5, 41.9, 32.4, -999],
               [81.8, 77.7, 75.1, 64.7, 56.8, 43.4, 40, 32.6, 27.1, 27, -999]]

    fit_pF(head, watcont_sedge, fig=True,percentage=True, kPa=True)
    fit_pF(head, [np.average(watcont_sedge[1:], axis=0, weights=~np.isin(watcont_sedge[1:], -999)*1)], fig=True,percentage=True, kPa=True)