Readme for Svarberget_meteo_2014_2016.csv

Kersti Haahti, Luke 2019-01-07

yyyy, mo, dd, hh, mm: datetime [UTC + 2.0]
Tsoil: Soil temperature at 5cm [degC]
  flag 0 (98.68%): SE-Svb_meteo: TS_avg_1_1
  flag 1 (1.32%): linearly interpolated
  flag 0 (0.00%): filled with nearest
Wliq: Soil volumetric moisture content at 5cm [%]
  flag 0 (96.91%): SE-Svb_eco: SWC_2_2_1_screened
  flag 1 (3.09%): linearly interpolated
  flag 0 (0.00%): filled with nearest
Tair: Air temperature [degC]
  flag 0 (95.42%): SE-Svb_meteo: Ta_1_1_1_screened
  flag 1 (4.58%): SE-Deg_meteo: Ta_1_1_1
  flag 2 (0.00%): linearly interpolated
  flag 1 (0.00%): filled with nearest
RH: Relative humidity [%]
  flag 0 (95.42%): SE-Svb_meteo: RH_1_1_1_screened
  flag 1 (4.58%): SE-Deg_meteo: RH_1_1_1
  flag 2 (0.00%): linearly interpolated
  flag 1 (0.00%): filled with nearest
P: Ambient pressure [hPa]
  flag 0 (94.70%): SE-Svb_meteo: Pa_1_1_1_screened
  flag 1 (5.29%): SE-Deg_meteo: Pa_1_1_1
  flag 2 (0.00%): linearly interpolated
  flag 1 (0.00%): filled with nearest
Rg: Global radiation i.e. incoming shortwave radiation [W/m2]
  flag 0 (89.93%): SE-Svb_meteo: Swin_1_1_1_screened
  flag 1 (5.49%): SE-Svb_meteo: Swin_1_2_1_screened
  flag 2 (4.59%): SE-Deg_meteo: Swin_1_1_1
  flag 3 (0.00%): SE-Deg_meteo: Swin_1_2_1
  flag 4 (0.00%): linearly interpolated
  flag 3 (0.00%): filled with nearest
LWin: Downwelling long wave radiation [W/m2]
  flag 0 (95.42%): SE-Svb_meteo: Lwin_1_2_1_screened
  flag 1 (4.58%): SE-Deg_meteo: Lwin_1_2_1
  flag 2 (0.00%): linearly interpolated
  flag 1 (0.00%): filled with nearest
LWout: Upwelling long wave radiation [W/m2]
  flag 0 (95.42%): SE-Svb_meteo: Lwout_1_2_1_screened
  flag 1 (4.58%): SE-Deg_meteo: Lwout_1_2_1
  flag 2 (0.00%): linearly interpolated
  flag 1 (0.00%): filled with nearest
U: Wind speed [m/s]
  flag 0 (89.73%): SE-Svb_fluxes: WS_1_1_1
  flag 1 (8.14%): SE-Deg_fluxes: WS_1_1_1
  flag 2 (2.13%): linearly interpolated
  flag 1 (0.00%): filled with nearest
Ustar: Friction velocity [m/s]
  flag 0 (89.52%): SE-Svb_fluxes: Ustar_1_1_1
  flag 1 (10.48%): Ustar = 0.17 * U
  flag 2 (0.00%): linearly interpolated
Prec: Precipitation  [mm/30min]
  flag 0 (78.10%): SE-Svb_meteo: P_1_1_1
  flag 1 (13.41%): SE-Deg_meteo: P_corrected to match daily Prec
  flag 2 (8.49%): Prec (mm)
  flag 3 (0.00%): filled with Prec = 0.0
