Readme for Svartberget_meteo_2019.csv

Kersti Haahti, Luke 2020-05-15

yyyy, mo, dd, hh, mm: datetime [UTC + 1.0]
Tsoil: Soil temperature at 5cm [degC]
  flag 0 (83.34%): TS_avg_1_1
  flag 1 (16.66%): linearly interpolated
  flag 0 (0.01%): filled with nearest
Wliq: Soil volumetric moisture content at 5cm [%]
  flag 0 (83.26%): SWC_avg_1_1
  flag 1 (16.73%): linearly interpolated
  flag 0 (0.01%): filled with nearest
Tair: Air temperature [degC]
  flag 0 (98.70%): Ta_1_1_1
  flag 1 (1.29%): linearly interpolated
  flag 0 (0.01%): filled with nearest
RH: Relative humidity [%]
  flag 0 (98.70%): RH_1_1_1
  flag 1 (1.29%): linearly interpolated
  flag 0 (0.01%): filled with nearest
P: Ambient pressure [hPa]
  flag 0 (97.27%): Pa_1_1_1
  flag 1 (2.72%): linearly interpolated
  flag 0 (0.01%): filled with nearest
Rg: Global radiation i.e. incoming shortwave radiation [W/m2]
  flag 0 (98.72%): Swin_1_1_1
  flag 1 (0.01%): Swin_1_2_1
  flag 2 (1.27%): linearly interpolated
  flag 1 (0.01%): filled with nearest
LWin: Downwelling long wave radiation [W/m2]
  flag 0 (98.72%): Lwin_1_2_1
  flag 1 (1.28%): linearly interpolated
  flag 0 (0.01%): filled with nearest
U: Wind speed [m/s]
  flag 0 (99.99%): wndspd_f_ms-1
  flag 1 (0.00%): linearly interpolated
  flag 0 (0.01%): filled with nearest
Ustar: Friction velocity [m/s]
  flag 0 (87.00%): Ustar_ms-1
  flag 1 (13.00%): Ustar = 0.19 * U
  flag 2 (0.00%): linearly interpolated
Prec: Precipitation [mm/30min]
  flag 0 (99.08%): P_1_1_1
  flag 1 (0.92%): filled with Prec = 0.0
CO2: Mixing ration of CO2 [ppm]
  flag 0 (84.72%): CO2_MixingRatio_umolmol-1
  flag 1 (15.27%): linearly interpolated
  flag 0 (0.01%): filled with nearest
