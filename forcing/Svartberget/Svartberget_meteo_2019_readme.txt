Readme for Svarberget_meteo_2019.csv

Kersti Haahti, Luke 2020-04-02

yyyy, mo, dd, hh, mm: datetime [UTC + 2.0]
Tsoil: Soil temperature at 5cm [degC]
  flag 0 (83.34%): TS_avg_1_1
  flag 1 (16.66%): linearly interpolated
  flag 0 (0.01%): filled with nearest
Wliq: Soil volumetric moisture content at 5cm [%]
  flag 0 (83.26%): SWC_avg_1_1
  flag 1 (16.73%): linearly interpolated
  flag 0 (0.01%): filled with nearest
Tair: Air temperature [degC]
  flag 0 (83.37%): Ta_1_1_1
  flag 1 (16.63%): linearly interpolated
  flag 0 (0.01%): filled with nearest
RH: Relative humidity [%]
  flag 0 (83.37%): RH_1_1_1
  flag 1 (16.63%): linearly interpolated
  flag 0 (0.01%): filled with nearest
P: Ambient pressure [hPa]
  flag 0 (83.38%): Pa_1_1_1
  flag 1 (16.61%): linearly interpolated
  flag 0 (0.01%): filled with nearest
Rg: Global radiation i.e. incoming shortwave radiation [W/m2]
  flag 0 (83.38%): Swin_1_1_1
  flag 1 (0.00%): Swin_1_2_1
  flag 2 (16.61%): linearly interpolated
  flag 1 (0.01%): filled with nearest
LWin: Downwelling long wave radiation [W/m2]
  flag 0 (83.38%): Lwin_1_2_1
  flag 1 (16.61%): linearly interpolated
  flag 0 (0.01%): filled with nearest
LWout: Upwelling long wave radiation [W/m2]
  flag 0 (83.38%): Lwout_1_2_1
  flag 1 (16.61%): linearly interpolated
  flag 0 (0.01%): filled with nearest
Prec: Precipitation  [mm/30min]
  flag 0 (83.03%): P_1_1_1
  flag 1 (16.97%): filled with Prec = 0.0
LWnet: Net longwave radiation [W m-2]
  flag 0 (83.38%): Lwnet_1_2_1
  flag 1 (16.62%): filled with LWnet = nan
SWnet: Net shortwave radiation [W m-2]
  flag 0 (83.38%): Swnet_1_2_1
  flag 1 (16.62%): filled with SWnet = nan
