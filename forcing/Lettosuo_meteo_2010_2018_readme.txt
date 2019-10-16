Readme for Lettosuo_meteo_2010_2018.csv

Kersti Haahti, Luke 2019-10-16

yyyy, mo, dd, hh, mm: datetime [UTC + 2.0]
Prec_ref: Jokioinen gapfilled precipitaion [mm/30min]
  flag 0 (63.52%): jokioinen_prec1: Prec [mm h-1]
  flag 1 (7.10%): jokioinen_meteo: Precipitation amount
  flag 2 (28.85%): jokioinen_prec2: Prec [mm h-1]
  flag 3 (0.08%): somero_meteo: Precipitation amount
  flag 4 (0.45%): filled with Prec_ref = 0.0
Prec_old: Lettosuo gapfilled precipitaion [mm/30min]
  flag 0 (34.74%): Letto1_metsanpohja: avg(Rain (mm)) !sections removed!
  flag 1 (32.24%): somero_meteo: Precipitation amount
  flag 2 (33.01%): jokioinen_prec: Prec
  flag 3 (0.00%): filled with Prec_old = 0.0
Prec: Lettosuo gapfilled precipitaion [mm/30min]
  flag 0 (34.98%): Letto1_metsanpohja: avg(Rain (mm)) !sections removed!
  flag 1 (32.24%): somero_meteo: Precipitation amount
  flag 2 (32.78%): jokioinen_prec: Prec
  flag 3 (0.00%): filled with Prec = 0.0
Prec_ref2: Somero gapfilled precipitaion [mm/30min]
  flag 0 (48.47%): somero_meteo: Precipitation amount
  flag 1 (51.53%): jokioinen_prec: Prec
  flag 2 (0.00%): filled with Prec_ref2 = 0.0
Tair: Air temperature [degC]
  flag 0 (91.44%): Letto1_meteo: avg(Temp (C))
  flag 1 (5.04%): somero_meteo: Air temperature
  flag 2 (0.03%): jokioinen_meteo: Air temperature
  flag 3 (3.44%): Letto1_meteo_gapfilled: PaikAirT T
  flag 4 (0.00%): linearly interpolated
  flag 3 (0.04%): filled with nearest
RH: Relative humidity [%]
  flag 0 (91.44%): Letto1_meteo: avg(RH (%))
  flag 1 (5.03%): somero_meteo: Relative humidity
  flag 2 (0.03%): jokioinen_meteo: Relative humidity
  flag 3 (3.45%): Letto1_meteo_gapfilled: Paik RH
  flag 4 (0.00%): linearly interpolated
  flag 3 (0.04%): filled with nearest
Rg: Global radiation [W/m2]
  flag 0 (91.44%): Letto1_meteo: avg(Glob (W/m2))
  flag 1 (4.56%): Letto1_meteo_gapfilled: PaikGlob2
  flag 2 (3.84%): jokioinen_rad: Global radiation
  flag 3 (0.11%): linearly interpolated
  flag 2 (0.04%): filled with nearest
U: Wind speed [m/s]
  flag 0 (87.27%): Letto1_EC: wind speed (m/s) !u > 10 removed!
  flag 1 (6.97%): somero_meteo: Wind speed
  flag 2 (3.74%): Derived from salo_kiikala_meteo: wind speed (0.57*U_ref + 0.55)
  flag 3 (1.97%): linearly interpolated
  flag 2 (0.05%): filled with nearest
Ustar: Friction velocity [m/s]
  flag 0 (85.46%): Letto1_EC: friction velocity (m/s)
  flag 1 (14.54%): Ustar = 0.2 * U
  flag 2 (0.00%): linearly interpolated
P: Ambient pressure [hPa]
  flag 0 (42.14%): Letto1_meteo: avg(Press (hPa))
  flag 1 (53.71%): Derived from salo_kiikala_meteo: Pressure (msl) (P_ref - 16.9)
  flag 2 (0.56%): linearly interpolated
  flag 1 (3.58%): filled with nearest
Snow_depth1: Jokioinen snow depth [cm]
  flag 0 (43.35%): jokioinen_meteo: Snow depth
  flag 1 (56.65%): filled with Snow_depth1 = nan
Snow_depth2: Somero snow depth [cm]
  flag 0 (54.23%): somero_meteo: Snow depth
  flag 1 (45.77%): filled with Snow_depth2 = nan
Snow_depth3: Salo-kiikala snow depth [cm]
  flag 0 (52.52%): salo_kiikala_meteo: Snow depth
  flag 1 (47.48%): filled with Snow_depth3 = nan
