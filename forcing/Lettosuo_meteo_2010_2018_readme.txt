Readme for Lettosuo_meteo_2010_2018.csv

Kersti Haahti, Luke 2018-06-01

yyyy, mo, dd, hh, mm: datetime [UTC + 2.0]
Prec_ref: Jokioinen gapfilled precipitaion [mm/30min]
  flag 0 (69.92%): jokioinen_prec1: Prec [mm h-1]
  flag 1 (29.56%): jokioinen_prec2: Prec [mm h-1]
  flag 2 (0.07%): somero_meteo: Precipitation amount
  flag 3 (0.45%): filled with Prec_ref = 0.0
Prec: Lettosuo gapfilled precipitaion [mm/30min]
  flag 0 (35.36%): Letto1_metsanpohja: avg(Rain (mm)) !sections removed!
  flag 1 (44.38%): jokioinen_prec1: Prec [mm h-1]
  flag 2 (19.89%): jokioinen_prec2: Prec [mm h-1]
  flag 3 (0.06%): somero_meteo: Precipitation amount
  flag 4 (0.31%): filled with Prec = 0.0
Tair: Air temperature [degC]
  flag 0 (93.31%): Letto1_meteo: avg(Temp (C))
  flag 1 (3.77%): somero_meteo: Air temperature
  flag 2 (0.04%): jokioinen_meteo: Air temperature
  flag 3 (2.87%): Letto1_meteo_gapfilled: PaikAirT T
  flag 4 (0.00%): linearly interpolated
RH: Relative humidity [%]
  flag 0 (93.31%): Letto1_meteo: avg(RH (%))
  flag 1 (3.76%): somero_meteo: Relative humidity
  flag 2 (0.04%): jokioinen_meteo: Relative humidity
  flag 3 (2.89%): Letto1_meteo_gapfilled: Paik RH
  flag 4 (0.00%): linearly interpolated
Rg: Global radiation [W/m2]
  flag 0 (93.31%): Letto1_meteo: avg(Glob (W/m2))
  flag 1 (4.31%): Letto1_meteo_gapfilled: PaikGlob2
  flag 2 (2.37%): jokioinen_rad: Global radiation
  flag 3 (0.00%): linearly interpolated
U: Wind speed [m/s]
  flag 0 (67.00%): Letto1_EC: wind speed (m/s) !u > 10 removed!
  flag 1 (27.54%): somero_meteo: Wind speed
  flag 2 (4.75%): Derived from salo_kiikala_meteo: wind speed (0.57*U_ref + 0.55)
  flag 3 (0.71%): linearly interpolated
Ustar: Friction velocity [m/s]
  flag 0 (65.44%): Letto1_EC: friction velocity (m/s)
  flag 1 (34.56%): Ustar = 0.2 * U
  flag 2 (0.00%): linearly interpolated
P: Ambient pressure [hPa]
  flag 0 (38.79%): Letto1_meteo: avg(Press (hPa))
  flag 1 (60.57%): Derived from salo_kiikala_meteo: Pressure (msl) (P_ref - 16.9)
  flag 2 (0.64%): linearly interpolated
  flag 1 (0.00%): filled with nearest
