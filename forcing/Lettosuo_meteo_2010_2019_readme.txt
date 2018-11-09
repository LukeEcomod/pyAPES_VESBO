Readme for Lettosuo_meteo_2010_2019.csv

Kersti Haahti, Luke 2018-11-09

yyyy, mo, dd, hh, mm: datetime [UTC + 2.0]
Prec_ref: Jokioinen gapfilled precipitaion [mm/30min]
  flag 0 (65.88%): jokioinen_prec1: Prec [mm h-1]
  flag 1 (4.58%): jokioinen_meteo: Precipitation amount
  flag 2 (26.28%): jokioinen_prec2: Prec [mm h-1]
  flag 3 (1.01%): somero_meteo: Precipitation amount
  flag 4 (2.25%): filled with Prec_ref = 0.0
Prec: Lettosuo gapfilled precipitaion [mm/30min]
  flag 0 (34.56%): Letto1_metsanpohja: avg(Rain (mm)) !sections removed!
  flag 1 (32.13%): somero_meteo: Precipitation amount
  flag 2 (13.50%): jokioinen_prec1: Prec [mm h-1]
  flag 3 (0.02%): jokioinen_meteo: Precipitation amount
  flag 4 (17.67%): jokioinen_prec2: Prec [mm h-1]
  flag 5 (2.13%): filled with Prec = 0.0
Tair: Air temperature [degC]
  flag 0 (90.69%): Letto1_meteo: avg(Temp (C))
  flag 1 (4.87%): somero_meteo: Air temperature
  flag 2 (0.03%): jokioinen_meteo: Air temperature
  flag 3 (2.55%): Letto1_meteo_gapfilled: PaikAirT T
  flag 4 (1.86%): linearly interpolated
RH: Relative humidity [%]
  flag 0 (90.69%): Letto1_meteo: avg(RH (%))
  flag 1 (4.86%): somero_meteo: Relative humidity
  flag 2 (0.03%): jokioinen_meteo: Relative humidity
  flag 3 (2.57%): Letto1_meteo_gapfilled: Paik RH
  flag 4 (1.86%): linearly interpolated
Rg: Global radiation [W/m2]
  flag 0 (90.69%): Letto1_meteo: avg(Glob (W/m2))
  flag 1 (3.84%): Letto1_meteo_gapfilled: PaikGlob2
  flag 2 (3.62%): jokioinen_rad: Global radiation
  flag 3 (1.86%): linearly interpolated
U: Wind speed [m/s]
  flag 0 (86.75%): Letto1_EC: wind speed (m/s) !u > 10 removed!
  flag 1 (6.90%): somero_meteo: Wind speed
  flag 2 (3.88%): Derived from salo_kiikala_meteo: wind speed (0.57*U_ref + 0.55)
  flag 3 (2.47%): linearly interpolated
Ustar: Friction velocity [m/s]
  flag 0 (84.91%): Letto1_EC: friction velocity (m/s)
  flag 1 (15.09%): Ustar = 0.2 * U
  flag 2 (0.00%): linearly interpolated
P: Ambient pressure [hPa]
  flag 0 (42.22%): Letto1_meteo: avg(Press (hPa))
  flag 1 (55.35%): Derived from salo_kiikala_meteo: Pressure (msl) (P_ref - 16.9)
  flag 2 (2.43%): linearly interpolated
  flag 1 (0.00%): filled with nearest
