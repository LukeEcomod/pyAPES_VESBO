Readme for Lettosuo_meteo_2010_2018.csv

Kersti Haahti, Luke 2019-04-05

yyyy, mo, dd, hh, mm: datetime [UTC + 2.0]
Prec_ref2: Somero gapfilled precipitaion [mm/30min]
  flag 0 (50.27%): somero_meteo: Precipitation amount
  flag 1 (23.04%): jokioinen_prec1: Prec [mm h-1]
  flag 2 (0.04%): jokioinen_meteo: Precipitation amount
  flag 3 (26.25%): jokioinen_prec2: Prec [mm h-1]
  flag 4 (0.40%): filled with Prec_ref2 = 0.0
Prec_ref: Jokioinen gapfilled precipitaion [mm/30min]
  flag 0 (65.88%): jokioinen_prec1: Prec [mm h-1]
  flag 1 (7.36%): jokioinen_meteo: Precipitation amount
  flag 2 (26.28%): jokioinen_prec2: Prec [mm h-1]
  flag 3 (0.08%): somero_meteo: Precipitation amount
  flag 4 (0.40%): filled with Prec_ref = 0.0
Prec: Lettosuo gapfilled precipitaion [mm/30min]
  flag 0 (35.07%): Letto1_metsanpohja: avg(Rain (mm)) !sections removed!
  flag 1 (33.46%): somero_meteo: Precipitation amount
  flag 2 (13.50%): jokioinen_prec1: Prec [mm h-1]
  flag 3 (0.03%): jokioinen_meteo: Precipitation amount
  flag 4 (17.67%): jokioinen_prec2: Prec [mm h-1]
  flag 5 (0.27%): filled with Prec = 0.0
Tair: Air temperature [degC]
  flag 0 (92.18%): Letto1_meteo: avg(Temp (C))
  flag 1 (5.23%): somero_meteo: Air temperature
  flag 2 (0.03%): jokioinen_meteo: Air temperature
  flag 3 (2.55%): Letto1_meteo_gapfilled: PaikAirT T
  flag 4 (0.00%): linearly interpolated
RH: Relative humidity [%]
  flag 0 (92.18%): Letto1_meteo: avg(RH (%))
  flag 1 (5.22%): somero_meteo: Relative humidity
  flag 2 (0.03%): jokioinen_meteo: Relative humidity
  flag 3 (2.57%): Letto1_meteo_gapfilled: Paik RH
  flag 4 (0.00%): linearly interpolated
Rg: Global radiation [W/m2]
  flag 0 (92.18%): Letto1_meteo: avg(Glob (W/m2))
  flag 1 (3.84%): Letto1_meteo_gapfilled: PaikGlob2
  flag 2 (3.98%): jokioinen_rad: Global radiation
  flag 3 (0.00%): linearly interpolated
U: Wind speed [m/s]
  flag 0 (88.28%): Letto1_EC: wind speed (m/s) !u > 10 removed!
  flag 1 (7.23%): somero_meteo: Wind speed
  flag 2 (3.88%): Derived from salo_kiikala_meteo: wind speed (0.57*U_ref + 0.55)
  flag 3 (0.61%): linearly interpolated
Ustar: Friction velocity [m/s]
  flag 0 (86.42%): Letto1_EC: friction velocity (m/s)
  flag 1 (13.58%): Ustar = 0.2 * U
  flag 2 (0.00%): linearly interpolated
P: Ambient pressure [hPa]
  flag 0 (43.71%): Letto1_meteo: avg(Press (hPa))
  flag 1 (55.71%): Derived from salo_kiikala_meteo: Pressure (msl) (P_ref - 16.9)
  flag 2 (0.58%): linearly interpolated
  flag 1 (0.00%): filled with nearest
