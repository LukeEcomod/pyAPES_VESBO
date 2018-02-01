# -*- coding: utf-8 -*-
"""
GENERAL PARAMETERS
"""

gpara = {
        'dt' : 3600.0,  # timestep in forcing data file [s]
        'start_time' : "2013-06-01",  # start time of simulation [yyyy-mm-dd]
        'end_time' : "2018-01-01",  # end time of simulation [yyyy-mm-dd]
        'forc_filename' : "FMI_jokioinen.csv"  ## forcing data file*
        }


#  output_folder: "results"
#  output_format: "netcdf"

#  * Forcing data file (comma delimited, .csv)
#    Headers:
#        yyyy,mo,dd,hh,mm - year, month, day, hour, minute
#        Prec - Precipitation [mm/dt]
#        Tair - Air temperature [degC]
#        U - Wind speed 10 min avg. [m/s]
#        RH - Relative humidity [%]
#        Rg - Global radiation [W/m2]
