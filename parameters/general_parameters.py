# -*- coding: utf-8 -*-
"""
GENERAL PARAMETERS
"""

gpara = {
        'dt' : 1800.0,  # timestep in forcing data file [s]
        'start_time' : "2010-01-01",  # start time of simulation [yyyy-mm-dd]
        'end_time' : "2010-12-31",  # end time of simulation [yyyy-mm-dd]
        'forc_filename' : "Hyde_2010_2016.csv"  #"FMI_jokioinen.csv"  ## forcing data file*
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
