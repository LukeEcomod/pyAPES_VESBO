# -*- coding: utf-8 -*-

"""
Main program for running CCFPeat

"""

import numpy as np
from matplotlib import pyplot as plt
from forcing.forc_utils import read_forcing
import canopy.canopy_model as cm
import soilprofile.soil_model as sm
import time  # For keeping track of computational time
running_time = time.time()

""" PARAMETERS """
# Import general parameters
from parameters.general_parameters import gpara
# Time step
dt0 = gpara['dt']
# Import canopy model parameters
from parameters.canopy_parameters import cpara
# Import soil model parameters
from parameters.soil_parameters import spara

""" FORCING DATA """
# Read forcing
Forc = read_forcing(gpara['forc_filename'],
                    gpara['start_time'],
                    gpara['end_time'],
                    dt=dt0)

""" INITIALIZE MODELS """
# Initialize canopy model
cmodel = cm.CanopyModel(cpara)
# Initialize soil model
smodel = sm.SoilModel(spara['z'], spara)

# Create results dictionary (not here?)
res = {
    'gwl': [],
#    'h': [],
#    'Wliq': [],
    'h_pond': [],
    'Infil': [],
    'Efloorw': [],
    'Transw': [],
    'Drain': [],
    'Roff': [],
#    'WSto': [],
    'Mbew': [],
    'Prec': [],
    'PotInf': [],
    'Interc': [],
    'Trfall': [],
    'CanEvap': [],
    'Efloor': [],
    'Transp': [],
    'SWE': [],
    'LAI': [],
    'Mbec': [],
    'Phenof': [],
    }

""" RUN MODELS """

for k in range(0, len(Forc)):
    print 'k = ' + str(k)

    # Soil moisture forcing for canopy model--- ROOTING ZONE???
    Rew = 1.0  # np.minimum((smodel.Wliq - smodel.Wp) / (smodel.Fc - smodel.Wp + eps), 1.0)???
    beta = 1.0  # smodel.Wliq / smodel.poros??

    """ Canopy and Snow """
    # run daily loop (phenology and seasonal LAI)
    if Forc['doy'].iloc[k] != Forc['doy'].iloc[k-1]:
        cmodel._run_daily(Forc['doy'].iloc[k], Forc['Tdaily'].iloc[k])
    # run timestep loop
    cm_flx, cm_state = cmodel._run_timestep(dt=dt0,
                                            forcing=Forc.iloc[k],
                                            beta=beta,
                                            Rew=Rew)

    """ Water and Heat """
    # potential infiltration and evaporation from ground surface
    ubc_w = {'Prec': cm_flx['PotInf'] / dt0, 'Evap': cm_flx['Efloor'] / dt0}
    # transpiration sink
    rootsink = np.zeros(smodel.Nlayers)
    rootsink[0] = cm_flx['Transp'] / dt0 / smodel.dz[0]  # ekasta layerista, ei väliä tasapainolaskennassa..
    # temperature above soil surface
    ubc_T = {'type': 'flux', 'value': None}
    # run soil water and heat flow
    sm_flx, sm_state = smodel._run(dt0, ubc_w, ubc_T, water_sink=rootsink)

    # -- update results (not here?)
    res['Prec'].append(Forc['Prec'].iloc[k] * dt0)
    res['PotInf'].append(cm_flx['PotInf'])
    res['Interc'].append(cm_flx['Interc'])
    res['Trfall'].append(cm_flx['Trfall'])
    res['CanEvap'].append(cm_flx['CanEvap'])
    res['Efloor'].append(cm_flx['Efloor'])
    res['Transp'].append(cm_flx['Transp'])
    res['SWE'].append(cm_state['SWE'])
    res['LAI'].append(cm_state['LAI'])
    res['Mbec'].append(cm_flx['MBE'])
    res['gwl'].append(sm_state['ground_water_level'])
    res['h_pond'].append(sm_state['pond_storage'])
    res['Infil'].append(sm_flx['infiltration'])
    res['Efloorw'].append(sm_flx['evaporation'])
    res['Transw'].append(sm_flx['transpiration'])
    res['Drain'].append(sm_flx['drainage'])
    res['Roff'].append(sm_flx['runoff'])
    res['Mbew'].append(sm_flx['water_closure'])
    res['Phenof'].append(cm_state['Phenof'])

#    if res['gwl'][-1] > 0.05 or any(np.isnan(res['h'][-1])) or res['h_pond'][-1] < 0.0 or abs(res['Mbew'][-1]) >1e-5:
#        break

print('--- running time %.2f seconds ---' % (time.time() - running_time))
