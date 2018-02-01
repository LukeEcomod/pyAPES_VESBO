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
Forc, Nsteps = read_forcing(gpara['forc_filename'],
                            gpara['start_time'],
                            gpara['end_time'])

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
    'LAIdecid': [],
    'Mbec': []
    }

""" RUN MODELS """
for k in range(0, Nsteps):
    print 'k = ' + str(k)

    # forcing canopy model
    doy = Forc['doy'].iloc[k]
    ta = Forc['Tair'].iloc[k]
    vpd = Forc['vpd'].iloc[k]
    rg = Forc['Rg'].iloc[k]
    par = Forc['Par'].iloc[k]
    prec = Forc['Prec'].iloc[k] * 1000
    u = Forc['U'].iloc[k]

    # Soil moisture information --- FOR ROOTING ZONE???
    Rew = 1.0  # np.minimum((smodel.Wliq - smodel.Wp) / (smodel.Fc - smodel.Wp + eps), 1.0)???
    beta = 1.0  # smodel.Wliq / smodel.poros??

    """ Canopy and Snow """
    # run canopy model
    cm_flx, cm_state = cmodel._run(doy, dt0, ta, prec, rg, par,
                                   vpd, U=u, beta=beta, Rew=Rew, P=101300.0)

    """ Water and Heat """
    # Potential infiltration and evaporation from ground surface
    ubc_w = {'Prec': float(cm_flx['PotInf']) / 1000 / dt0, 'Evap': float(cm_flx['Efloor']) / 1000 / dt0}
    # Transpiration sink
    rootsink = np.zeros(smodel.Nlayers)
    rootsink[0] = float(cm_flx['Transp']) / 1000 / dt0 / smodel.dz[0]  # ekasta layerista, ei väliä nyt kun tasapaino..
    # Temperature above soil surface
    ubc_T = {'type': 'flux', 'value': ta}
    # run soil water and heat flow
    sm_flx, sm_state = smodel._run(dt0, ubc_w, ubc_T, water_sink=rootsink)

    # -- update results (not here?)
    res['Prec'].append(prec * dt0)
    res['PotInf'].append(cm_flx['PotInf'])
    res['Interc'].append(cm_flx['Interc'])
    res['Trfall'].append(cm_flx['Trfall'])
    res['CanEvap'].append(cm_flx['CanEvap'])
    res['Efloor'].append(cm_flx['Efloor'])
    res['Transp'].append(cm_flx['Transp'])
    res['SWE'].append(cm_state['SWE'])
    res['LAI'].append(cm_state['LAI'])
    res['LAIdecid'].append(cm_state['LAIdecid'])
    res['Mbec'].append(cm_flx['MBE'] / 1000)
    res['gwl'].append(sm_state['ground_water_level'])
    res['h_pond'].append(sm_state['pond_storage'])
    res['Infil'].append(1000 * sm_flx['infiltration'])
    res['Efloorw'].append(sm_flx['evaporation'])
    res['Transw'].append(sm_flx['transpiration'])
    res['Drain'].append(1000 * sm_flx['drainage'])
    res['Roff'].append(1000 * sm_flx['runoff'])
#    res['WSto'].append(state['column_water_content'])
    res['Mbew'].append(sm_flx['water_closure'])
#    
#    if res['gwl'][-1] > 0.05 or any(np.isnan(res['h'][-1])) or res['h_pond'][-1] < 0.0 or abs(res['Mbew'][-1]) >1e-5:
#        break
    
plt.figure()
plt.subplot(5,1,1)
plt.plot(Forc.index, np.cumsum(res['Prec']), label='Prec')
plt.plot(Forc.index, np.cumsum(res['Transp']), label='Transpiration (c)')
plt.plot(Forc.index, 1000 * np.cumsum(res['Transw']), label='Transpiration (s)')
plt.plot(Forc.index, np.cumsum(res['Efloor']), label='Efloor')
plt.plot(Forc.index, np.cumsum(res['PotInf']), label='PotInf')
plt.plot(Forc.index, np.cumsum(res['Trfall']), label='Trfall')
plt.plot(Forc.index, np.cumsum(res['Roff']), label='Surface runoff')
plt.plot(Forc.index, np.cumsum(res['Drain']), label='Drainage runoff')
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.legend(loc=2, fontsize=8)
plt.subplot(5,1,2)
plt.plot(Forc.index, res['Transp'], 'g', linewidth=0.5, label='Transp')
plt.plot(Forc.index, res['Efloor'], linewidth=0.5, label='Efloor')
plt.plot(Forc.index, np.add(res['Drain'], res['Roff']) / dt0 * 3600, '-r', label='Total')
plt.ylabel('(mm h$^{-1}$)', fontsize=8)
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.subplot(5,1,3)
plt.plot(Forc.index, res['gwl'], label='gwl')
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.ylim([-1.0, 0.0])
plt.subplot(5,1,4)
plt.plot(Forc.index, res['SWE'], label='SWE')
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])
plt.subplot(5,1,5)
plt.plot(Forc.index, res['Mbec'], label='Mbec')
plt.plot(Forc.index, res['Mbew'], label='Mbew')
plt.legend(loc=2, fontsize=8)
plt.xlim([Forc.index[0], Forc.index[-1]])

##days = np.arange(0, k) * dt0 / 3600 /24
#plt.figure(2)
#plt.plot(res['gwl'])
#plt.title('Ground water level (m)', fontsize=14)
#plt.figure(3)
#plt.plot(res['Mbew'])
#plt.title('Water model mass balance error (m)', fontsize=14)
#plt.figure(4)
#plt.plot(1000 * np.add(res['Drain'], res['Roff']) / dt0 * 3600)
#plt.title('Total runoff (mm/h)', fontsize=14)
#
#print 'Water balance:'
#print('Precipition:\t %.2f mm (%.1f)' % (sum(res['Prec'])*1000, 100*sum(res['Prec'])/sum(res['Prec'])))
#print('Evaporation:\t %.2f mm (%.1f)' % (sum(res['Evap'])*1000, 100*sum(res['Evap'])/sum(res['Prec'])))
#print('Transpition:\t %.2f mm (%.1f)' % (sum(res['Trans'])*1000, 100*sum(res['Trans'])/sum(res['Prec'])))
#print('Drainage:\t %.2f mm (%.1f)' % (sum(res['Drain'])*1000, 100*sum(res['Drain'])/sum(res['Prec'])))
#print('Surface runoff:\t %.2f mm (%.1f)' % (sum(res['Roff'])*1000, 100*sum(res['Roff'])/sum(res['Prec'])))
#print('Storage change:\t %.2f mm (%.1f)' % ((res['WSto'][-1]-res['WSto'][0])*1000, 100*(res['WSto'][-1]-res['WSto'][0])/sum(res['Prec'])))

print('--- running time %.2f seconds ---' % (time.time() - running_time))
