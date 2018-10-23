#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
.. module: bryotype
    :synopsis: APES-model component
.. moduleauthor:: Antti-Jussi Kieloaho

Bryotype describes structural and functional properties and processes of
moss/lichen (bryophytes) species/groups at the forest bottom layer. Original
implementation is done in MatLab by Samuli Launiainen.

Created on Tue Mar 14 08:15:50 2017

TODO:
    only water_balance calculation.

References:


CHANGES (13.-14.7.2017 SL):
    * BryoType: added state variables 'hydraulic_conductivity',
        'thermal_conductivity', 'carbon_pool'
    * changes names of some functions ('estimate_' --> 'moss_'),
        some simplifications of scrips
    * energy_and_water_balance: returns separate dicts 'fluxes' & 'states'
    * arg 'time' --> 'dt' in all calls
    * _run_timestep(args): simple entry-point to energy_and_water_balance
        and carbon_exchange
"""

from sys import version_info
import numpy as np

from heat_and_water import heat_and_water_exchange, evaporation_through
from heat_and_water import convert_hydraulic_parameters
from carbon import carbon_exchange

from canopy.constants import WATER_DENSITY, MOLAR_MASS_H2O, MOLAR_MASS_C

EPS = np.finfo(float).eps  # machine epsilon


class Bryophyte(object):
    r""" Represents bryophyte community-soil-atmosphere interactions.

    Characteristics of BryoType object are stored in 'properties' dictionary.
    These describes physical characteristics of a bryophyte layer.
    """

    # pylint: disable=too-many-instance-attributes
    # instance attributes are necessary to functioning of BryoModel

    def __init__(self, properties, initial_conditions=None):
        r""" Initialises a bryophyte object by using bryophyte's properties and
        initial states.

        Volumetric water content, relative water content is assumed to
        be equal to maximal retention capacity at field capacity.

        Leaf area index (*LAI*) is calculated as follows

        .. math::
            LAI = 1\\mathrm{e}^{3} \\frac{m_{dry} SLA}{1\\mathrm{e}^{4}}

        where :math:`m_{dry}` is dry mass and *SLA* is specific leaf area.

        Args:
            properties (dict):
                'ground_covarege':  [-]
                'height': [m]
                'roughness_height': [m]
                'leaf_area_index': [m\ :sup:`2` m :sup:`-2`\ ]
                'specific_leaf_area': [m\ :sup:`3` m :sup:`-3`\ ]
                'dry_mass': [kg m\ :sup:`-2`]
                'bulk_density': [kg m\ :sup:`-3`]
                'max_water_content': [g g\ :sup:`-1`\ ]
                'min_water_content': [g g\ :sup:`-1`\ ]
                'water_retention' (dict):
                    'theta_s': saturated water content [m\ :sup:`3` m :sup:`-3`\ ]
                    'theta_r': residual water content [m\ :sup:`3` m :sup:`-3`\ ]
                    'alpha': air entry suction [cm\ :sup:`-1`]
                    'n': pore size distribution [-]
                    'saturated conductivity': [m s\ :sup:`-1`]
                    'pore connectivity': (l) [-]
                    'compressability': 
                'porosity': [m\ :sup:`3` m\ :sup:`-3`\ ]
                'photosynthesis' (list):
                    0. Amax [\ :math:`\mu`\ mol m\ :sup:`-1`\ :sub:`leaf` s\ :sup:`-1`]
                    1. b [mol mol\ :sup:`-1`]
                'respiration' (list):
                    0. Q10 [-]
                    1. R10 [\ :math:`\mu`\ mol m\ :sup:`-1`\ :sub:`leaf` s\ :sup:`-1`]
                'optical_properties' (dict):
                    'albedo_PAR': [-] photosynthetically active radiation (PAR)
                    'albedo_NIR': [-] near infrared radiation (NIR)
                    'emissivity': [-] 

            initial_conditions (dict):
                initial_temperature: [\ :math:`^{\circ}`\ C]
                initial_water_content: [g g\ :sup:`-1`\ ]
        """

        dry_mass = properties['bulk_density'] * properties['height']
        properties['dry_mass'] = dry_mass

        residual_water_content = (properties['min_water_content']
                                  / WATER_DENSITY
                                  * properties['bulk_density'])

        field_capacity = (properties['max_water_content']
                          / WATER_DENSITY
                          * properties['bulk_density'])

        saturated_water_content = properties['porosity']

        if 'water_retention' not in properties:
            water_retention = {}

            water_retention['theta_r'] = residual_water_content
            water_retention['theta_s'] = saturated_water_content
            water_retention['field_capacity'] = field_capacity

            fraction = (
                    (saturated_water_content - field_capacity)
                    / (saturated_water_content - residual_water_content))

            water_retention['fraction'] = fraction
            water_retention['saturated_conductivity'] = properties['saturated_conductivity']

            water_retention['alpha'] = properties['alpha']
            water_retention['n'] = properties['n']

            if 'pore_connectivity' in properties:
                water_retention['pore_connectivity'] = properties['pore_connectivity']

            if 'compressibility' in properties:
                water_retention['compressability'] = properties['compressability']

            properties['water_retention'] = water_retention

        else:
            water_retention = properties['water_retention']
            water_retention['field_capacity'] = field_capacity

        self.properties = properties

        # set initial conditions
        if initial_conditions is not None:
            #: [:math:`^{\circ}`\ C]
            self.temperature = initial_conditions['temperature']

            #: [g g\ :sup:`-1`\ ]
            if initial_conditions['water_content'] <= properties['max_water_content']:
                self.water_content = initial_conditions['water_content']

            else:
                self.water_content = properties['max_water_content']

        else:
            #: [:math:`^{\circ}`\ C]
            self.temperature = 10.

            #: [g g\ :sup:`-1`\ ]
            self.water_content = (
                    properties['max_water_content']
                    + properties['min_water_content']) / 2.0

        self.water_storage = self.water_content * properties['dry_mass']

        #: [m\ :sup:`3` m\ :sup:`-3`\ ]
        self.volumetric_water = (
            self.water_content / WATER_DENSITY * properties['bulk_density'])

        #: [m]
        self.water_potential = convert_hydraulic_parameters(
                self.volumetric_water,
                self.properties['water_retention'],
                'volumetric_water')

        #: [kg C m-2], 'free carbon pool'
        self.carbon_pool = 0.0
        self.coverage = properties['ground_coverage']

        self.old_carbon_pool = self.carbon_pool
        self.old_water_content = self.water_content
        self.old_water_storage = self.water_storage
        self.old_volumetric_water = self.volumetric_water
        self.old_water_potential = self.water_potential
        self.old_temperature = self.temperature

    def update(self):
        """ Updates old states to states after iteration.
        """

        self.old_carbon_pool = self.carbon_pool
        self.old_water_content = self.water_content
        self.old_water_storage = self.water_storage
        self.old_volumetric_water = self.volumetric_water
        self.old_water_potential = self.water_potential
        self.old_temperature = self.temperature

    def restore(self):
        """ Restores new states back to states before iteration.
        """

        self.carbon_pool = self.old_carbon_pool
        self.water_content = self.old_water_content
        self.water_storage = self.old_water_storage
        self.volumetric_water = self.old_volumetric_water
        self.water_potential = self.old_water_potential
        self.temperature = self.old_temperature

    def run(self, dt, forcing, solver=False):
        r""" Calculates one timestep and updates states of BryoModel instance.

        Args:
            dt: timestep [s]
            forcing (dict): states of microclimate
                'throughfall': [mm s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ]
                'lwdn': [W m\ :sup:`-2`\ ]
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'precipitation_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_depth': [m]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [m]
                'soil_hydraulic_conductivity': [m s\ :sup:`-1`\ ]
                'soil_thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ]
                'nsteps' number of steps in odesolver

        Returns:
            fluxes (dict)
            states (dict)
        """

        # calculate moss energy and water balance and new state
        fluxes, states = heat_and_water_exchange(self.properties,
                                                 self.old_temperature,
                                                 self.old_water_content,
                                                 dt,
                                                 forcing,
#                                                 solver=solver
                                                 )

        # update state variables
        self.temperature = states['temperature']
        self.water_content = states['water_content']
        self.water_storage = states['water_storage']
        self.volumetric_water = states['volumetric_water']
        self.water_potential = states['water_potential']

        # solve photosynthesis and respiration

        # [umol m-2(ground) s-1]
        cflx = carbon_exchange(self.properties,
                               self.water_content,  # old
                               self.temperature,  # old
                               forcing['par'])

        nee = -cflx['photosynthesis_rate'] + cflx['respiration_rate']
        fluxes.update({'photosynthesis_rate': cflx['photosynthesis_rate'],
                       'respiration_rate': cflx['respiration_rate'],
                       'nee': nee})

        # update bryophyte free carbon pool (g C m-2) of bryophyte
        self.carbon_pool = self.old_carbon_pool + 1e3 * MOLAR_MASS_C * 1e-6 * nee
        states['carbon_pool'] = self.carbon_pool

        # compute soil evaporation through moss layer

        # [mol mol-1] -> [Pa]
        h2o = forcing['h2o'] * forcing['air_pressure']

        # [mol m-2 s-1]
        soil_evaporation = evaporation_through(
            self.properties,
            self.volumetric_water,  # old
            self.temperature,  # old
            forcing['air_temperature'],
            h2o,
            forcing['wind_speed'],
            forcing['air_pressure'],
            forcing['soil_temperature'],
            forcing['soil_water_potential'],
            forcing['soil_hydraulic_conductivity'],
            forcing['depth'])

        # unit conversion: 1000 kg m-2 s-1 = mm s-1

        if version_info.major < 3:
            soil_evaporation = {key: value * MOLAR_MASS_H2O for (key, value) in soil_evaporation.items()}

        else:
            soil_evaporation = {key: value * MOLAR_MASS_H2O for (key, value) in soil_evaporation}

        fluxes.update(soil_evaporation)

        return fluxes, states

class MossLayer():
    """ Simple moss layer for testing
    """

    def __init__(self, para):
        """ Moss layer interception, evaporation and CO2 exchange model
        """
        self.f_cover = para['ground_coverage']  # fraction of moss ground coverage [-]
        self.LAI = para['LAI']  # leaf area index
        self.Amax = para['Amax']  # max photo rate [umolm-2s-1]
        self.b = self.Amax / (2.0 * para['qeff'])  # half-saturation par
        self.R10 = para['R10']  # base respiration at 10degC
        self.Q10 = para['Q10']  # temperature sensitivity [-]
        
        self.zr = para['zr']  # roughness height m
        self.Mdry = para['Mdry']
        self.Wmax = para['Mdry']*para['Wmax']
        self.Wmin = para['Mdry']*para['Wmin']

        self.W = para['Wmax']*para['Mdry']      # current water content
        self.Wold = self.W

    def waterbalance(self, dt, Prec, U, T, H2O, P=101300.0):
        """
        Moss layer interception, evaporation and water balance.
        Args:
            dt - timestep [s]
            Prec - precipitation [mm]
            U - wind speed [m s-1]
            T - air temperature [degC]
            H2O - mixing ratio [mol mol-1]
            P - ambient pressure [Pa]
        Returns:
            Trfall - trfall rate below moss layer [mm]
            Evap - evaporation rate [mm/s]
            updates self.W
        """
        # VPD at air temperature; neglect condensation conditions
        es, _ = e_sat(T)
        D = np.maximum(0.0, es / P - H2O)  # mol / mol

        # initial water content
        Wo = self.Wold

        # interception and throughfall rate, new storage
        Ir = np.maximum(0.0, np.minimum(Prec, self.Wmax - Wo))
        Trfall = Prec - Ir  # mm

        W = Wo + Ir  # intermediate storage mm

        # evaporation from moss layer: actual conductance is boundary layer x
        # correction for internal resistance
        grel = np.minimum(0.1285 * W / self.Mdry - 0.1285, 1.0)
        gb = grel * self._boundary_layer_conductance(U)

        erate = gb * D  # mol m-2 s-1
        # rate = 1.26*eq_evap(Rn, T, units='mol')  # unrestricted rate
        Evap = np.minimum(erate * MOLAR_MASS_H2O * dt, W - self.Wmin)  # mm
        self.W = W - Evap  # mm

        Mbe = (self.W - Wo) - (Prec - Evap - Trfall)
        # print('Mbe', Mbe)

        return Evap/dt, Trfall, Mbe

    def co2_exchange(self, Par, T):
        """
        moss photosynthetic rate umolm-2s-1
        Args:
            Par (umolm-2s-1)
            T (degC)
        Returns:
            net photosynthetic rate (umolm-2s-1)
        """
        # Williams and Flanagan (1996),Oecologia 108, 38-46. Frolking et al. 1996 GCB
        a = [6.4355, -14.0605, 9.1867, -0.8720]
        b = [-4.3e-5, -8.3e-4, 0.08, 0.1]

        wn = self.W / self.Wmax

        # moisture response, always keep least 5% of capacity
        fW = np.maximum(0.05, a[3] + a[2]*wn + a[1]*wn**2.0 + a[0]*wn**3.0)

        # temperature response
        fT = b[0]*T**3.0 + b[1]*T**2.0 + b[2]*T + b[3]

        # compute photosynthetic rate [umol m-2 s-1]. Slice LAI into 10 layers, attenuate Par
        # exponentially and sum up layerwise photos. rates
        L = np.linspace(0, self.LAI, 10)
        dL = L[1] - L[0]
        Qp = Par*np.exp(-0.7*L)
        Ab = - fW * fT * np.sum(dL * (Qp / (Qp + self.b)))

        del fT, fW

        # respiration rate [umol m-2 s-1]
        if self.W <= 7.0:
            fW = -0.45 + 0.4*self.W - 0.0273*self.W**2.0
        else:
            fW = -0.04*self.W + 1.38

        fW = np.maximum(0.01, np.minimum(1.0, fW))

        Rb = self.R10 * self.Q10 ** ((T - 10.0) / 10.0) * fW

        return Ab + Rb

    def _boundary_layer_conductance(self, U):
        """ Moss boundary layer conductance as in Rice et al. 2001 eq. 1
        Args:
            zr - roughness lenght scale [m]
            U - mean wind speed [m s-1]
        Returns:
            gb - boundary layer conductance for H2O [mol m-2 s-1]
        """

        Dv = 24e-6  # m2s-1  molecular diffusitity at 20degC
        mu = 15.1e-6  # m2s-1 viscosity of air
        Sc = mu / Dv  # 0.63  # [-] ratio of viscosity to diffusivity
        rhoa = 41.6  # molm-3, density of air

        Re = U*self.zr / mu  # [-], Reynolds number

        gb = rhoa * 10**(-3.18) * Re**1.61 * Dv / self.zr * Sc**(0.33)  # m s-1

        return gb + eps

    def update(self):
        self.Wold = self.W


def test_bryomodel(fstep, nstep, param, forcing, odesteps=500, solver=False):
    """ this is for testing BryoModel stand-alone

    needs to access soilprofile to calculate:
        - soil thermal conductivity
        - soil hydraulic conductivity
        - soil water potential

    """

    import pandas as pd
    import soilprofile.soil_water as sw
    import soilprofile.soil_heat as sh

    from heat_and_energy import saturation_vapor_pressure

    columns = ['carbon_pool',
               'hydraulic_conductivity',
               'temperature',
               'thermal_conductivity',
               'volumetric_water_content',
               'water_content',
               'water_potential',
               'net_radiation_balance',
               'latent_heat_flux',
               'sensible_heat_flux',
               'ground_heat_flux',
               'emitted_longwave_radiation',
               'water_storage_change',
               'heat_storage_change',
               'interception',
               'throughfall_rate',
               'capillary_rise',
               'water_closure',
               'energy_closure']

    bryo_results = pd.DataFrame(index=forcing.index, columns=columns)

    dt = 1800.0

    result_list = []

    bryo = BryoModel(param)

    print("Wind speed is set to be 5% of forcing value!")

    pond_water_potential = 0.0  #1

    for k in range(fstep, fstep + nstep):

        wliq = forcing.iloc[k]['Wh']
#        wliq = 0.8889

        soil_thermal_conductivity = sh.thermal_conductivity_deVries(
            poros=0.89,
            wliq=wliq,
            T=forcing.iloc[k]['Tsh'],
            vOrg=0.11)

        soil_hydraulic_conductivity = sw.hydraulic_conductivity(
            pF={'alpha': 4.556640738735543,
                'n': 1.3112324995868292,
                'ThetaR': 0.074,
                'ThetaS': 0.91},
            x=wliq,
            var='Th',
            Ksat=2.42e-05)

        soil_water_potential = sw.wrc(
            pF={'alpha': 4.556640738735543,
                'n': 1.3112324995868292,
                'ThetaR': 0.074,
                'ThetaS': 0.91},
            x=wliq,
            var='Th')

        # compute H2O from relative humidity

#        if 'RH' in forcing.columns:
#            relative_humidity = forcing['RH'].iloc[k]
#
#        else:
#            relative_humidity = (
#                    forcing['h2o'].iloc[k]
#                    * 101300.0
#                    / saturation_vapor_pressure(forcing['Ta'].iloc[k]))

#            relative_humidity = h2o * air_pressure / svp
#            h_atm = (GAS_CONSTANT * (forc['air_temperature'] + DEG_TO_KELVIN)
#                     * np.log(rh) / (MOLAR_MASS_H2O*GRAVITY))


        par = forcing['diffPar'].iloc[k] + forcing['dirPar'].iloc[k]
        nir = forcing['diffNir'].iloc[k] + forcing['dirNir'].iloc[k]
        throughfall = forcing['Prec'].iloc[k]
        lwdn = forcing['LWin'].iloc[k]
        wind_speed = forcing['U'].iloc[k] * 0.05

        bryo_forcing = {
            'throughfall': throughfall,
            'air_temperature': forcing['Ta'].iloc[k],
            'soil_temperature': forcing['Tsh'].iloc[k],
            'soil_water_potential': soil_water_potential,
            'soil_depth': -0.01,
            'soil_hydraulic_conductivity': soil_hydraulic_conductivity,
            'soil_thermal_conductivity': soil_thermal_conductivity[0],
            'par': par,
            'nir': nir,
            'lwdn': lwdn,
            'wind_speed': wind_speed,
            'air_pressure': 101300.0,
            'h2o': forcing['H2O'].iloc[k],
            'nsteps': odesteps,
            'pond_water_potential': pond_water_potential
            }

        # compute bryophyte water, energy and carbon balances
        bryo_flx, bryo_state = bryo.run(dt=dt,
                                        forcing=bryo_forcing,
                                        solver=solver)

        bryo_state.update(bryo_flx)
        result_list.append(bryo_state)
        new_state = pd.Series(bryo_state)
        bryo_results.iloc[k] = new_state


#        pond_water_potential = max(pond_water_potential
#                                - bryo_state['pond_recharge'],
#                                0.0)

    # combine results into pandas dataframe

    df = pd.DataFrame.from_dict(result_list)
    df = df.set_index(forcing.index)

    return bryo_results, df
