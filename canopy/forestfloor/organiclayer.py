#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module: organiclayer
    :synopsis: APES-model component
.. moduleauthor:: Antti-Jussi Kieloaho and Samuli Launiainen

organiclayer describes structural and functional properties and processes of
organic layer (porous media) at forest floor bottom layer.
It can represnt moss/lichen (bryophyte) species/groups or litter layer.

Created on Tue Mar 14 08:15:50 2017.

Last edit 15.1.2020 / SL

Todo:
    - check behavior of water_exchange
    - handle conditions under snowcover
    - organic layer freezing/melting
    - solve surface temperature from energy balance (virtual 2-layer)
    - surface conductance parameterization: now based on classic log wind profile

    - respiration of litter layer: parameter values, dry mass or surface area-based?
        Add moisture response
    - Farquhar for bryophyte photosynthesis: surface-area or dry-mass based?
"""

import numpy as np

from canopy.constants import WATER_DENSITY, MOLAR_MASS_H2O, MOLAR_MASS_C, LATENT_HEAT, \
                             STEFAN_BOLTZMANN, DEG_TO_KELVIN, SPECIFIC_HEAT_AIR, \
                             SPECIFIC_HEAT_ORGANIC_MATTER, SPECIFIC_HEAT_H2O, GAS_CONSTANT, \
                             MOLECULAR_DIFFUSIVITY_CO2, MOLECULAR_DIFFUSIVITY_H2O, \
                             THERMAL_DIFFUSIVITY_AIR, AIR_DENSITY, AIR_VISCOSITY, GRAVITY

from .carbon import BryophyteFarquhar, OrganicRespiration

#: machine epsilon
EPS = np.finfo(float).eps

import logging
logger = logging.getLogger(self.__name__)

class OrganicLayer(object):
    r"""
    Universal Forestfloor Organic Object (a.k.a. combined Bryotype and Litter)
    """

    def __init__(self, properties):
        r""" Initialises object

        Volumetric water content, relative water content is assumed to
        be equal to maximal retention capacity at field capacity.

#        Leaf area index (*LAI*) is calculated as follows
#
#        .. math::
#            LAI = 1\\mathrm{e}^{3} \\frac{m_{dry} SLA}{1\\mathrm{e}^{4}}
#
#        where :math:`m_{dry}` is dry mass and *SLA* is specific leaf area.

        Args:
            properties (dict):
                'name': str
                'layer_type': 'bryophyte' or 'litter'
                'ground_coverage':  [-]
                'height': [m]
                'roughness_height': [m]
                #'leaf_area_index': [m\ :sup:`2` m :sup:`-2`\ ]
                #'specific_leaf_area': [m\ :sup:`3` m :sup:`-3`\ ]
                #'dry_mass': [kg m\ :sup:`-2`]
                'bulk_density': [kg m\ :sup:`-3`]
                'porosity': [m3 m-3]
                'max_water_content': [g g\ :sup:`-1`\ ]
                'min_water_content': [g g\ :sup:`-1`\ ]
                'water_retention' (dict):
                    #'theta_s': saturated water content [m\ :sup:`3` m :sup:`-3`\ ]
                    #'theta_r': residual water content [m\ :sup:`3` m :sup:`-3`\ ]
                    'alpha': air entry suction [cm\ :sup:`-1`]
                    'n': pore size distribution [-]
                    'saturated_conductivity': [m s\ :sup:`-1`]
                    'pore_connectivity': (l) [-]
                'porosity': [m\ :sup:`3` m\ :sup:`-3`\ ]
                'photosynthesis' (dict): : only if layer_type == 'bryophyte'
                    if farquhar-model as now:
                        'Vcmax', 'Jmax', 'Rd', 'alpha', 'theta', 'beta',
                        'gmax', 'wopt', 'a0', 'a1', 'CAP_desic', 'tresp'
                'respiration' (dict): only if layer_type == 'litter'
                    'q10' [-]
                    'r10' [\ :math:`\mu`\ mol m\ :sup:`-1`\ :sub:`leaf` s\ :sup:`-1`]
                'optical_properties' (dict):
                    'albedo' (dict)
                        'PAR', 'NIR'
                    'emissivity': [-]
                'initial_conditions' (dict):
                    temperature: [\ :math:`^{\circ}`\ C]
                    water_content: [g g\ :sup:`-1`\ ]
        """
        self.name = properties['name']
        self.layer_type = properties['layer_type']

        self.coverage = properties['coverage']
        self.height = properties['height']
        self.roughness_height = properties['roughness_height']
        #self.LAI = properties['leaf_area_index']
        #self.SLA = properties['specific_leaf_area']

        self.dry_mass = properties['bulk_density'] * properties['height']
        self.bulk_density = properties['bulk_density']
        self.porosity = properties['porosity']

        # hydraulic properties
        #self.max_water_content = properties['max_water_content']
        #self.min_water_content = properties['min_water_content']

        self.max_water_content = properties['max_water_content']
        self.max_symplast_water = self.max_water_content * properties['water_content_ratio']

        #self.max_symplast_water = properties['max_symplast_water_content']
        self.min_water_content = properties['min_water_content']

        self.water_retention = properties['water_retention']

        self.water_retention['theta_r'] = (self.max_symplast_water
                    / WATER_DENSITY * self.bulk_density)
        self.water_retention['theta_s'] = (self.max_water_content
                    / WATER_DENSITY * self.bulk_density)

        # optical properties at full saturation
        self.optical_properties = properties['optical_properties']

        # --- carbon exchange model
        if self.layer_type == 'bryophyte':
            # photosynthesis & respiration
            self.Carbon = BryophyteFarquhar(properties['photosynthesis'], carbon_pool=0.0)
        elif self.layer_type == 'litter':
            # only respiring
            self.Carbon = OrganicRespiration(properties['respiration'], carbon_pool=0.0)
        else:
            raise ValueError('Organiclayer init!: layer_type '
                + 'should be bryophyte or litter')

        # --- initial conditions
        initial_conditions = properties['initial_conditions']

        #: [:math:`^{\circ}`\ C]
        self.temperature = initial_conditions['temperature']
        self.surface_temperature = initial_conditions['temperature']
        #: [g g\ :sup:`-1`\ ]
        self.water_content = min(initial_conditions['water_content'], self.max_water_content)

        # [kg m-2]
        self.water_storage = self.water_content * self.dry_mass

        #: [m\ :sup:`3` m\ :sup:`-3`\ ]
        self.volumetric_water = (self.water_content / WATER_DENSITY * self.bulk_density)

        #: [m]
        self.water_potential = water_retention_curve(
            self.water_retention,
            theta=self.volumetric_water
        )
#        self.water_potential = convert_hydraulic_parameters(
#                self.volumetric_water,
#                self.water_retention,
#                'volumetric_water')

        self.albedo = reflectance(
            water_content=self.water_content,
            max_water_content=self.max_water_content,
            albedo=self.optical_properties['albedo']
        )
    
        self.emissivity = self.optical_properties['emissivity']

        #-- dict for temporal storage of object state after iteration
        self.iteration_results = None

    def update_state(self):
        """
        Updates object states after converged iteration. In pyAPES, called from canopy.canopy
        """
        self.temperature = self.iteration_results['temperature']
        self.surface_temperature = self.iteration_results['surface_temperature']
        self.water_storage = self.iteration_results['water_storage']
        self.water_content = self.iteration_results['water_content']
        self.volumetric_water = self.iteration_results['volumetric_water']
        self.water_potential = self.iteration_results['water_potential']

        #self.Carbon.carbon_pool = self.iteration_results['carbon_pool']

        self.albedo = reflectance(water_content=self.water_content,
                                  max_water_content=self.max_water_content,
                                  albedo=self.optical_properties['albedo'])

    def run(self, dt, forcing, parameters, controls):
        r""" Calculates one timestep and updates states of Bryophyte instance.

        Args:
            dt: timestep [s]
            forcing (dict):
                'throughfall': [kg m\ :sup:`-2`\ s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ] if energy_balance is True
                'lw_dn': [W m\ :sup:`-2`\ ] if energy_balance is True
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [m]
                'soil_pond_storage': [kg m-2]
                'snow_water_equivalent': [kg m-2]
            parameters (dict):
                'reference_height' [m]
                'soil_depth': [m]
                'soil_hydraulic_conductivity': [m s\ :sup:`-1`\ ]
                'soil_thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ]
                    if energy_balance is True
            controls (dict):
                'energy_balance': boolean
                'logger_info': str

        Returns:
            fluxes (dict)
            states (dict)
        """

        if controls['energy_balance']:
            # calculate moss energy and water balance
            fluxes, states = self.heat_and_water_exchange(
                                dt=dt,
                                forcing=forcing,
                                parameters=parameters,
                                sub_dt=60.0,
                                logger_info=controls['logger_info']
                                )

        else:
            # only water balance
            fluxes, states = self.water_exchange(
                                dt=dt,
                                forcing=forcing,
                                parameters=parameters
                                )

        #-- solve c02 exchange
        if self.layer_type == 'bryophyte':
            cflx = self.Carbon.co2_exchange(forcing['par'],
                                            forcing['co2'],
                                            states['temperature'],
                                            states['water_content']
                                            )

        else: # 'litter'
            cflx = self.Carbon.respiration(states['water_content'],
                                           states['temperature'])
            cflx.update({'internal_co2': -999.0,
                         'conductance_co2': -999.0
                         })

        # append to outputs
        fluxes.update(cflx)

        # kg m-2
        states['carbon_pool'] = self.Carbon.carbon_pool - 1e-6 * cflx['net_co2']* MOLAR_MASS_C

        # store iteration results
        self.iteration_results = states

        #-- compute soil evaporation through moss layer: diffusion through moss, then turbulent transport
        conductance_to_air = surface_atm_conductance(wind_speed=forcing['wind_speed'],
                                                     zref=parameters['reference_height'],
                                                     zom=self.roughness_height,
                                                     dT=0.0)
        # [mol m-2 s-1]
        soil_evaporation = evaporation_through_organic_layer(forcing,
                                                             conductance_to_air['h2o'],
                                                             states['volumetric_water'],
                                                             self.porosity,
                                                             self.height,
                                                             parameters)

        #soil_evaporation = 0.0
        # [kg m-2 s-1]
        fluxes['soil_evaporation'] = MOLAR_MASS_H2O * soil_evaporation

        return fluxes, states

    def heat_and_water_exchange(self,
                                dt,
                                forcing,
                                parameters,
                                sub_dt=60.0,
                                logger_info=''):
        r""" Solves organic layer coupled water and energy balance
        by Forward Euler integration

        Args:
            self: object
            dt: [s], timestep
            forcing (dict):
                'precipitation': [kg m\ :sup'-2' s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ]
                'lw_dn': [W m\ :sup:`-2`\ ]
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [Pa]
            parameters (dict):
                'soil_depth': [m]
                'soil_hydraulic_conductivity'
                'soil_thermal_conductivity'
                'reference_height' [m]
            sub_dt (float): internal timestep [s]
            logger_info: str

        Returns:
            fluxes (dict):
                'net_radiation' [W m\ :sup:`-2`\ ]
                'latent_heat' [W m\ :sup:`-2`\ ]
                'sensible_heat' [W m\ :sup:`-2`\ ]
                'ground_heat' [W m\ :sup:`-2`\ ] (negative towards soil)
                'heat_advection' [W m\ :sup:`-2`\ ]
                'water_closure' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'energy_closure' [W m\ :sup:`-2`\ ]
                'evaporation' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'interception' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'pond_recharge' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'capillary_rise' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'throughfall' [kg m\ :sup:`-2`\ s\ :sup`-1`\]

            states (dict):
                'temperature': [\ :math:`^{\circ}`\ C]
                'volumetric_water': [m\ :sup:`3` m\ :sup:`-3`\ ]
                'water_potential': [m]
                'water_content': [g g\ :sup:`-1`\ ]
                'water_storage': [kg m\ :sup:`-2`\ ]
                'hydraulic_conductivity': [m s\ :sup:`-1`\]
                'thermal_conductivity': [Wm-1K-1]
        """

        # initial conditions
        u = np.zeros(12) # 11
        u[0] = self.temperature # [degC]
        u[1] = self.water_storage # [kg m-2]

        # -- time loop
        t = 0.0
        while t < dt:
            # compute derivatives
            dudt, surface_temperature = self.water_heat_tendencies(u, sub_dt, forcing, parameters)
            # integrate in time
            u = u + dudt * sub_dt

            # advance in time
            t = t + sub_dt
            #sub_dt = min(t-dt, sub_dt)

        # unpack variables
        # [deg C]
        temperature = u[0]

        # [kg m-2]
        water_storage = u[1]

        # water fluxes [kg m-2 s-1]
        pond_recharge_rate = u[2] / dt
        capillary_rise = u[3] / dt
        interception_rate = u[4] / dt
        evaporation_rate = u[5] / dt

        # energy fluxes [W m-2] or [J m-2 s-1]
        net_radiation = u[6] / dt
        sensible_heat_flux = u[7] / dt
        latent_heat_flux = u[8] / dt
        conducted_heat_flux = u[9] / dt
        heat_advection = u[10] / dt
        ground_heat_flux = u[11] / dt

        # water balance closure [kg m-2 s-1] or [mm s-1]
        water_storage_change = (water_storage - self.water_content * self.dry_mass) / dt

        water_fluxes = (pond_recharge_rate
                        + capillary_rise
                        - evaporation_rate
                        + interception_rate)

        water_closure = -water_storage_change + water_fluxes

        # energy closure [W m-2] or [J m-2 s-1]
        # heat capacities [J m-2 K-1]
        C_old = (SPECIFIC_HEAT_ORGANIC_MATTER * self.dry_mass +
                 SPECIFIC_HEAT_H2O * self.water_content * self.dry_mass)
        C_new = (SPECIFIC_HEAT_ORGANIC_MATTER * self.dry_mass +
                 SPECIFIC_HEAT_H2O * water_storage)

        heat_storage_change = (C_new * temperature - C_old * self.temperature) / dt

        heat_fluxes = (net_radiation
                       - sensible_heat_flux
                       - latent_heat_flux
                       - ground_heat_flux
                       + heat_advection)

        energy_closure = -heat_storage_change + heat_fluxes

        # --- State variables of bryophyte layer
        # [g g-1]
        water_content = (water_storage / self.dry_mass)

        # [m3 m-3]
        volumetric_water = (water_content / WATER_DENSITY * self.bulk_density)
        # [m]
#        matrix_potential = convert_hydraulic_parameters(volumetric_water,
#                                                        self.water_retention,
#                                                        'volumetric_water')

        # must constrain theta as: theta = max(volumetric_water, pF['theta_r'])
        matrix_potential = water_retention_curve(self.water_retention, theta=volumetric_water)
        # [m s-1]
        Kliq = hydraulic_conductivity(self.water_retention, volumetric_water,)

        # [W m-1 K-1]
        Lambda = thermal_conductivity(volumetric_water)

        # return fluxes and state variables
        fluxes = {
            'net_radiation': net_radiation,  # [W m-2]
            'latent_heat': latent_heat_flux,  # [W m-2]
            'sensible_heat': sensible_heat_flux,  # [W m-2]
            'conducted_heat': conducted_heat_flux, # [W m-2]
            'ground_heat': ground_heat_flux,  # [W m-2]
            'heat_advection': heat_advection,  # [W m-2]
            'water_closure': water_closure,  # [mm s-1]
            'energy_closure': energy_closure,  # [W m-2]
            'interception': interception_rate,  # [mm s-1]
            'evaporation': evaporation_rate,  # [mm s-1]
            'pond_recharge': pond_recharge_rate,  # [mm s-1]
            'capillary_rise': capillary_rise,  # [mm s-1]
            'throughfall': forcing['precipitation'] - interception_rate  # [mm s-1]
            }

        states = {
            'volumetric_water': volumetric_water,  # [m3 m-3]
            'water_potential': matrix_potential,  # [m]
            'water_content': water_content,  # [g g-1]
            'water_storage': water_storage,  # [kg m-2] or [mm]
            'temperature': temperature,  # [degC]
            'surface_temperature': surface_temperature, # [degC]
            'hydraulic_conductivity': Kliq,  # [m s-1]
            'thermal_conductivity': Lambda,  # [W m-1 K-1]
            }

        return fluxes, states


    def water_heat_tendencies(self, y, dt, forcing, parameters):
        """
        Solves coupled water and heat balance over timestep dt.
        Returns temperature and water storage tendensies and
        water & energy fluxes.
        Args:
            y (array): initial values: returns derivatives
                0. temperature [degC] --> computes (du/dt)
                1. water storage [kg m-2] --> (du/dt)
                2. pond_recharge
                3. capillary_rise
                4. interception
                5. evaporation/condensation
                6. radiative_heat
                7. sensible_heat
                8. latent_heat
                9. conducted_heat
                10 advected_heat
                11. ground_heat

            dt (float): time step
            forcing
            parameters
        Returns:
            derivatives (du/dt) of y
            surface_temperature [deg C]

        Units:  temperature [K s-1]
                water content and water fluxes [kg m-2 s-1 = mm s-1]
                energy fluxes [J m-2 s-1 = W m-2]
        """
        logger = logging.getLogger(self.__name__)
        
        dudt = np.zeros(12)

        if dt == 0.0:
            dt = dt + EPS

        zm = 0.5 * self.height
        zs = abs(parameters['soil_depth'])
        zref = parameters['reference_height']

        max_storage = self.max_water_content * self.dry_mass
        min_storage = self.min_water_content * self.dry_mass
        symplast_storage = self.max_symplast_water * self.dry_mass

        # --- compute water balance ---
        # nomenclature:
        #   water_content [g g-1], water_storage [kg m-2 s-1 == mm]
        #   volumetric water content [-], water potential [m]
        #   all fluxes [kg m-2 s-1]

        # initial state

        # [kg m-2] or [mm]
        water_storage = min(max_storage, y[1])
        water_storage = max(min_storage, water_storage)

        # [g g-1]
        water_content = (water_storage / self.dry_mass)

        # [m m-3]
        volumetric_water = (water_content / WATER_DENSITY
                            * self.bulk_density)

        # [m]
        # theta is restricted above or equal Theta_r in water_retention_curve
        water_potential = water_retention_curve(self.water_retention,
                                                theta=volumetric_water
                                                )

        # --- constraints for recharge & evaporation
        max_recharge = max(max_storage - water_storage, 0.0)
        max_recharge = min(max_storage, max_recharge)

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max_recharge / dt

        # [kg m-2 s-1] or [mm s-1]
        max_evaporation_rate = (y[1] - (min_storage + EPS)) / dt

        if np.isinf(max_evaporation_rate) or max_evaporation_rate < 0.0:
            max_evaporation_rate = 0.0

        # [kg m-2 s-1] or [mm s-1]
        max_condensation_rate = -((max_storage - EPS) - y[1]) / dt

        if np.isinf(max_condensation_rate) or max_condensation_rate > 0.0:
            max_condensation_rate = 0.0

        #--- compute surface temperature from surface energy balance
        Ta = forcing['air_temperature']

        # take reflectivities from previous timestep
        # radiation balance # [J m-2 s-1] or [W m-2]

        # [J m-2 s-1] or [W m-2]
        SWabs = (1.0 - self.albedo['PAR']) * forcing['par'] + \
                (1.0 - self.albedo['NIR']) * forcing['nir']

        # conductance for heat from surface layer to moss
        # [W m-2 K-1]
        moss_thermal_conductivity = thermal_conductivity(volumetric_water)
        gms = moss_thermal_conductivity / zm

        # conductances from surface to air
        conductance_to_air = surface_atm_conductance(wind_speed=forcing['wind_speed'],
                                                     zref=zref,
                                                     zom=self.roughness_height,
                                                     dT=0.0)
        ga = conductance_to_air['heat'] # heat [mol m-2 s-1]
        gav = conductance_to_air['h2o'] #  [mol m-2 s-1]

        # linear decrease of evaporation when capillary water has been evaporated
        # relative_conductance = min(1.0, (0.1285 * y[1] / self.dry_mass - 0.1285))
        relative_conductance = min(1.0, (y[1] / symplast_storage) + EPS)
        gv = gav * relative_conductance # h2o, [mol m-2 s-1]

        # -- solve surface temperature Ts iteratively
        err = 999.0
        itermax = 50
        iter_no = 0

        #wo = 0.8 # weight of old Ts
        #
        ## SL 28.7.20
        #Ts = 0.5 * (Ta + self.temperature)

        Ts = Ta  # initial guess

        while err > 0.01 and iter_no < itermax:

            # evaporation demand and supply --> latent heat flux
            es = saturation_vapor_pressure(Ts) / forcing['air_pressure']
            LEdemand = LATENT_HEAT * gv * (es - forcing['h2o'])
            if LEdemand > 0:
                LE = min(LEdemand, LATENT_HEAT * max_evaporation_rate)
                if iter_no > 100:
                    LE = 0.1 * LATENT_HEAT * max_evaporation_rate
            else: # condensation
                LE = max(LEdemand, LATENT_HEAT * max_condensation_rate)

            Told = np.copy(Ts)

            # --- find Ts: Long-wave term is linearized as in Campbell & Norman 1998 Ch 12.
            #Te = Ta
            #gr = 4 * self.emissivity * STEFAN_BOLTZMANN * (Te + DEG_TO_KELVIN)**3 / SPECIFIC_HEAT_AIR
            #a = Rni - LE
            #b = SPECIFIC_HEAT_AIR * (ga + gr)
            #Ts = (a + b * Ta + gms * y[0]) / (b + gms)
            
            #LW_up linearized against Told (instead of Ta): eoT_s^4 ~= eoT_old^4 + 4eoT_old^3*Ts
            gr = 4 * self.emissivity * STEFAN_BOLTZMANN * (Told + DEG_TO_KELVIN)**3 / SPECIFIC_HEAT_AIR
            Rn = (SWabs
                  + self.emissivity * forcing['lw_dn'] - self.emissivity * STEFAN_BOLTZMANN * (Told + DEG_TO_KELVIN)**4)

            Ts = (Rn + SPECIFIC_HEAT_AIR * (gr * Told + ga * Ta) - LE + gms * y[0]) / (
                SPECIFIC_HEAT_AIR * (ga + gr) + gms)

            err = abs(Ts - Told)

            iter_no += 1

            if iter_no == itermax:
                #logger.debug('Maximum number of iterations reached: Ts = %.2f, err = %.2f', np.mean(Ts), err)

                Ts = Ta
                es = saturation_vapor_pressure(Ts) / forcing['air_pressure']
                LEdemand = LATENT_HEAT * gv * (es - forcing['h2o'])
                if LEdemand > 0:
                    LE = min(LEdemand, LATENT_HEAT * max_evaporation_rate)
                else: # condensation
                    LE = max(LEdemand, LATENT_HEAT * max_condensation_rate)

        # -- energy fluxes  [J m-2 s-1] or [W m-2]

        LWup = self.emissivity * STEFAN_BOLTZMANN * (Ts + DEG_TO_KELVIN)**4
        net_radiation = SWabs +  self.emissivity * forcing['lw_dn'] - LWup
        sensible_heat_flux = SPECIFIC_HEAT_AIR * ga * (Ts - Ta)
        conducted_heat_flux = gms * (Ts - y[0])
        latent_heat_flux = LE

        # -- evaporation rate [kg m-2 s-1]
        evaporation_rate = LE / LATENT_HEAT * MOLAR_MASS_H2O

        max_recharge_rate = max(max_recharge_rate + evaporation_rate, 0.0)
        max_recharge = max_recharge_rate * dt

        # -- recharge from rainfall interception and/or from pond storage
        # assumptotic function
        interception = (max_recharge *
                        (1.0 - np.exp(-(1.0 / max_storage)
                         * forcing['precipitation'] * dt)))

        # [kg m-2 s-1] or [mm s-1]
        interception_rate = interception / dt

        # [kg m-2] or [mm]
        max_recharge = max(max_recharge - interception, 0.0)

        pond_recharge = min(max_recharge - EPS, forcing['soil_pond_storage'])

        # [kg m-2 s-1] or [mm s-1]
        pond_recharge_rate = pond_recharge / dt

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max(max_recharge - pond_recharge, 0.0) / dt

        # --- compute capillary rise from soil [ kg m-2 s-1 = mm s-1]

        # Kh: soil-moss hydraulic conductivity assuming two resistors in series
        Km = hydraulic_conductivity(self.water_retention, volumetric_water)
        Ks = parameters['soil_hydraulic_conductivity']

        # conductance of layer [s-1]
        g_moss = Km / zm
        g_soil = Ks / zs

        # [m s-1]
        Kh = (g_moss * g_soil / (g_moss + g_soil)) * (zm + zs)

        capillary_rise = WATER_DENSITY * max(0.0, - Kh * ((water_potential - forcing['soil_water_potential'])
                                                          / (zm + zs) + 1.0))

        # [kg m-2 s-1] or [mm s-1]
        capillary_rise = min(capillary_rise, max_recharge_rate)

        # --- calculate mass balance of water

        # [kg m-2 s-1] or [mm s-1]
        dy_water = (
                interception_rate
                + pond_recharge_rate
                + capillary_rise
                - evaporation_rate
                )

        # --- calculate change in moss heat content

        # heat conduction between moss and soil
        # [W m-2 K-1]
        moss_thermal_conductivity = thermal_conductivity(volumetric_water)

        # thermal conductance [W m-2 K-1]
        # assume the layers act as two resistors in series
        g_moss = moss_thermal_conductivity / zm
        g_soil = parameters['soil_thermal_conductivity'] / zs

        thermal_conductance = (g_moss * g_soil) / (g_moss + g_soil)

        # [J m-2 s-1] or [W m-2]
        ground_heat_flux = thermal_conductance *(y[0] - forcing['soil_temperature'])

        # heat lost or gained with liquid water removing/entering
        # [J m-2 s-1] or [W m-2]
        heat_advection = SPECIFIC_HEAT_H2O * (
                        interception_rate * forcing['air_temperature']
                        + capillary_rise * forcing['soil_temperature']
                        + pond_recharge_rate * forcing['soil_temperature']
                        )

        # heat capacities
        # [J K-1]
        heat_capacity_old = (SPECIFIC_HEAT_ORGANIC_MATTER
                         * self.dry_mass
                         + SPECIFIC_HEAT_H2O * y[1])

        heat_capacity_new = (SPECIFIC_HEAT_ORGANIC_MATTER
                             * self.dry_mass
                             + SPECIFIC_HEAT_H2O * (y[1] + dy_water * dt))

        # calculate new temperature from heat balance
        heat_fluxes = (
                + conducted_heat_flux
                + heat_advection
                - ground_heat_flux
                )

        new_temperature = (heat_fluxes * dt + heat_capacity_old * y[0]) / heat_capacity_new

        # -- return tendencies and fluxes
        # [K s-1]
        dudt[0] = (new_temperature - y[0]) / dt
        # [kg m-2 s-1] or [mm s-1]
        dudt[1] = dy_water

        # water fluxes
        # [kg m-2 s-1] or [mm s-1]
        dudt[2] = pond_recharge_rate
        dudt[3] = capillary_rise
        dudt[4] = interception_rate
        dudt[5] = evaporation_rate

        # energy fluxes
        # [J m-2 s-1] or [W m-2]
        dudt[6] = net_radiation
        dudt[7] = sensible_heat_flux
        dudt[8] = latent_heat_flux
        dudt[9] = conducted_heat_flux
        dudt[10] = heat_advection
        dudt[11] = ground_heat_flux

        return dudt, Ts

    def water_exchange(self,
                       dt,
                       forcing,
                       parameters):
        """
        Args:
            dt: [s], timestep
            forcing (dict):
                'precipitation': [kg m\ :sup'-2' s\ :sup:`-1`\ ]
                'wind_speed': [m s\ :sup:`-1`\ ]
                'friction_velocity': [m s\ :sup:`-1`\ ]
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_pond_storage': ??
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [Pa]
            parameters (dict):
                'soil_depth': [m]
                'soil_hydraulic_conductivity'
                'soil_thermal_conductivity'
                'reference_height' [m]

        Returns:
            fluxes (dict):
                'latent_heat' [W m\ :sup:`-2`\ ]
                'sensible_heat' [W m\ :sup:`-2`\ ]
                'ground_heat' [W m\ :sup:`-2`\ ] (negative towards soil)
                'evaporation' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'interception' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'pond_recharge' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'capillary_rise' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'throughfall' [kg m\ :sup:`-2`\ s\ :sup`-1`\]

            states (dict):
                'temperature': [\ :math:`^{\circ}`\ C]
                'volumetric_water': [m\ :sup:`3` m\ :sup:`-3`\ ]
                'water_potential': [m]
                'water_content': [g g\ :sup:`-1`\ ]
                'water_storage': [kg m\ :sup:`-2`\ ]
                'hydraulic_conductivity': [m s\ :sup:`-1`\]
                'thermal_conductivity': [Wm-1K-1]
        """
        # moss is assumed to be at air temperature
        temperature = forcing['air_temperature']

        # [kg m-2] or [mm]
        water_storage = self.water_storage
        max_storage = self.max_water_content * self.dry_mass
        min_storage = self.min_water_content * self.dry_mass
        symplast_storage = self.max_symplast_water * self.dry_mass

        # change in water storage during dt [kg m-2]
        d_water_storage = 0.0

        # --- evaporation during dt [kg m-2]

        # boundary layer conductances for H2O, heat and CO2
        # [mol m-2 s-1]
        conductance_to_air = surface_atm_conductance(
            wind_speed=forcing['wind_speed'],
            zref=parameters['reference_height'],
            zom=self.roughness_height,
            dT=0.0)

        # water vapor conductance from moss to air
        relative_conductance = min(1.0, (water_storage / symplast_storage) + EPS)

        # [mol m-2 s-1]
        conductance_to_air_h2o = (
            conductance_to_air['h2o'] * relative_conductance
        )

        # [kg m-2]
        evaporation = (
            conductance_to_air_h2o * (611.0 / forcing['air_pressure']
            * np.exp(17.502 * temperature / (forcing['air_temperature'] + 240.0))
             - forcing['h2o']) * MOLAR_MASS_H2O * dt
        )
        # restricted by space available for condensation and water available for evaporation
        evaporation = min(max(water_storage - max_storage, evaporation),
                          water_storage - min_storage)
        d_water_storage -= evaporation

        #--- interception of rainfall and recharge from ponding water during dt
        # [mm] or [kg m-2]
        # interception = min(max_storage - (water_storage + d_water_storage),
        #                     forcing['precipitation'] * dt)
        # assumptotic function
        interception = ((max_storage - (water_storage + d_water_storage)) *
                        (1.0 - np.exp(-(1.0 / max_storage)
                          * forcing['precipitation'] * dt)))
        d_water_storage += interception

        pond_recharge = min(max_storage - (water_storage + d_water_storage),
                            forcing['soil_pond_storage'])
        d_water_storage += pond_recharge

        #--- capillary rise from underlying soil during dt [kg m-2]
        # water potential[m]
        water_potential = water_retention_curve(self.water_retention, self.volumetric_water)

        # hydraulic conductivity from soil to moss [m s-1]
        Km = hydraulic_conductivity(self.water_retention, self.volumetric_water)
        Ks = parameters['soil_hydraulic_conductivity']

        zm = 0.5 * self.height
        zs = abs(parameters['soil_depth'])

        # conductance of layer [s-1]
        g_moss = Km / zm
        g_soil = Ks / zs

        # [m s-1]
        Kh = (g_moss * g_soil / (g_moss + g_soil)) * (zm + zs)

        #[kg m-2]
        capillary_rise = (max(0.0, - Kh * (
            (water_potential - forcing['soil_water_potential'])
            / (zm + zs) + 1.0)) * WATER_DENSITY * dt)

        capillary_rise = min(max_storage - (water_storage + d_water_storage),
                             capillary_rise)
        d_water_storage += capillary_rise

        #--- new state
        water_storage = (water_storage + d_water_storage)  #[kg m-2]
        water_content = water_storage / self.dry_mass  # [g g-1]
        volumetric_water = (water_content / WATER_DENSITY * self.bulk_density)  # [m3 m-3]

        water_potential = water_retention_curve(self.water_retention, volumetric_water)

        Kliq = hydraulic_conductivity(self.water_retention, volumetric_water)

        # --- Heat exchange dummy variables ---

        # heat conduction between moss and soil [W m-2 K-1]
        moss_thermal_conductivity = thermal_conductivity(volumetric_water)

        # thermal conductance [W m-2 K-1]; assume the layers act as two resistors in series
        g_moss = moss_thermal_conductivity / zm
        g_soil = parameters['soil_thermal_conductivity'] / zs

        thermal_conductance = (g_moss * g_soil) / (g_moss + g_soil)

        # [J m-2 s-1] or [W m-2]
        ground_heat = thermal_conductance *(forcing['air_temperature'] - forcing['soil_temperature'])

        latent_heat = (LATENT_HEAT * (evaporation / dt + EPS) / MOLAR_MASS_H2O)

        sensible_heat = 0.0

        fluxes = {
            'evaporation': evaporation / dt,  # [mm s-1]
            'capillary_rise': capillary_rise / dt,  # [mm s-1]
            'pond_recharge': pond_recharge / dt,  # [mm s-1]
            'throughfall': forcing['precipitation'] - interception / dt,  # [mm s-1]
            'interception': interception / dt,
            'ground_heat': ground_heat,  # [W m-2]
            'latent_heat': latent_heat,  # [W m-2]
            'sensible_heat': sensible_heat,  # [W m-2]
        }

        states = {
            'volumetric_water': volumetric_water,  # [m3 m-3]
            'water_potential': water_potential,  # [m]
            'water_content': water_content,  # [g g-1]
            'water_storage': water_storage,  # [kg m-2] or [mm]
            'hydraulic_conductivity': Kliq,  # [m s-1]
            'thermal_conductivity': moss_thermal_conductivity,  # [W m-1 K-1]
            'temperature': temperature,  # [degC]
            'surface_temperature': temperature # [degC]
        }

        return fluxes, states

def reflectance(water_content, max_water_content, albedo):
        """
        Water-content depended bryophyte shortwave albedo [-] for PAR and NIR wavebands

        The effect of water content of spectral properties are based on studies by
        Vogelmann and Moss (1993) and Fernandes (1999)on Sphagnum cuspidatum and
        Pleurozium schreberi, respectively.

        The albedo is scaled specific reflectance for PAR (400-700 nm) or
        NIR (750-1400) regions. The scaling coefficient is common for both
        PAR and NIR and scales relative to fully hydrated bryophyte albedo

        References:
            Vogelmann and Moss (1993)
                Remote Sensing of Environment 45:273-279.
            Fernandes (1999)
                PhD thesis entitled: 'Scale influences of surface
                parametrization on modelled boreal carbon and water budgets'
        Arg:
            object
            water_content [g g-1]
            max_water_content [g g-1]
            albedo (dict): 'PAR','NIR' [-]
        Returns (dict):
            'PAR'
            'NIR'
        """

        x = water_content / max_water_content

        scaling_coefficient = 1.0 + (4.5 - 1.0) / (1.00 + np.power(10, 4.53 * x))

        albedo_nir = albedo['NIR'] * scaling_coefficient
        albedo_par = albedo['PAR'] * scaling_coefficient

        return {'PAR': albedo_par, 'NIR': albedo_nir}

def thermal_conductivity(volumetric_water, method='odonnel'):
    r""" Estimates thermal conductivity (km) of bryophyte layer.

    By default organic matter heat conductivity is calculated by using equation
    by O'Donnel et al. (2009, Soil Sci. 174).

    Args:
        volumetric_water: [m\ :sup:`3` m\ :sup:`-3`\ ]
        flag (optional):
            optional methods are:
                * 'campbell' (Campbell et al. (1985))
                * 'constant' (0.25)
                * 'lauren' (Lauren (1999, table 4))

    Returns:
        float: heat conductivity in [W m\ :sup:`-1` K\ :sup:`-1`\ ]
    """

    method = method.lower()

    heat_conductivity = None

    if method == 'odonnel':  # O'Donnell
        heat_conductivity = np.minimum(0.6,
                                       3.0e-2 + 5.0e-1 * volumetric_water)

        heat_conductivity = np.maximum(4.0e-2, heat_conductivity)

    elif method == 'campbell':
        heat_conductivity = (
            0.4 + 0.5
            * volumetric_water(0.4 - 0.06)
            * np.exp(-(1. * volumetric_water)**4))

    elif method == 'constant':
        heat_conductivity = 0.25

    elif method == 'lauren':
        heat_conductivity = -0.004 + 0.609 * volumetric_water

    # [W/(m K)]
    return heat_conductivity


def surface_atm_conductance(wind_speed, zref, friction_velocity=None, dT=0.0, zom=0.01, b=1.1e-3):
    """
    Soil surface - atmosphere transfer conductance for scalars. Two paralell
    mechanisms: forced and free convection
    Args:
        wind_speed - wind speed (m/s) at zref
        zref - reference height (m). Log-profile assumed below zref.
        zom - roughness height for momentum (m), ~ 0.1 x canopy height
        ustar - friction velocity (m/s) at log-regime. if ustar not given,
                it is computed from Uo, zref and zom
        b - parameter for free convection. b=1.1e-3 ... 3.3e-3 from smooth...rough surface
    Returns:
        conductances for CO2, H2O and heat (mol m-2 s-1), dict
    References:
        Schuepp and White, 1975:Transfer Processes in Vegetation by Electrochemical Analog,
        Boundary-Layer Meteorol. 8, 335-358.
        Schuepp (1977): Turbulent transfer at the ground: on verification of
        a simple predictive model. Boundary-Layer Meteorol., 171-186
        Kondo & Ishida, 1997: Sensible Heat Flux from the Earth’s Surface under
        Natural Convective Conditions. J. Atm. Sci.
    """

    Sc_v = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_H2O
    Sc_c = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_CO2
    Pr = AIR_VISCOSITY / THERMAL_DIFFUSIVITY_AIR
    kv = 0.4  # von Karman constant (-)
    d = 0.0 # displacement height

    if friction_velocity == None:
        friction_velocity = wind_speed * kv / np.log((zref - d) / zom)

    delta = MOLECULAR_DIFFUSIVITY_H2O / (kv*friction_velocity + EPS)

    gb_h = (kv * friction_velocity) / (Pr - np.log(delta / zref))
    gb_v = (kv * friction_velocity) / (Sc_v - np.log(delta / zref))
    gb_c = (kv * friction_velocity) / (Sc_c - np.log(delta / zref))

    # free convection as parallel pathway, based on Condo and Ishida, 1997.
    #b = 1.1e-3 #ms-1K-1 b=1.1e-3 for smooth, 3.3e-3 for rough surface
    dT = np.maximum(dT, 0.0)

    gf_h = b * dT**0.33  # ms-1

    # mol m-2 s-1
    gb_h = (gb_h + gf_h) * AIR_DENSITY
    gb_v = (gb_v + Sc_v / Pr * gf_h) * AIR_DENSITY
    gb_c = (gb_c + Sc_c / Pr * gf_h) * AIR_DENSITY

    return {'co2': gb_c, 'h2o': gb_v, 'heat': gb_h}

def evaporation_through_organic_layer(forcing,
                                      boundary_layer_conductance,
                                      volumetric_water,
                                      porosity,
                                      moss_height,
                                      parameters
                                      ):
    r""" Estimates soil evaporation rate through bryophyte layer.

    Evaporation in bryophyte layer is limited either by atmospheric demand
    and transport or soil supply of water.

    Water vapor flow from soil to air must overcome two resistances
    in series: 1. molecular diffusion through porous living moss, 2. molecular
    diffusion from moss canopy to 1st caluclation node in the atomosphere. The
    2nd resistance is assumed to be equal to conductance from wet moss canopy.

    Args:
        volumetric_water: [m\ :sup:`3` m\ :sup:`-3`\ ]
        temperature [degC]
        forcing (dict)
            air_temperature: [\ :math:`^{\circ}`\ C]
            h2o: [Pa]
            wind_speed: [m s\ :sup:`-1`\ ]
            air_pressure: [Pa]
            soil_temperature: [\ :math:`^{\circ}`\ C]
                from 1st calculation node
            soil_water_potential: [m]
                from 1st calculation node
        parameters (dict)
            reference_height [m]
            soil_hydraulic_conductivity: [m s\ :sup:`-1`\ ]
                from 1st calculation node
            soil_depth: [m]

    Returns:
        float:
            evaporation in [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
    """
    # [mol mol-1] -> [Pa]
    h2o = forcing['h2o'] * forcing['air_pressure']
    Ta = forcing['air_temperature']
    Ts = forcing['soil_temperature']
    Pamb = forcing['air_pressure']
    afp = porosity - volumetric_water

    #-- conductance for h2o "from soil surface through porous media"

    # [mol m-3], air molar density
    cair = Pamb / (GAS_CONSTANT * (Ta + DEG_TO_KELVIN))

    # D/Do, diffusivity in porous media relative to that in free air,
    # Millington and Quirk (1961)
    relative_diffusivity = (np.power(Ta / 293.16, 1.75) * np.power(afp, 10.0/3.0) / porosity**2)

    g_molecular = cair * MOLECULAR_DIFFUSIVITY_H2O * relative_diffusivity / moss_height

    # [mol m-2 s-1], two resistors in series
    conductance_h2o = (g_molecular * boundary_layer_conductance / (g_molecular + boundary_layer_conductance))

    # Assuming soil rh = 1, calculate the maximum evaporation rate     # [mol/(m2 s)]
    # atmospheric evaporative demand
    evaporative_demand = max(0.0, conductance_h2o *
                                 (saturation_vapor_pressure(Ts) - h2o) / Pamb
                            )

    #  soil supply
    # [- ]
    rh_air = min(1.0, h2o / saturation_vapor_pressure(Ta))

    # [m], in equilibrium with atmospheric relative humidity
    atm_hydraulic_head = (GAS_CONSTANT * (Ta + DEG_TO_KELVIN) * np.log(rh_air)
                        / (MOLAR_MASS_H2O * GRAVITY))

    # E = -Kh * [(ha - hs) / zs + 1.0]
    evaporative_supply = max(0.0,
        - WATER_DENSITY / MOLAR_MASS_H2O * parameters['soil_hydraulic_conductivity']
        * ((atm_hydraulic_head - forcing['soil_water_potential']) / abs(parameters['soil_depth']) + 1.0))

    soil_evaporation =  min(evaporative_demand, evaporative_supply)

    return soil_evaporation


def water_retention_curve(pF, theta=None, psi=None):
    """
    Water retention curve vanGenuchten - Mualem
    Args:
        pF - parameter dict
        theta - vol. water content (m3m-3)
        psi - matrix potential (m), <=0
    Returns:
        theta or psi
    """

    Ts = np.array(pF['theta_s'], ndmin=1)
    Tr = np.array(pF['theta_r'], ndmin=1)
    alfa = np.array(pF['alpha'], ndmin=1)
    n = np.array(pF['n'], ndmin=1)
    m = 1.0 - np.divide(1.0, n)

    def theta_psi(x):
        # converts water content (m3m-3) to potential (m)
        # checks limits
        x = np.minimum(x, Ts)
        x = np.maximum(x, Tr)
        s = (Ts - Tr) / ((x - Tr) + EPS)
        Psi = -1e-2 / alfa*(s**(1.0 / m) - 1.0)**(1.0 / n)  # m
        Psi[np.isnan(Psi)] = 0.0
        return Psi

    def psi_theta(x):
        # converts water potential (m) to water content (m3m-3)
        x = 100*np.minimum(x, 0)  # cm
        Th = Tr + (Ts - Tr) / (1 + abs(alfa*x)**n)**m
        return Th

    if theta is not None:
        return theta_psi(theta)
    else:
        return psi_theta(psi)

def hydraulic_conductivity(pF, wliq):
    r""" Unsaturated liquid-phase hydraulic conductivity following
    vanGenuchten-Mualem -model.

    Args:
        pF (dict):
            'theta_s' (float/array): saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'theta_r' (float/array): residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'alpha' (float/array): air entry suction [cm\ :sup:`-1`]
            'n' (float/array): pore size distribution [-]
            'pore_connectivity': pore connectivity parameter
        wliq (float or array): liquid water content
    Returns:
        Kh (float or array): hydraulic conductivity (if Ksat ~=1 then in [units], else relative [-])

    """
    w = np.array(wliq)

    # water retention parameters
    l = pF['pore_connectivity']
    n = np.array(pF['n'])
    m = 1.0 - np.divide(1.0, n)
    Ksat = pF['saturated_conductivity']

    # effectve saturation [-]
    S = np.minimum(1.0, (w - pF['theta_r']) / (pF['theta_s'] - pF['theta_r']) + EPS)

    Kh = Ksat * S**l * (1 - (1 - S**(1/m))**m)**2
    Kh = np.maximum(10*EPS, Kh)

    return Kh

def saturation_vapor_pressure(temperature):
    r""" Calculates saturation vapor pressure over water surface.

    Args:
        temperature: [\ :math:`^{\circ}`\ C]

    Returns:
        float: saturation vapor pressure in [Pa]
    """

    # [Pa]
    return 611.0 * np.exp((17.502 * temperature) / (temperature + 240.97))

#%%
""" alternative functions to estimate surface - air conductance for scalars
    not used in code at the moment!
"""

def moss_atm_conductance(wind_speed, roughness_height, dT=0.0, atten_factor=0.25):
    r""" Estimates boundary layer conductance of bryophyte canopy for paralell
    forced and free convection.

    Wind speed should represent vertical wind speed at ca. 20 cm above moss
    canopy (parametrsization derived from wind-tunnel studies). Roughness
    lengths scale is equal to the 'characteristic vertical height of
    individual bryophyte shoots' (typically order of 3-10 mm).

    Estimated CO\ :sub:`2` and heat conductances are related to water vapor
    by the ratio of molecular diffusivities
    (Campbell and Norman, 1998, eq. 7.29-7.33).

    References:
        Rice et al., 2001.
            Significance of variation in bryophyte canopy structure.
            Amer. J. Bot. 88:1568-1576.
        Rice, 2006.
            Towards an integrated undestanding of Bryophyte performance:
            the dimensions of space and time.
            Lindbergia 31:41-53.
        Schuepp, 1980.
            Observations on the use of analytical and numerical models for the
            description of transfer to porous surface vegetation such as
            lichen.
            Boundary-Layer Meteorol. 29: 59-73.
        Kondo & Ishida, 1997

    Args:
        wind_speed: [m s\ :sup:`-1`\ ]
        roughness_height: [m]
        deltaT: [degC], moss-air temperature difference
        atten_factor: [-] dimensionless attenuation factor for continuous moss carpets

    Returns:
        dictionary:
            boundary layer conductances for
                * 'co2': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'h2o': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'heat': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]

    """

    Sc_v = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_H2O
    Sc_c = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_CO2
    Pr = AIR_VISCOSITY / THERMAL_DIFFUSIVITY_AIR
    Re = wind_speed * roughness_height / AIR_VISCOSITY

    # Rice et al. (2001) eq. 1 ShSc**0.33 = CRe**n, where C=6.6e-4 and n=1.61.
    # C = 6.6e-4
    # n = 1.61

    #... however, exponent n for individual species is <1.53 so use median values of model
    # fitted to individual species by SL.
    C = 0.0067
    n = 1.27

    Sh_v = atten_factor * C*Re**n * Sc_v**0.33 # Sherwood numbner for H2O

    conductance_h2o = Sh_v * MOLECULAR_DIFFUSIVITY_H2O / roughness_height # ms-1

    # free convection as paralel pathway, based on Condo and Ishida, 1997.
    b = 2.2e-3 #ms-1K-1 b=1.1e-3 for smooth, 3.3e-3 for rough surface
    dT = np.maximum(dT, 0.0)
    gfree = Sc_v / Pr * b * dT**0.33  # mol m-2 s-1

    # [mol m-2 s-1]
    conductance_h2o = (conductance_h2o + gfree) * AIR_DENSITY


    # [mol m-2 s-1], differ from H2O by ratio of molecular diffusivities
    conductance_co2 = Sc_c / Sc_v * conductance_h2o

    conductance_heat = Pr / Sc_v * conductance_h2o

    return {
        'co2': conductance_co2,
        'h2o': conductance_h2o,
        'heat': conductance_heat
        }

def soil_boundary_layer_conductance(u, z, zo, Ta, dT, P=101300.):
    """
    Computes soil surface boundary layer conductance (mol m-2 s-1)
    assuming velocity profile logarithmic between z and z0.
    INPUT: u - mean velocity (m/s)
           z - height of mean velocity u (m)
           zo - soil surface roughness length for momentum (m)
           Ta - ambient temperature (degC)
           dT - soil surface-air temperature difference (degC)
           P - pressure(Pa)
    OUTPUT: boundary-layer conductances (mol m-2 s-1)
        gb_h - heat (mol m-2 s-1)
        gb_c- CO2 (mol m-2 s-1)
        gb_v - H2O (mol m-2 s-1)
    Based on Daamen & Simmons (1996). Note: gb decreases both in
    unstable and stable conditions compared to near-neutral;
    nonfeasible?
    Samuli Launiainen, 18.3.2014
    to python by Kersti
    """

    u = np.maximum(u, EPS)

    rho_air = 44.6*(P / 101300.0)*(273.15 / (Ta + 273.13))  # molar density of air [mol/m3]

    delta = 5.0 * GRAVITY * z * dT / ((Ta + 273.15) * u**2)
    if delta > 0:
        d = -0.75
    else:
        d = -2
    rb = (np.log(z/zo))**2 / (0.4**2*u)*(1 + delta)**d

    gb_h = rho_air * 1 / rb
    gb_v = MOLECULAR_DIFFUSIVITY_H2O / THERMAL_DIFFUSIVITY_AIR * gb_h
    gb_c = MOLECULAR_DIFFUSIVITY_CO2 / THERMAL_DIFFUSIVITY_AIR * gb_h

    return gb_h, gb_c, gb_v

    def water_heat_tendencies_bulk(self, y, dt, forcing, parameters):
        """
        Solves coupled water and heat balance over timestep dt.
        Returns temperature and water storage tendensies and
        water & energy fluxes.
        NOTE: solves layer bulk temperature directly!

        Args:
            y (array): initial values: returns derivatives
                0. temperature [degC] --> computes (du/dt)
                1. water storage [kg m-2] --> (du/dt)
                2. pond_recharge
                3. capillary_rise
                4. interception
                5. evaporation/condensation
                6. radiative_heat
                7. sensible_heat
                8. latent_heat
                9. conducted_heat
                10 advected_heat
            dt (float): time step
            forcing
            parameters
        Returns:
            derivatives (du/dt) of y.

        Units:  temperature [K s-1]
                water content and water fluxes [kg m-2 s-1 = mm s-1]
                energy fluxes [J m-2 s-1 = W m-2]
        """

        dudt = np.zeros(11)

        if dt == 0.0:
            dt = dt + EPS

        zm = 0.5 * self.height
        zs = abs(parameters['soil_depth'])
        zref = parameters['reference_height']

        max_storage = self.max_water_content * self.dry_mass
        min_storage = self.min_water_content * self.dry_mass

        # --- compute water balance ---
        # nomenclature:
        #   water_content [g g-1], water_storage [kg m-2 s-1 == mm]
        #   volumetric water content [-], water potential [m]
        #   all fluxes [kg m-2 s-1]

        # initial state

        # [kg m-2] or [mm]
        water_storage = min(max_storage, y[1])
        water_storage = max(min_storage, water_storage)

        # [g g-1]
        water_content = (water_storage / self.dry_mass)

        # [m m-3]
        volumetric_water = (water_content / WATER_DENSITY
                            * self.bulk_density)

        # [m]
#        water_potential = convert_hydraulic_parameters(
#                volumetric_water,
#                self.water_retention,
#                'volumetric_water')

        water_potential = water_retention_curve(self.water_retention,
                                                theta=volumetric_water)

        # --- constraints for recharge & evaporation
        max_recharge = max(max_storage - water_storage, 0.0)
        max_recharge = min(max_storage, max_recharge)

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max_recharge / dt

        # [kg m-2 s-1] or [mm s-1]
        max_evaporation_rate = (y[1] - (min_storage + EPS)) / dt

        if np.isinf(max_evaporation_rate) or max_evaporation_rate < 0.0:
            max_evaporation_rate = 0.0

        # [kg m-2 s-1] or [mm s-1]
        max_condensation_rate = -((max_storage - EPS) - y[1]) / dt

        if np.isinf(max_condensation_rate) or max_condensation_rate > 0.0:
            max_condensation_rate = 0.0

        # --- compute  evaporation/condensation rate [kg m-2 s-1] ---

        # boundary layer conductances for H2O, heat and CO2  [mol m-2 s-1]
        #dT = y[0] - forcing['air_temperature']

        conductance_to_air = surface_atm_conductance(wind_speed=forcing['wind_speed'],
                                                     zref=zref,
                                                     #friction_velocity=forcing['friction_velocity'],
                                                     zom=self.roughness_height,
                                                     dT=0.0)

        # Relative conductance is from Williams and Flanagan (1996), Oecologia.
        relative_conductance = min(1.0,
                                   (0.1285 * y[1]
                                    / self.dry_mass - 0.1285)) #* 0.5

        # [mol m-2 s-1]
        conductance_to_air_h2o = (
            conductance_to_air['h2o'] * relative_conductance)

        # [kg m-2 s-1]
        evaporation_demand = (conductance_to_air_h2o
                            * (611.0 / forcing['air_pressure']
                                * np.exp(17.502 * y[0] / (y[0] + 240.0))
                                - forcing['h2o'])
                            * MOLAR_MASS_H2O)

        evaporation_rate = min(evaporation_demand, max_evaporation_rate)
        evaporation_rate = max(evaporation_rate, max_condensation_rate)

        max_recharge_rate = max(max_recharge_rate + evaporation_rate, 0.0)
        max_recharge = max_recharge_rate * dt

        # -- recharge from rainfall interception and/or from pond storage
        # interception = min(max_recharge - EPS, forcing['precipitation'] * dt)
        # assumptotic function
        interception = (max_recharge *
                        (1.0 - np.exp(-(1.0 / max_storage)
                         * forcing['precipitation'] * dt)))

        # [kg m-2 s-1] or [mm s-1]
        interception_rate = interception / dt

        # [kg m-2] or [mm]
        max_recharge = max(max_recharge - interception, 0.0)

        pond_recharge = min(max_recharge - EPS, forcing['soil_pond_storage'])

        # [kg m-2 s-1] or [mm s-1]
        pond_recharge_rate = pond_recharge / dt

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max(max_recharge - pond_recharge, 0.0) / dt

        # --- compute capillary rise from soil [ kg m-2 s-1 = mm s-1]

        # Kh: soil-moss hydraulic conductivity assuming two resistors in series
        Km = hydraulic_conductivity(self.water_retention, volumetric_water)
        Ks = parameters['soil_hydraulic_conductivity']

        # conductance of layer [s-1]
        g_moss = Km / zm
        g_soil = Ks / zs

        # [m s-1]
        Kh = (g_moss * g_soil / (g_moss + g_soil)) * (zm + zs)

        capillary_rise = WATER_DENSITY * max(0.0, - Kh * ((water_potential - forcing['soil_water_potential'])
                                                          / (zm + zs) + 1.0))

        # [kg m-2 s-1] or [mm s-1]
        capillary_rise = min(capillary_rise, max_recharge_rate)

        # calculate mass balance of water

        # [kg m-2 s-1] or [mm s-1]
        dy_water = (
                interception_rate
                + pond_recharge_rate
                + capillary_rise
                - evaporation_rate
                )

        #--- compute energy balance
        # take reflectivities from previous timestep
        # radiation balance # [J m-2 s-1] or [W m-2]

        # [-]
        #albedo = self.reflectivity(water_content)
        #emissivity = self.optical_properties['emissivity']

        # [J m-2 s-1] or [W m-2]
        net_shortwave_radiation = (forcing['par'] * (1.0 - self.albedo['PAR'])
                                   + forcing['nir'] * (1.0 - self.albedo['NIR'])
                                  )


        net_longwave_radiation = self.emissivity * (forcing['lw_dn']
                                - STEFAN_BOLTZMANN * (y[0] + DEG_TO_KELVIN)**4.0)


        # [J m-2 s-1] or [W m-2]
        net_radiation = net_shortwave_radiation + net_longwave_radiation

        # [J m-2 s-1] or [W m-2]
        sensible_heat_flux = (SPECIFIC_HEAT_AIR
                              * conductance_to_air['heat']
                              * (y[0] - forcing['air_temperature']))

        # [J m-2 s-1] or [W m-2]
        latent_heat_flux = LATENT_HEAT / MOLAR_MASS_H2O * (evaporation_rate + EPS)

        # heat conduction between moss and soil
        # [W m-2 K-1]
        moss_thermal_conductivity = thermal_conductivity(volumetric_water)

        # thermal conductance [W m-2 K-1]
        # assume the layers act as two resistors in series
        g_moss = moss_thermal_conductivity / zm
        g_soil = parameters['soil_thermal_conductivity'] / zs

        thermal_conductance = (g_moss * g_soil) / (g_moss + g_soil)

        # [J m-2 s-1] or [W m-2]
        heat_conduction = thermal_conductance *(y[0] - forcing['soil_temperature'])

        # heat lost or gained with liquid water removing/entering

        # [J m-2 s-1] or [W m-2]
        heat_advection = SPECIFIC_HEAT_H2O * (
                        interception_rate * forcing['air_temperature']
                        + capillary_rise * forcing['soil_temperature']
                        + pond_recharge_rate * forcing['soil_temperature']
                        )

        # heat capacities
        # [J K-1]
        heat_capacity_old = (SPECIFIC_HEAT_ORGANIC_MATTER
                         * self.dry_mass
                         + SPECIFIC_HEAT_H2O * y[1])

        heat_capacity_new = (SPECIFIC_HEAT_ORGANIC_MATTER
                             * self.dry_mass
                             + SPECIFIC_HEAT_H2O * (y[1] + dy_water * dt))

        # calculate new temperature from heat balance

        heat_fluxes = (
                net_radiation
                - sensible_heat_flux
                - latent_heat_flux
                - heat_conduction
                + heat_advection
                )

        new_temperature = (heat_fluxes * dt + heat_capacity_old * y[0]) / heat_capacity_new

        # -- tendencies
        # [K m-2 s-1]
        dudt[0] = (new_temperature - y[0]) / dt
        # [kg m-2 s-1] or [mm s-1]
        dudt[1] = dy_water

        # water fluxes
        # [kg m-2 s-1] or [mm s-1]
        dudt[2] = pond_recharge_rate
        dudt[3] = capillary_rise
        dudt[4] = interception_rate
        dudt[5] = evaporation_rate

        # energy fluxes
        # [J m-2 s-1] or [W m-2]
        dudt[6] = net_radiation
        dudt[7] = sensible_heat_flux
        dudt[8] = latent_heat_flux
        dudt[9] = heat_conduction
        dudt[10] = heat_advection

        return dudt

#def bryophyte_shortwave_albedo(water_content, properties=None):
#    r""" Bryophyte albedo for PAR and NIR regions as function of water content
#
#    The effect of water content of spectral properties are based on studies by
#    Vogelmann and Moss (1993) and Fernandes (1999)on Sphagnum cuspidatum and
#    Pleurozium schreberi, respectively.
#
#    The albedo is scaled specific reflectance for PAR (400-700 nm) or
#    NIR (750-1400) regions. The scaling coefficient is common for both
#    PAR and NIR and it is based on relationship between normalized
#    reflectaces and hydration status. The species specific albedo is
#    assumed to represent a reflectance when bryophyte is in full hydration.
#
#    If bryophyte's properties are not given, estimation is based on generic
#    fit of water content against reflectances separately for PAR and NIR.
#    Fits are based on studies by Vogelmann and Moss (1993), and
#    Fernandes (1999) on Sphagnum cuspidatum and Pleurozium schreberi,
#    respectively.
#
#    References:
#        Vogelmann and Moss (1993)
#            Remote Sensing of Environment 45:273-279.
#        Fernandes (1999)
#            PhD thesis entitled: 'Scale influences of surface
#            parametrization on modelled boreal carbon and water budgets'
#
#    Args:
#        water_content (float): [g g\ :sup:`-2`]
#        max_water_content (float): [g g\ :sup:`-2`]
#
#    Returns (dict)
#    """
#
#    if properties is None:
#
#        def reflectance(x, a, b):
#            r""" Describes relationship between water content and reflectance
#
#            Args:
#                x (float): water_content
#                a (float): fitting parameter
#                b (float): fitting parameter
#
#            Returns:
#                percentage (float)
#            """
#            return np.abs(a * np.power(x, b))
#
#        # Original reflectance in percents
#        if isinstance(water_content, float):
#            if water_content < 0.3:
#                albedo_nir = 68.41/100.0
#            else:
#                albedo_nir = reflectance(water_content, 47.66, -0.3002) / 100.0
#            albedo_par = reflectance(water_content, 8.84, -0.1303) / 100.0
#        else:
#            albedo_nir = np.empty(np.shape(water_content))
#            albedo_nir = reflectance(water_content, 47.66, -0.3002) / 100.0
#            albedo_nir[water_content < 0.3] = 68.41 / 100.0
#
#            albedo_par = reflectance(water_content, 8.84, -0.1303) / 100.0
#
#        return {'PAR': albedo_par, 'NIR': albedo_nir}
#
#    else:
#        albedo_nir = properties['optical_properties']['albedo_NIR']
#        albedo_par = properties['optical_properties']['albedo_PAR']
#        normalized_water_content = water_content / properties['max_water_content']
#
#        scaling_coefficient = 1.0 + (4.5 - 1.0) / (1.00 + np.power(10, 4.53 * normalized_water_content))
#
#        albedo_par = scaling_coefficient * albedo_par
#        albedo_nir = scaling_coefficient * albedo_nir
#
#        return {'PAR': albedo_par, 'NIR': albedo_nir}

#def convert_hydraulic_parameters(value, water_retention, input_unit):
#    r""" Converts between volumetric water content and water potential.
#
#    Note:
#        In case of heterogenous poresystem, linear interpolation is used
#        to convert effective saturation to water potential. Water retention
#        curve with 500 points with logarithmic spacing in which water potential
#        is ranging from -1e-6 to -1e2. In case of effective saturation is out
#        of interpolation bounds, water potential is set to 0.0 m in lower bound
#        or -1e2 m in upper bound.
#
#    For conversion, VanGenuchten-model for water retention is used with
#    following parameters:
#            - saturated water content (theta_s, [cm\ :sup:`3` cm\ :sup:`-3`\ ])
#            - residual water content (theta_r, [cm\ :sup:`3` cm\ :sup:`-3`\ ])
#            - air entry suction (alpha, [cm\ :sup:`-1`\ ])
#            - pore size distribution (n, [-])
#
#    Volumetric water content is restricted to be between residual water
#    content and saturated water content.
#
#    Water potential is restricted to be between -10000 m and 0.0 m.
#
#    Args:
#        value:
#            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is VOLUMETRIC_WATER
#            * [m] if input is WATER_POTENTIAL
#        water_retention (dict): water retension parameters of BryoType object
#        input_unit:
#                * VOLUMETRIC_WATER
#                * WATER_POTENTIAL
#    Returns:
#        float:
#            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input was WATER_POTENTIAL
#            * [m] if input was VOLUMETRIC_WATER
#    """
#
#    if isinstance(input_unit, str):
#        input_unit.lower()
#
#    theta_s = water_retention['theta_s']  # sat. water content [m3m-3]
#    theta_r = water_retention['theta_r']  # residual water content [m3m-3]
#
#    if input_unit == 'water_potential':
#        eff_sat = effective_saturation(value,
#                                       water_retention,
#                                       input_unit)
#
#        return eff_sat * (theta_s - theta_r) + theta_r
#
#    elif input_unit == 'volumetric_water':
#
#        alpha = water_retention['alpha']  # micropore air entry suction [cm-1]
#        n = water_retention['n']  # micropore shape parameter [-]
#
#        eff_sat = effective_saturation(value, water_retention, input_unit)
#        inv_effs = 1.0 / eff_sat
#        m = np.divide(n, n - 1.0)
#
#        # [m]
#        psi = 1e-2 * 1.0 / -alpha * (inv_effs ** m - 1.0) ** (1.0 / n)
#
#        if isinstance(psi, list):
#            psi = np.array(list)
#
#        if isinstance(psi, np.ndarray):
#            psi[psi <= -1000.0] = -1000.0
#            #psi = np.where(psi <= -100.0, -100.0, psi)
#        else:
#            if psi <= -1000.0:
#                return -1000.0
#
#        return psi
#
#
#def hydraulic_conductivity(water_potential, water_retention, method=None):
#    r""" Estimates hydraulic conductivity in porous media.
#
#    Hydraulic conductivity is based on Voortman et al. (2014) and
#    Voortman et al. (2015) for xerophilious mosses.
#
#    The implementation is valid for unimodal and multimodal poresystems.
#    Multimodal Mualem-VanGenuchten is used as introduced by
#    Priesack and Durner (2006).
#
#    09.08.2018 AJK:
#    At the moment hydraulic conduction used only for micropore
#    system (capillary pore system) and it is only needed to calculate capillary
#    rise. Macropore system does not retain water and it is assumed that
#    macropore system does not restrict hydraulic conduction.
#
#    References:
#        Voortman et al. (2014) Hydrological Processes, 28, 6251-6264
#        Priesack and Durner (2006) Vadose Zone Journal, 5, 121-124
#
#    Lauren:
#    Estimation is based on Lauren (1999, table 4) and parametrization is valid
#    from -4 to -80 kPa. Assumption is that in the dry states
#    (lower than -80 kPa) there are no hydraulic conductivity.
#
#    Args:
#        water_potential: [m]
#        water_retention (list):
#            parameters of water retention curve
#                0. saturated water content (theta_s) [m\ :sup:`3` m :sup:`-3`\ ]
#                1. residual water content (theta_r) [m\ :sup:`3` m :sup:`-3`\ ]
#                2. air entry suction (alpha) [cm\ :sup:`-1`]
#                3. pore size distribution (n) [-]
#                4. saturated conductivity (K_s) [m s\ :sup:`-1`]
#                5. pore connectivity (l) [-]
#
#    Returns:
#        float: hydraulic conductivity in [m s\ :sup:`-1`\ ]
#    """
#
#    if isinstance(method, str):
#        method.lower()
#
#    saturated_conductivity = water_retention['saturated_conductivity']
#
#    if method == 'lauren':
#        # kPa (standard gravity 10 m/s2)
#        saturated_conductivity = water_retention['saturated_conductivity']
#        water_potential = -10.0 * water_potential
#
#        # Assuming that in dry states there is no hydraulic conductivity
#        if water_potential > 80.0:
#            return 0.0
#
#        if water_potential < 4.0:
#            water_potential = 4.0
#
#        # relative conductivity respect to 4 kPa
#        conductivity = (
#            np.power(10.0, 1.0/(-0.62 + 0.26 * np.log10(water_potential)))
#            / np.power(10.0, 1.0/(-0.62 + 0.26 * np.log10(4)))
#            )
#        # [m/s]
#        return conductivity * saturated_conductivity
#
#    else:
#        # Possibility to add more poresystems
#
#        #psi = 100.0 * np.minimum(water_potential, 0.0)
#        psi = np.minimum(water_potential, 0.0)
#
#        micropores = {'alpha': water_retention['alpha'],
#                      'n': water_retention['n']}
#
#        poresystems = [micropores]
#
#        coefficients = []
#        denominators = []
#        nominators = []
#
#        for idx in range(len(poresystems)):
#            m = 1.0 - np.divide(1.0, poresystems[idx]['n'])
#            alpha = poresystems[idx]['alpha']
#
#            sat_eff = effective_saturation(psi,
#                                           poresystems[idx],
#                                           'water_potential')
#
#            coefficients.append(sat_eff)
#
#            denom = 1.0 - np.power(sat_eff, 1.0 / m)
#            denom = 1.0 - np.power(denom, m)
#            denominators.append(alpha * denom)
#
#            nominators.append(alpha)
#
#        pore_connectivity = water_retention['pore_connectivity']
#
#        coefficient = np.power(np.sum(coefficients, axis=0), pore_connectivity)
#
#        denominator = np.power(np.sum(denominators, axis=0), 2.0)
#
#        nominator = np.power(np.sum(nominators, axis=0), 2.0)
#
#        return saturated_conductivity * coefficient * (denominator / nominator)

#def effective_saturation(value, water_retention, input_unit):
#    r""" Calculates effective saturation from volumetric water
#    content (theta, [m3 m-3]) or from water potential (psi, [m])
#
#    Args:
#        value (float):
#            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is VOLUMETRIC_WATER
#            * [m] if input is WATER_POTENTIAL
#        water_retention (dict):
#            if input_unit == 'water_potential'
#                'alpha'
#                'n'
#            if input_unit == VOLUMETRIC_WATER
#                'theta_r'
#                'theta_s'
#        input_unit (str): VOLUMETRIC_WATER or WATER_POTENTIAL
#    Returns:
#        float: [-]
#    """
#
#    options = set(['water_potential', 'volumetric_water'])
#
#    if isinstance(input_unit, str):
#        input_unit = input_unit.lower()
#    else:
#        raise ValueError("Input unit has to be string")
#
#    if input_unit not in options:
#        raise ValueError("Input unit options are: ".format(*options))
#
#
#    if input_unit == 'water_potential':
#        n = water_retention['n']
#        m = 1.0 - np.divide(1.0, n)
#        psi = 100.0 * np.minimum(value, 0.0)
#
#        eff_sat = (1.0 + (water_retention['alpha'] * abs(psi)) ** n) ** -m
#
#    elif input_unit == 'volumetric_water':
#
#        theta_r = water_retention['theta_r']
#        theta_s = water_retention['theta_s']
#
#        theta = np.minimum(value, theta_s)
#        theta = np.maximum(theta, theta_r)
#
#        eff_sat = np.divide(((theta - theta_r) + EPS), (theta_s - theta_r) + EPS)
#
#    return eff_sat
