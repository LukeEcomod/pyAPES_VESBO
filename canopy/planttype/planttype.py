# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
.. module: planttype
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes plant seasonal cycle and dry leaf gas exchange.
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Oct 02 09:04:05 2018

Note:
    migrated to python3
    - absolute imports

References:
Launiainen, S., Katul, G.G., Lauren, A. and Kolari, P., 2015. Coupling boreal
forest CO2, H2O and energy flows by a vertically structured forest canopy –
Soil model with separate bryophyte layer. Ecological modelling, 312, pp.385-405.

"""

import numpy as np
from .photo import leaf_interface
from .phenology import Photo_cycle, LAI_cycle
from .rootzone import RootUptake
from canopy.constants import PAR_TO_UMOL, EPS

class PlantType(object):
    r""" Contains plant-specific properties, state variables and phenological
    functions.
    """

    def __init__(self, z, p, dz_soil, ctr, loc):
        r""" Initialises a planttype object and submodel objects
        using given parameters.

        Args:
            z (array): canopy model nodes, height from soil surface (= 0.0) [m]
            p (dict):
                'name' (str): name of planttype
                'LAImax' (float): leaf area index
                'lad' (array): normalized leaf area density profile
                'phenop' (dict): parameters for seasonal cycle of phenology
                    'Xo': initial delayed temperature [degC]
                    'fmin': minimum photocapacity [-]
                    'Tbase': base temperature [degC]
                    'tau': time constant [days]
                    'smax': threshold for full acclimation [degC]
                'laip' (dict): parameters for LAI seasonal dynamics
                    'lai_min': minimum LAI, fraction of annual maximum [-]
                    'lai_ini': initial LAI fraction, if None lai_ini = Lai_min * LAImax
                    'DDsum0': degreedays at initial time [days]
                    'Tbase': base temperature [degC]
                    'ddo': degreedays at bud burst [days]
                    'ddur': duration of recovery period [days]
                    'sso': start doy of decrease, based on daylength [days]
                    'sdur': duration of decreasing period [days]
                'photop' (dict): leaf gas-exchange parameters
                    'Vcmax': maximum carboxylation velocity [umolm-2s-1]
                    'Jmax': maximum rate of electron transport [umolm-2s-1]
                    'Rd': dark respiration rate [umolm-2s-1]
                    'alpha': quantum yield parameter [mol/mol]
                    'theta': co-limitation parameter of Farquhar-model
                    'La': stomatal parameter (Lambda, m, ...) depending on model
                    'g1':
                    'g0': residual conductance for CO2 [molm-2s-1]
                    'kn':
                    'beta':  co-limitation parameter of Farquhar-model
                    'drp':
                    'tresp' (dict): temperature sensitivity parameters
                        'Vcmax': [Ha, Hd, Topt]; activation energy [kJmol-1], deactivation energy [kJmol-1], optimum temperature [degC]
                        'Jmax': [Ha, Hd, Topt];
                        'Rd': [Ha]; activation energy [kJmol-1)]
                'leafp' (dict): leaf properties
                    'lt': leaf lengthscale [m]
                    'par_alb': leaf Par albedo [-]
                    'nir_alb': leaf Nir albedo [-]
                    'emi': leaf emissivity [-]
                'rootp' (dict): root zone properties
                    'root_depth': root depth [m]
                    'beta': shape parameter for root distribution model
                    'RAI_LAI_multiplier': multiplier for total fine root area index (RAI = 2*LAImax)
                    'fine_radius': fine root radius [m]
                    'radial_K': maximum bulk root membrane conductance in radial direction [s-1]
            dz_soil (array): thickness of soilprofile layers, needed for rootzone [m]
            ctr (dict): switches and specifications for computation
                'WaterStress' (str): account for water stress using 'Rew', 'PsiL' or 'None'
                'seasonal_LAI' (bool): account for seasonal LAI dynamics
                'pheno_cycle' (bool): account for phenological cycle
        Returns:
            self (object):
                .name (str)
                .pheno_state(float): phenology state [0...1]
                .relative_LAI (float): LAI relative to annual maximum [0...1]
                .lad (array): leaf area density [m2 m-3]
                .lad_normed (array): normalized leaf area density [-]
                ...
                .Switch_pheno (bool): account for phenological cycle
                .Switch_lai (bool): account for seasonal LAI dynamics
                .Switch_WaterStress (bool): account for water stress in planttypes
                .Pheno_Model (object): model for phenological cycle
                .LAI_Model (object): model for seasonal development of LAI
                .Roots (object): root properties
        """

        self.Switch_pheno = ctr['pheno_cycle']  # include phenology
        self.Switch_lai = ctr['seasonal_LAI']  # seasonal LAI
        self.Switch_WaterStress = ctr['WaterStress']  # water stress affects stomata
        self.StomaModel = 'MEDLYN_FARQUHAR' # stomatal model

        self.name = p['name']

        # phenology model
        if self.Switch_pheno:
            self.Pheno_Model = Photo_cycle(p['phenop'])  # phenology model instance
            self.pheno_state = self.Pheno_Model.f  # phenology state [0...1]
        else:
            self.pheno_state = 1.0

        # dynamic LAI model
        if self.Switch_lai:
            # seasonality of leaf area
            self.LAI_Model = LAI_cycle(p['laip'], loc)  # LAI model instance
            self.relative_LAI = self.LAI_Model.f  # LAI relative to annual maximum [0...1]
        else:
            self.relative_LAI = 1.0

        # physical structure
        self.LAImax = p['LAImax']  # maximum annual 1-sided LAI [m2m-2]
        self.LAI = self.LAImax * self.relative_LAI  # current LAI
        self.lad_normed = p['lad']  # normalized leaf-area density [m-1]
        self.lad = self.LAI * self.lad_normed  # current leaf-area density [m2m-3]

        # root properties
        self.Roots = RootUptake(p['rootp'], dz_soil, self.LAImax)

        self.dz = z[1] - z[0]

        # leaf gas-exchange parameters
        self.photop0 = p['photop']   # A-gs parameters at pheno_state = 1.0 (dict)
        self.photop = self.photop0.copy()  # current A-gs parameters (dict)
        # leaf properties
        self.leafp = p['leafp']  # leaf properties (dict)

    def update_daily(self, doy, T, PsiL=0.0, Rew=1.0):
        r""" Updates planttype pheno_state, gas-exchange parameters, LAI and lad.

        Args:
            doy (float): day of year [days]
            Ta (float): mean daily air temperature [degC]
            PsiL (float): leaf water potential [MPa] --- CHECK??
            Rew (float): relatively extractable water (-)

        Note: Call once per day
        """

        if self.Switch_pheno:
            self.pheno_state = self.Pheno_Model.run(T, out=True)

        if self.Switch_lai:
            self.relative_LAI =self.LAI_Model.run(doy, T, out=True)
            self.LAI = self.relative_LAI * self.LAImax
            self.lad = self.lad_normed * self.LAI

        # scale photosynthetic capacity using vertical N gradient
        f = 1.0
        if 'kn' in self.photop0:
            kn = self.photop0['kn']
            Lc = np.flipud(np.cumsum(np.flipud(self.lad*self.dz)))
            Lc = Lc / np.maximum(Lc[0], EPS)
            f = np.exp(-kn*Lc)
        # preserve proportionality of Jmax and Rd to Vcmax
        self.photop['Vcmax'] = f * self.pheno_state * self.photop0['Vcmax']
        self.photop['Jmax'] =  f * self.pheno_state * self.photop0['Jmax']
        self.photop['Rd'] =  f * self.pheno_state * self.photop0['Rd']

# TEST
        self.photop['Vcmax'] = f * self.pheno_state * self.relative_LAI * self.photop0['Vcmax']
        self.photop['Jmax'] =  f * self.pheno_state * self.relative_LAI * self.photop0['Jmax']
        self.photop['Rd'] =  f * self.pheno_state * self.relative_LAI * self.photop0['Rd']

        # water stress responses
        if self.Switch_WaterStress == 'Rew':
            # drought responses from Hyde pine shoot chambers, 2006
            # for 'Medlyn - model'
            b = self.photop['drp']
            fm = np.minimum(1.0, (Rew / b[0])**b[1])
            self.photop['g1'] = fm * self.photop0['g1']

            # apparent Vcmax decrease with Rew
            fv = np.minimum(1.0, (Rew / b[2])**b[3])
            self.photop['Vcmax'] *= fv
            self.photop['Jmax'] *= fv
            self.photop['Rd'] *= fv

        if self.Switch_WaterStress == 'PsiL':
#            PsiL = np.maximum(np.minimum(-1e-5, PsiL),-3)
            PsiL = np.minimum(-1e-5, PsiL)
            b = self.photop0['drp']
            # medlyn g1-model, decrease with decreasing Psi
            self.photop['g1'] = self.photop0['g1'] * np.maximum(0.05, np.exp(b*PsiL))

            # Vmax and Jmax responses to leaf water potential. Kellomäki & Wang, 1996.
            # (Huom! artikkelissa kertoimet väärinpäin, tarkistettu kuvista)
            fv = 1.0 / (1.0 + (PsiL / - 2.04)**2.78)  # vcmax
            fj = 1.0 / (1.0 + (PsiL / - 1.56)**3.94)  # jmax
            fr = 1.0 / (1.0 + (PsiL / - 2.53)**6.07)  # rd
            self.photop['Vcmax'] *= fv
            self.photop['Jmax'] *= fj
            self.photop['Rd'] *= fr

    def run(self, forcing, parameters, controls):
        r"""Computes dry leaf gas-exchange for shaded and sunlit leaves.

        Args:
            forcing (dict):
                'h2o'
                'co2'
                'air_temperature'
                'air_pressure'
                'wind_speed'
                'nir'
                'par'
                'nir'
            'parameters' (dict):
                'sunlit_fraction'
                'dry_leaf_fraction'
            'controls' (dict):
                'energy_balance'

            f_sl (array): sunlit fraction [-]
            H2O (array): water vapor mixing ratio [mol mol-1]
            CO2 (array): carbon dioxide mixing ratio [ppm]
            T (array): ambient air temperature [degC]
            U (array): mean wind speed [m s-1]
            P: ambient pressure [Pa]
            Q_sl1, Q_sh1 (arrays): incident PAR at sunlit and shaded leaves [Wm-2]
            SWabs_sl, SWabs_sh (arrays): absorbed SW (PAR + NIR) at sunlit and shaded leaves [Wm-2]
            LWl (array): leaf net long-wave radiation [Wm-2]
            df (array): dry leaf fraction [-]
            Ebal (bool): solve dry leaf energy balance
            Tl_ave (array): average leaf temperature used in LW computation [degC]
            gr (array): radiative conductance [mol m-2 s-1]
        """
        df = parameters['dry_leaf_fraction']

        controls.update({
            'photo_model': self.StomaModel
        })

        leaf_forcing = {
            'h2o': forcing['h2o'],
            'co2': forcing['co2'],
            'air_temperature': forcing['air_temperature'],
            'wind_speed': forcing['wind_speed'],
            'air_pressure': forcing['air_pressure'],
            'average_leaf_temperature': forcing['leaf_temperature']
        }

        # --- compute sunlit leaves
        if controls['energy_balance']:
            sw_absorbed_sl = (forcing['par']['sunlit']['absorbed']
                        + forcing['nir']['sunlit']['absorbed'])

            leaf_forcing.update({
                'radiative_conductance': forcing['lw']['radiative_conductance'],
                'lw_net': df * forcing['lw']['net_leaf'],
                'sw_absorbed': df * sw_absorbed_sl,
                'leaf_temperature': self.Tl_sl
            })

        logger_info = controls['logger_info']
        logger_info = logger_info + ' leaftype: sunlit'


        Qp_sl = forcing['par']['sunlit']['incident'] * PAR_TO_UMOL
        leaf_forcing.update({'par_incident': Qp_sl}) # for sunlit

        sl = leaf_interface(
            photop=self.photop,
            leafp=self.leafp,
            forcing=leaf_forcing,
            controls=controls,
            df=df,
            logger_info=logger_info)

        # --- compute shaded leaves
        logger_info = controls['logger_info']
        logger_info = logger_info + ' leaftype: shaded'

        if controls['energy_balance']:
            sw_absorbed_sh = (forcing['par']['shaded']['absorbed']
                        + forcing['nir']['shaded']['absorbed'])

            leaf_forcing.update({'sw_absorbed': df * sw_absorbed_sh,  # for shaded
                                 'leaf_temperature': self.Tl_sh})  # for shaded

        Qp_sh = forcing['par']['shaded']['incident'] * PAR_TO_UMOL
        leaf_forcing.update({'par_incident': Qp_sh}) # for sunlit

        sh = leaf_interface(
            photop=self.photop,
            leafp=self.leafp,
            forcing=leaf_forcing,
            controls=controls,
            df=df,
            logger_info=logger_info)

        if controls['energy_balance']:
            self.Tl_sh = sh['Tl'].copy()
            self.Tl_sl = sl['Tl'].copy()

        # integrate water and C fluxes over all leaves in PlantType, store resuts
        f_sl = parameters['sunlit_fraction']
        pt_stats, layer_stats = self._integrate(sl, sh, f_sl)

        # --- sink/source terms
        sources = {'h2o': layer_stats['transpiration'],  # [mol m-3 s-1]
                   'co2': -layer_stats['net_co2'],  # [umol m-3 s-1]
                   'sensible_heat': layer_stats['sensible_heat'],  # [W m-3]
                   'latent_heat': layer_stats['latent_heat'],  # [W m-3]
                   'fr': layer_stats['fr']}  # [W m-3]

        return pt_stats, sources

    def _integrate(self, sl, sh, f_sl):
        """
        integrates layerwise statistics (per unit leaf area) to plant level
        Arg:
            sl, sh - dict of leaf-level outputs for sunlit and shaded leaves:

            x = {'An': An, 'Rd': Rd, 'E': fe, 'H': H, 'Fr': Fr, 'Tl': Tl, 'Ci': Ci,
                 'Cs': Cs, 'gs_v': gsv, 'gs_c': gs_opt, 'gb_v': gb_v}
        Returns:
            y - plant level statistics
        """

        # plant fluxes, weight factors are sunlit and shaded LAI at each layer
        f1 = f_sl * self.lad
        f2 = (1.0 - f_sl) * self.lad

        keys = ['net_co2', 'dark_respiration', 'transpiration', 'latent_heat', 'sensible_heat', 'fr']

        pt_stats = {k: np.sum((sl[k]*f1 + sh[k]*f2) * self.dz) for k in keys}  # flux per m-2(ground)

        layer_stats = {k: sl[k]*f1 + sh[k]*f2 for k in keys}  # flux per m-3

        keys = ['stomatal_conductance', 'boundary_conductance','leaf_internal_co2', 'leaf_surface_co2']
        pt_stats.update({k: np.sum((sl[k]*f1 + sh[k]*f2)) / np.sum(self.lad + EPS) for k in keys})

        # leaf temperatures
        pt_stats.update({'Tleaf': f_sl * sl['Tl'] + (1.0 - f_sl) * sh['Tl']})  # dry leaf temperature
        pt_stats.update({'Tleaf_sl': sl['Tl'] * self.lad})
        pt_stats.update({'Tleaf_sh': sh['Tl'] * self.lad})

        return pt_stats, layer_stats

