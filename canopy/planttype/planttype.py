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
eps = np.finfo(float).eps  # machine epsilon
from .photo import leaf_interface
from .phenology import Photo_cycle, LAI_cycle
from .rootzone import RootUptake
from canopy.constants import LATENT_HEAT

class PlantType(object):
    r""" Contains plant-specific properties, state variables and phenological
    functions.
    """

    def __init__(self, z, p, dz_soil, ctr):
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
                'laip' (dict): parameters forleaf-area seasonal dynamics
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
                    'm':
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
                'StomaModel' (str): stomatal model e.g. 'MEDLYN_FARQUHAR'
                'WaterStress' (bool): account for water stress in planttypes --- TRUE NOT SUPPORTED!
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
            self.LAI_Model = LAI_cycle(p['laip'])  # LAI model instance
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
        self.StomaModel = ctr['StomaModel']

    def update_daily(self, doy, T, PsiL=0.0):
        r""" Updates planttype pheno_state, gas-exchange parameters, LAI and lad.

        Args:
            doy (float): day of year [days]
            Ta (float): mean daily air temperature [degC]
            PsiL (float): leaf water potential [MPa] --- CHECK??

        Note: Call once per day
        """
        PsiL = np.minimum(-1e-5, PsiL)

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
            Lc = Lc / np.maximum(Lc[0], eps)
            f = np.exp(-kn*Lc)
        # preserve proportionality of Jmax and Rd to Vcmax
        self.photop['Vcmax'] = f * self.pheno_state * self.photop0['Vcmax']
        self.photop['Jmax'] =  f * self.pheno_state * self.photop0['Jmax']
        self.photop['Rd'] =  f * self.pheno_state * self.photop0['Rd']

        if self.Switch_WaterStress:
            b = self.photop0['drp']
            if 'La' in self.photop0:
                # lambda increases with decreasing Psi as in Manzoni et al., 2011 Funct. Ecol.
                self.photop['La'] = self.photop0['La'] * np.exp(-b*PsiL)
            if 'm' in self.photop0:  # medlyn g1-model, decrease with decreasing Psi
                self.photop['m'] = self.photop0['m'] * np.maximum(0.05, np.exp(b*PsiL))

    def run(self, forcing, parameters, controls):
        r"""Computes dry leaf gas-exchange shale and sunlit leaves.

        Args:
            forcing (dict):
                'h2o'
                'co2'
                'air_temperature'
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
        # --- sunlit leaves
        sl = leaf_interface(
            self.photop,
            self.leafp,
            self.Tl_sl,
            forcing,
            controls,
            leaftype='sunlit',
            model=self.StomaModel
        )

        # --- shaded leaves
        sh = leaf_interface(
            self.photop,
            self.leafp,
            self.Tl_sh,
            forcing,
            controls,
            leaftype='shaded',
            model=self.StomaModel)

        if controls['energy_balance']:
            self.Tl_sh= sh['Tl'].copy()
            self.Tl_sl = sl['Tl'].copy()

        # integrate water and C fluxes over all leaves in PlantType, store resuts
        f_sl = forcing['radiation']['PAR']['sunlit']['fraction']
        pt_stats, layer_stats = self._integrate(sl, sh, f_sl, df)

        # --- sink/source terms
        sources = {'h2o': layer_stats['transpiration'],  # [mol m-3 s-1]
                   'co2': -layer_stats['net_co2'],  # [umol m-3 s-1]
                   'sensible_heat': layer_stats['sensible_heat'],  # [W m-3]
                   'latent_heat': layer_stats['transpiration'] * LATENT_HEAT,  # [W m-3]
                   'fr': layer_stats['fr']}  # [W m-3]

        return pt_stats, sources

    def _integrate(self, sl, sh, f_sl, df):
        """
        integrates layerwise statistics (per unit leaf area) to plant level
        Arg:
            sl, sh - dict of leaf-level outputs for sunlit and shaded leaves:

            x = {'An': An, 'Rd': Rd, 'E': fe, 'H': H, 'Fr': Fr, 'Tl': Tl, 'Ci': Ci,
                 'Cs': Cs, 'gs_v': gsv, 'gs_c': gs_opt, 'gb_v': gb_v}
        Returns:
            y - plant level statistics
        """
        # plant fluxes, weight factors is sunlit and shaded LAI at each layer
        f1 = df * f_sl * self.lad
        f2 = df * (1.0 - f_sl) * self.lad

        keys = ['net_co2', 'dark_respiration', 'transpiration', 'sensible_heat', 'fr']
        pt_stats = {k: np.sum((sl[k]*f1 + sh[k]*f2) * self.dz) for k in keys}  # flux per m-2(ground)

        layer_stats = {k: sl[k]*f1 + sh[k]*f2 for k in keys}  # flux per m-3

        # leaf temperatures
        pt_stats.update({'Tleaf': f_sl * sl['Tl'] + (1.0 - f_sl) * sh['Tl']})  # dry leaf temperature
        pt_stats.update({'Tleaf_sl': sl['Tl'] * self.lad})
        pt_stats.update({'Tleaf_sh': sh['Tl'] * self.lad})

        return pt_stats, layer_stats

