# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:01:50 2017

@author: slauniai

******************************************************************************
CanopyModel:

Gridded canopy and snow hydrology model for SpatHy -integration
Uses simple schemes for computing water flows and storages within vegetation
canopy and snowpack at daily or sub-daily timesteps.

(C) Samuli Launiainen, 2016-2017
last edit: 1.11.2017: Added CO2-response to dry_canopy_et
******************************************************************************

"""
import numpy as np
eps = np.finfo(float).eps


class CanopyModel():
    def __init__(self, cpara):
        """
        initializes CanopyModel -object

        Args:
            cpara (dict): see canopy_parameters.py
        Returns:
            self (object):
                parameters/state variables:
                    Switch_MLM (boolean): control for multilayer computation
                    location (dict): 
                        'lat'(float): latitude
                        'lon'(float): longitude
                    LAI (float): total leaf area index [m2 m-2]
                    lad (float): total leaf area density [m2 m-3]
                    hc (float): canopy heigth [m]
                    cf (float): canopy closure [-]
                    gsref (float): LAI-weigthed average gsref of canopy [m s-1]
                    pheno_state (float): LAI-weigthed average phenological state of canopy (0..1)
                grid (only for multilayer computation):
                    z, dz, Nlayers (floats)
                objects:
                    Ptypes (.Pheno_Model, .LAI_model), Forestfloor (.Evap_Model, ...)
                    Radi_Model, Aero_Model, Interc_model, Canopy_Tr, Snow_Model
        """

        from radiation import Radiation
        from micromet import Aerodynamics
        from interception import Interception
        from evapotranspiration import Canopy_Transpiration
        from snow import Snowpack

        # --- control switches ---
        self.Switch_MLM = cpara['ctr']['multilayer_model']['ON']

        # --- site location ---
        self.location = cpara['loc']

        # --- grid ---
        if self.Switch_MLM:
            self.Nlayers = cpara['grid']['Nlayers']  # number of layers
            self.z = np.linspace(0, cpara['grid']['zmax'], self.Nlayers)  # grid [m] above ground
            self.dz = self.z[1] - self.z[0]  # gridsize [m]
        else:
            self.z = None

        # --- Plant types (with phenoligical models) ---
        ptypes = []
        stand_lai, stand_lad = 0.0, 0.0
        # para = [pinep, shrubp]
        for k in range(len(cpara['plant_types'])):
            p = cpara['plant_types'][k]
            for n in range(len(p['LAImax'])):
                pp = p.copy()
                pp['LAImax'] = p['LAImax'][n]
                if self.Switch_MLM:
                    pp['lad'] = p['lad'][:, n]
                else:
                    pp['lad'] = 0.0
                ptypes.append(PlantType(self.z, pp,
                                        Switch_pheno=cpara['ctr']['pheno_cylcle'],
                                        Switch_lai=cpara['ctr']['seasonal_LAI']))
                stand_lai += pp['LAImax']
                stand_lad += pp['lad']
        self.Ptypes = ptypes
        del p, pp, k, ptypes

        # --- stand characteristics ---
        self.LAI = stand_lai  # canopy total 1-sided leaf area index [m2m-2]
        self.lad = stand_lad  # canopy total 1-sided leaf area density [m2m-3]
        if self.Switch_MLM:
            f = np.where(self.lad > 0)[0][-1]
            self.hc = self.z[f].copy()  # canopy height [m]
        else:
            self.hc = cpara['canopy_para']['hc']  # canopy height [m]
            self.cf = cpara['canopy_para']['cf']  # canopy closure [-]

            # --- initialize canopy evapotranspiration computation ---
            self.Canopy_Tr = Canopy_Transpiration(cpara['phys_para'])

        # --- initialize radiation model ---
        self.Radi_Model = Radiation(cpara['radi'])

        # --- initialize aerodynamic resistances computation ---
        self.Aero_Model = Aerodynamics(cpara['aero'])

        # --- initialize interception model---
        self.Interc_Model = Interception(cpara['interc_snow'], self.LAI)

        # --- initialize snow model ---
        self.Snow_Model = Snowpack(cpara['interc_snow'], self.cf)

        # --- forest floor (models for evaporation...)---
        self.ForestFloor = ForestFloor({'soilrp': 300.0, 'f': cpara['phys_para']['f']})  ### --- soilrp to parameters?

    def _run_daily(self, doy, Ta, PsiL=0.0):
        """
        Computatations done at daily timestep. Updates 
        LAI and phenology
        Args:
            doy: day of year [days]
            Ta: mean daily air temperature [degC]
            PsiL: leaf water potential [MPa]
        Returns:
            None
        """

        """ update physiology and leaf area of planttypes & canopy"""
        stand_lai, stand_lad, gsref, pheno_state = 0.0, 0.0, 0.0, 0.0
        for pt in self.Ptypes:
            pt._update_daily(doy, Ta, PsiL=PsiL)  # updates pt properties
            stand_lad += pt.lad
            stand_lai += pt.LAI
            gsref += pt.LAI * pt.gsref
            pheno_state += pt.LAI * pt.pheno_state
        self.lad = stand_lad
        self.LAI = stand_lai
        self.gsref = gsref / stand_lai
        self.pheno_state = pheno_state / stand_lai

    def _run_timestep(self, dt, forcing, Rew=1.0, beta=1.0):
        """
        Runs CanopyModel instance for one timestep
        Args:
            dt: timestep [s]
            forcing (dataframe): meteorological forcing data
                'Prec': precipitation rate [m s-1]
                'Tair': air temperature [degC]
                'Rg': global radiation [W m-2]
                'vpd': vapor pressure deficit [kPa]
                'Par': fotosynthetically active radiation [W m-2]
                'U': wind speed [m s-1]
                'CO2': atm. CO2 mixing ratio [ppm]
                'P': pressure [Pa]
            Rew - relative extractable water [-], scalar or matrix
            beta - term for soil evaporation resistance (Wliq/FC) [-]
        Returns:
            fluxes (dict):
            state (dict):
        """

        Ta = forcing['Tair']
        Prec = forcing['Prec']
        Rg = forcing['Rg']
        Par = forcing['Par']
        VPD = forcing['vpd']
        U = forcing['U']
        P = forcing['P']
        CO2 = forcing['CO2']

        """ --- radiation --- """
        # estimate net radiation based on global radiation [W/m2]
        # Launiainen et al. 2016 GCB, fit to Fig 2a
        Rn = np.maximum(2.57 * self.LAI / (2.57 * self.LAI + 0.57) - 0.2, 0.55) * Rg
        # Rnet available at canopy and ground [W/m2]
        Rnet_c, Rnet_gr = self.Radi_Model.layerwise_Rnet(
                LAI=self.LAI,
                Rnet=Rn)

        """ --- aerodynamic conductances --- """
        ra, rb, ras, ustar, Uh, Ug = self.Aero_Model._run(
                LAI=self.LAI,
                hc=self.hc,
                Uo=U)

        """ --- interception and interception storage evaporation --- """
        Trfall_rain, Trfall_snow, Interc, Evap, MBE_interc = self.Interc_Model._run(
                dt=dt,
                LAI=self.LAI,
                cf=self.cf,
                T=Ta,
                Prec=Prec,
                AE=Rnet_c,
                VPD=VPD,
                Ra=ra,
                U=Uh)

        """--- snowpack ---"""
        PotInf, MBE_snow = self.Snow_Model._run(
                dt=dt,
                T=Ta,
                Trfall_rain=Trfall_rain,
                Trfall_snow=Trfall_snow)

        """--- canopy transpiration --- """  # DRY ??
        Transpi = self.Canopy_Tr._run(
                dt=dt,
                LAI=self.LAI,
                gsref=self.gsref,
                T=Ta,
                AE=Rnet_c,
                Qp=Par,
                VPD=VPD,
                Ra=ra,
                CO2=CO2,
                fPheno=self.pheno_state,
                Rew=Rew)

        """ --- forest floor (only in absence of snowpack) --- """
        if self.Snow_Model.swe <= eps:
            Efloor = self.ForestFloor._run(
                    dt=dt,
                    T=Ta,
                    AE=Rnet_gr,
                    VPD=VPD,
                    Ra=ras,
                    beta=beta)
        else:
            Efloor = 0.0

        # return state and fluxes in dictionary
        state = {'SWE': self.Snow_Model.swe,
                 'LAI': self.LAI,
                 'Phenof': self.pheno_state
                 }
        fluxes = {'PotInf': PotInf,
                  'Trfall': Trfall_rain + Trfall_snow,
                  'Interc': Interc,
                  'CanEvap': Evap,
                  'Transp': Transpi,
                  'Efloor': Efloor,
                  'MBE_interc': MBE_interc,
                  'MBE_snow': MBE_snow
                  }

        return fluxes, state

class PlantType():
    """
    PlantType -class.
    Contains plant-specific properties, state variables and phenology functions
    """

    def __init__(self, z, p, Switch_pheno=True, Switch_lai=True, Switch_WaterStress=True):
        """
        Creates PlantType
        Args:
            z - grid, evenly spaced, np.array
            p - parameters (dict)
            Switch_x - controls
        Returns:
            PlantType instance
        """

        from phenology import Photo_cycle, LAI_cycle

        self.Switch_pheno = Switch_pheno  # include phenology
        self.Switch_lai = Switch_lai  # seasonal LAI
        self.Switch_WaterStress = Switch_WaterStress  # water stress affects stomata

        self.z = z  # grid [m]

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
            self.relative_LAI = self.LAI_Model.f  # LAI relative to annual maximum [..1]
        else:
            self.relative_LAI = 1.0

        # physical structure
        self.LAImax = p['LAImax']  # maximum annual 1-sided LAI [m2m-2]
        self.LAI = self.LAImax * self.relative_LAI  # current LAI
        self.lad_normed = p['lad']  # normalized leaf-area density [m-1]
        self.lad = self.LAI * self.lad_normed  # current leaf-area density [m2m-3]

        if z != None:
            # plant height [m]
            f = np.where(self.lad_normed > 0)[0][-1]
            self.hc = z[f]
            # leaf gas-exchange parameters
            self.photop0 = p['photop']   # A-gs parameters at pheno_state = 1.0 (dict)
            self.photop = self.photop0.copy()  # current A-gs parameters (dict)
            # leaf properties
            self.leafp = p['leafp']  # leaf properties (dict)
            self.gsref = 0.0
        else:
            self.gsref = p['gsref']
 
    def _update_daily(self, doy, T, PsiL=0.0):
        """
        Updates PlantType pheno_state, gas-exchange parameters, LAI & lad
        Args:
            doy - day of year
            T - daily air temperature [degC]
            Psi_leaf - leaf (or soil) water potential, <0 [MPa]
        NOTE: CALL ONCE PER DAY
        """
        PsiL = np.minimum(-1e-5, PsiL)

        if self.Switch_pheno:
            self.pheno_state = self.Pheno_Model._run(T, out=True)

        if self.Switch_lai:
            self.relative_LAI =self.LAI_Model._run(doy, T, out=True)
            self.LAI = self.relative_LAI * self.LAImax
            self.lad = self.lad_normed * self.LAI
        """
        # scale photosynthetic capacity using vertical N gradient
        f = 1.0
        if 'kn' in self.photop0:
            kn = self.photop0['kn']
            Lc = np.flipud(np.cumsum(np.flipud(self.lad*self.dz)))
            Lc = Lc / Lc[0]
            f = np.exp(-kn*Lc)
            # print f
            # plt.plot(f, z, 'k', Lc, z, 'g-')
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
        """

class ForestFloor():
    """
    Forest floor
    """
    def __init__(self, p):

        from evapotranspiration import ForestFloor_Evap

        self.Evap = ForestFloor_Evap(p)

    def _run(self, dt, T, AE, VPD, Ra, beta):
        """
        Args:
            T: air temperature [degC]
            AE: available energy at forest floor [W m-2]
            VPD: vapor pressure deficit in [kPa]
            Ra: soil aerodynamic resistance [s m-1]
            beta: 
        Returns:
            Efloor: forest floor evaporation [m]
        """

        Efloor = self.Evap._run(dt, T, AE, VPD, Ra, beta)

        # Interception?

        return Efloor