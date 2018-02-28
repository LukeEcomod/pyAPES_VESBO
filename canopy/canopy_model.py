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
import canopy_utils as cu
import phenology as pheno
eps = np.finfo(float).eps


class CanopyModel():
    def __init__(self, cpara):
        """
        initializes CanopyModel -object

        Args:
            cpara (dict): see canopy_parameters.py
        Returns:
            self (object):
                parameters:
                    LAI - leaf area index [m2 m-2]
                    hc - canopy heigth [m]
                    cf - closure [-]
                    physpara - physiologic parameters dict
                    phenopara - phenology param. dict
                    cpara - copy of inputs, for testing snow model
                    Wmax - maximum water storage capacity [m]
                    WmaxSnow -  maximum snow storage capacity [m]
                    Kmelt - degree day snowmelt factor [m degC-1 s-1]
                    Kfreeze - degree day freezing factor [m degC-1 s-1]
                    R - snow water retention coeff [-]
                state variables:
                    X - 'phenological state', i.e. delayed temperature [degC]
                    W - canopy water or snow storage [m]
                    SWE - snow water equivalent [m]
                    SWEi - snow water as ice [m]
                    SWEl - snow water as liquid [m]
                    ddsum - degree-day sum [degC]
        """

        # --- control switches ---
        self.Switch_MLM = cpara['ctr']['multilayer_model']

        # --- site location ---
        self.Lat = cpara['loc']['lat']
        self.Lon = cpara['loc']['lon']

        # --- grid ---
        if self.Switch_MLM:
            self.Nlayers = cpara['grid']['Nlayers']  # number of layers
            self.z = np.linspace(0, cpara['grid']['zmax'], self.Nlayers)  # grid [m] above ground
            self.dz = self.z[1] - self.z[0]  # gridsize [m]
        else:
            self.z = None

        # --- create PlantTypes ----
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

        self.physpara = cpara['phys_para']

        # interception and snow model parameters
        self.snowpara = cpara['interc_snow']
        self.snowpara.update({'kp': self.physpara['kp']})

        # --- for computing aerodynamic resistances
        self.zmeas = 2.0
        self.zground = 0.5
        self.zo_ground = 0.01
        self.soilrp = 300.0

        # --- state variables
        self.X = 0.0
        self.W = np.minimum(cpara['interc_snow']['w_ini'], cpara['interc_snow']['wmax'] * self.LAI)
        self.SWE = {'SWE': cpara['interc_snow']['swe_ini'],
                    'SWEi': cpara['interc_snow']['swe_ini'],
                    'SWEl': 0.0}

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
                'U': wind speed [m s-1] - if not given U = 2.0
                'CO2': atm. CO2 mixing ratio [ppm] - if not given C02 = 380
                'P': pressure [Pa] - if not given P = 101300.0
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
        if 'U' in forcing: 
            U = forcing['U']
        else:
            U = 2.0
        if 'P' in forcing:
            P = forcing['P']
        else:
            P = 101300.0
        if 'CO2' in forcing:
            CO2 = forcing['CO2']
        else:
            CO2 = 380.0

        # Rn = 0.7 * Rg #net radiation
        Rn = np.maximum(2.57 * self.LAI / (2.57 * self.LAI + 0.57) - 0.2,
                        0.55) * Rg  # Launiainen et al. 2016 GCB, fit to Fig 2a

        if self.SWE['SWE'] > 0.0:  # in case of snow, neglect forest floor evap
            f = 0.0
        else:
            f = self.physpara['f']

        """ --- aerodynamic conductances --- """
        Ra, Rb, Ras, ustar, Uh, Ug = cu.aerodynamics(self.LAI,
                                                     self.hc,
                                                     U,
                                                     w=0.01,
                                                     zm=self.zmeas,
                                                     zg=self.zground,
                                                     zos=self.zo_ground)

        """ --- interception, evaporation and snowpack --- """
        self.W, self.SWE, PotInf, Trfall, Evap, Interc, MBE, erate, unload, fact = \
            cu.canopy_water_snow(W=self.W,
                                 SWE=self.SWE,
                                 LAI=self.LAI,
                                 cf=self.cf,
                                 snowpara=self.snowpara,
                                 dt=dt,
                                 T=Ta,
                                 Prec=Prec,
                                 AE=Rn,
                                 VPD=VPD,
                                 Ra=Ra,
                                 U=U)

        """--- dry-canopy evapotranspiration [m s-1] --- """
        Transpi, Efloor, Gc = cu.dry_canopy_et(LAI=self.LAI,
                                               gsref=self.gsref,
                                               physpara=self.physpara,
                                               soilrp=self.soilrp,
                                               D=VPD,
                                               Qp=Par,
                                               AE=Rn,
                                               Ta=Ta,
                                               Ra=Ra,
                                               Ras=Ras,
                                               CO2=CO2,
                                               f=f,
                                               Rew=Rew,
                                               beta=beta,
                                               fPheno=self.pheno_state)

        Transpi = Transpi * dt
        Efloor = Efloor * dt
        ET = Transpi + Efloor

        # return state and fluxes in dictionary
        state = {"SWE": self.SWE['SWE'],
                 "LAI": self.LAI,
                 "Phenof": self.pheno_state
                 }
        fluxes = {'PotInf': PotInf,
                  'Trfall': Trfall,
                  'Interc': Interc,
                  'CanEvap': Evap,
                  'Transp': Transpi, 
                  'Efloor': Efloor,
                  'MBE': MBE
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
        self.Switch_pheno = Switch_pheno  # include phenology
        self.Switch_lai = Switch_lai  # seasonal LAI
        self.Switch_WaterStress = Switch_WaterStress  # water stress affects stomata

        self.z = z  # grid [m]

        # phenology model
        if self.Switch_pheno:
            self.Pheno_Model = pheno.Photo_cycle(p['phenop'])  # phenology model instance
            self.pheno_state = self.Pheno_Model.f  # phenology state [0...1]
        else:
            self.pheno_state = 1.0

        # dynamic LAI model
        if self.Switch_lai:
            # seasonality of leaf area
            self.LAI_Model = pheno.LAI_cycle(p['laip'])  # LAI model instance
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
            self.relative_LAI =self.LAI_Model._run(doy, T,  out=True)
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