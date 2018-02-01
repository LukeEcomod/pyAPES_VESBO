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
eps = np.finfo(float).eps


class CanopyModel():
    def __init__(self, cpara, lai_conif=None, lai_decid=None, cf=None, hc=None,  cmask=np.ones(1), outputs=False):
        """
        initializes CanopyModel -object

        Args:
            cpara - parameter dict:
                {'spatial': no, 'wmaxsnow': 1.0, 'lon': 24.21, 'cf': 0.7,
                'wmax': 0.5, 'albedo': 0.1, 'swe': 40.0,pmax': 20.0,
                'q50': 30.0, 'lai': 4.0, 'kfreeze': 5.79e-06, 'elev': 100.0,
                'lat': 60.38, 'hc': 22.0, 'dt': 86400.0, 'kmelt': 2.8934e-05,
                'ka': 0.6, 'f': 0.4, 'emi': 0.98, 'kp': 0.6, 'gsref': 0.0018,
                'w': 0.0, 'clump': 0.7}
            lai - (conifer) leaf area index grid; if spatial: yes
            lai_decid - (deciduous) leaf area index grid;            
            cf - canopy closure grid
            hc - canopy height grid
            cmask - catchment mask (when lai=None)
            outputs - True saves output grids to list at each timestep

        Returns:
            self - object
        self.
            LAI - leaf area index [m2m-2]
            hc - canopy heigth [m]
            cf - closure [-]
            physpara - physiologic parameters dict
            phenopara - phenology param. dict
            cpara - copy of inputs, for testing snow model
            Wmax - maximum water storage capacity [mm = kgm-2(ground)]
            WmaxSnow -  maximum snow storage capacity [mm]
            Kmelt - degree day snowmelt factor [mm d-1 K-1]
            Kfreeze - degree day freezing factor [mm d-1 K-1]
            R - snow water retention coeff [-]

            State variables:
            X - 'phenological state', i.e. delayed temperature [degC]
            W - canopy water or snow storage [mm]
            SWE - snow water equivalent [mm]
            SWEi - snow water as ice [mm]
            SWEl - snow water as liquid [mm]
            ddsum - degree-day sum [degC]
        """
        epsi = 0.01

        self.Lat = cpara['lat']
        self.Lon = cpara['lon']

        # physiology: transpi + floor evap
        self.physpara = {'q50': cpara['q50'],
                         'gsref_conif': cpara['gsref_conif'],
                         'gsref_decid': cpara['gsref_decid'],
                         'kp': cpara['kp'],
                         'f': cpara['f'], 'rw': cpara['rw'],
                         'rwmin': cpara['rwmin'], 'ga': cpara['ga']}

        # phenology
        self.phenopara = {'Smax': cpara['smax'], 'Xo': cpara['xo'],
                          'tau': cpara['tau'], 'fmin': cpara['fmin'],
                          'ddo': cpara['ddo'], 'ddur': cpara['ddur'],
                          'sdl': cpara['sdl'], 'sdur': cpara['sdur'],
                          'lai_decid_min': cpara['lai_min']}

        if lai_conif is not None:
            # print '**** CanopyModel - stands from lai data *******'
            self._LAIconif = lai_conif + epsi  # m2m-2

            if lai_decid is not None:
                self._LAIdecid_max = lai_decid + epsi
            else:
                self._LAIdecid_max = np.zeros(np.shape(lai_conif)) + epsi

            self.hc = hc + epsi
            self.cf = cf + epsi
            # self.cf = 0.1939 * ba / (0.1939 * ba + 1.69) + epsi
            # canopy closure [-] as function of basal area ba m2ha-1;
            # fitted to Korhonen et al. 2007 Silva Fennica Fig.2

        else:  # spatially constant properties are used, given in cpara
            # print '**** CanopyModel - stands with constant values *******'
            self._LAIconif = cpara['lai_conif'] * cmask
            self._LAIdecid_max = cpara['lai_decid'] * cmask
            self.hc = cpara['hc'] * cmask
            self.cf = cpara['cf'] * cmask

        # current deciduous and total LAI
        self._LAIdecid = self.phenopara['lai_decid_min'] * self._LAIdecid_max
        self.LAI = self._LAIconif + self._LAIdecid

        gridshape = np.shape(self.LAI)
        # self.cf = np.ones(gridshape)*0.9
        # deciduous leaf growth stage
        self._growth_stage = 0.0
        self._senesc_stage = 0.0

        # senescence starts at first doy when daylength < self.phenopara['sdl']
        doy = np.arange(1, 366)
        dl = cu.daylength(self.Lat, self.Lon, doy)

        ix = np.max(np.where(dl > self.phenopara['sdl']))
        self.phenopara['sso'] = doy[ix]  # this is onset date for senescence
        del ix

        # self.cpara = cpara  # added new parameters self.cpara['kmt'],
        # self.cpara['kmr'] here for testing radiation-based snow melt model
        self.wmax = cpara['wmax']
        self.wmaxsnow = cpara['wmaxsnow']
        self.Kmelt = cpara['kmelt']
        self.Kfreeze = cpara['kfreeze']
        self.R = 0.05  # max fraction of liquid water in snow

        # --- for computing aerodynamic resistances
        self.zmeas = 2.0
        self.zground = 0.5
        self.zo_ground = 0.01
        self.soilrp = 300.0

        # --- state variables
        self.X = 0.0
        self.W = np.minimum(cpara['w'], self.wmax*self.LAI)
        self.SWE = cpara['swe']*np.zeros(gridshape)
        self.SWEi = self.SWE
        self.SWEl = np.zeros(gridshape)
        self.DDsum = 0.0

        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'PotInf': [], 'Trfall': [], 'Interc': [], 'Evap': [], 'ET': [],
                            'Transpi': [], 'Efloor': [], 'SWE': [], 'LAI': [],
                            'LAIfract': [], 'Mbe': [], 'LAIdecid': [], 'erate': [], 'Unload': [],
                            'fact': []}

    def _run(self, doy, dt, Ta, Prec, Rg, Par, VPD, U=2.0, CO2=380.0, Rew=1.0, beta=1.0, P=101300.0):
        """
        Runs CanopyModel instance for one timestep
        IN:
            doy - day of year
            dt - timestep [s]
            Ta - air temperature  [degC], scalar or (n x m) -matrix
            prec - precipitatation rate [mm/s]
            Rg - global radiation [Wm-2], scalar or matrix
            Par - photos. act. radiation [Wm-2], scalar or matrix
            VPD - vapor pressure deficit [kPa], scalar or matrix
            U - mean wind speed at ref. height above canopy top [ms-1], scalar or matrix
            CO2 - atm. CO2 mixing ratio [ppm]
            Rew - relative extractable water [-], scalar or matrix
            beta - term for soil evaporation resistance (Wliq/FC) [-]
            P - pressure [Pa], scalar or matrix
        OUT:
            updated CanopyModel instance state variables
            flux grids PotInf, Trfall, Interc, Evap, ET, MBE [mm]
        """

        # Rn = 0.7 * Rg #net radiation
        Rn = np.maximum(2.57 * self.LAI / (2.57 * self.LAI + 0.57) - 0.2,
                        0.55) * Rg  # Launiainen et al. 2016 GCB, fit to Fig 2a


        f = self.physpara['f'] * np.ones(np.shape(self.SWE))
        f[self.SWE > 0] = eps  # in case of snow, neglect forest floor evap

        """ --- update phenology: self.ddsum & self.X ---"""
        self._degreeDays(Ta, doy)
        fPheno = self._photoacclim(Ta)

        """ --- update deciduous leaf area index --- """
        laifract = self._lai_dynamics(doy)

        """ --- aerodynamic conductances --- """
        Ra, Rb, Ras, ustar, Uh, Ug = cu.aerodynamics(self.LAI, self.hc, U, w=0.01, zm=self.zmeas,
                                                  zg=self.zground, zos=self.zo_ground)

        """ --- interception, evaporation and snowpack --- """
        PotInf, Trfall, Evap, Interc, MBE, erate, unload, fact = self.canopy_water_snow(dt, Ta, Prec, Rn, VPD, Ra=Ra)

        """--- dry-canopy evapotranspiration [mm s-1] --- """
        Transpi, Efloor, Gc = self.dry_canopy_et(VPD, Par, Rn, Ta, Ra=Ra, Ras=Ras, CO2=CO2, f=f,
                                                 Rew=Rew, beta=beta, fPheno=fPheno)

        Transpi = Transpi * dt
        Efloor = Efloor * dt
        ET = Transpi + Efloor

        # append results to lists; use only for testing small grids!
        if hasattr(self, 'results'):
            self.results['PotInf'].append(PotInf)
            self.results['Trfall'].append(Trfall)
            self.results['Interc'].append(Interc)
            self.results['Evap'].append(Evap)
            self.results['ET'].append(ET)
            self.results['Transpi'].append(Transpi)
            self.results['Efloor'].append(Efloor)
            self.results['SWE'].append(self.SWE)
            self.results['LAI'].append(self.LAI)
            self.results['LAIfract'].append(laifract)
            self.results['Mbe'].append(np.nanmax(MBE))
            self.results['LAIdecid'].append(self._LAIdecid)
            self.results['erate'].append(erate)
            self.results['Unload'].append(unload)
            self.results['fact'].append(fact)

        # return state and fluxes in dictionary
        state = {"SWE": self.SWE,
                 "LAI": self.LAI,
                 "LAIdecid": self._LAIdecid,
                 }
        fluxes = {'PotInf': PotInf,
                  'Trfall': Trfall,
                  'Interc': Interc,
                  'CanEvap': Evap,
                  'Transp': Transpi, 
                  'Efloor': Efloor,
                  'MBE': np.nanmax(MBE)
                  }

        return fluxes, state

    def _degreeDays(self, T, doy):
        """
        Calculates and updates degree-day sum from the current mean Tair.
        INPUT:
            T - daily mean temperature (degC)
            doy - day of year 1...366 (integer)
        """
        To = 5.0  # threshold temperature
        if doy == 1:  # reset in the beginning of the year
            self.DDsum = 0.
        else:
            self.DDsum += np.maximum(0.0, T - To)

    def _photoacclim(self, T):
        """
        computes new stage of temperature acclimation and phenology modifier.
        Peltoniemi et al. 2015 Bor.Env.Res.
        IN: object, T = daily mean air temperature
        OUT: fPheno - phenology modifier [0...1], updates object state
        """

        self.X = self.X + 1.0 / self.phenopara['tau'] * (T - self.X)  # degC
        S = np.maximum(self.X - self.phenopara['Xo'], 0.0)
        fPheno = np.maximum(self.phenopara['fmin'],
                            np.minimum(S / self.phenopara['Smax'], 1.0))
        return fPheno

    def _lai_dynamics(self, doy):
        """
        Seasonal cycle of deciduous leaf area

        Args:
            self - object
            doy - day of year

        Returns:
            none, updates state variables self.LAIdecid, self._growth_stage,
            self._senec_stage
        """
        lai_min = self.phenopara['lai_decid_min']
        ddo = self.phenopara['ddo']
        ddur = self.phenopara['ddur']
        sso = self.phenopara['sso']
        sdur = self.phenopara['sdur']

        # growth phase
        if self.DDsum <= ddo:
            f = lai_min
            self._growth_stage = 0.
            self._senesc_stage = 0.
        elif self.DDsum > ddo:
            self._growth_stage += 1.0 / ddur
            f = np. minimum(1.0, lai_min + (1.0 - lai_min) * self._growth_stage)

        # senescence phase
        if doy > sso:
            self._growth_stage = 0.
            self._senesc_stage += 1.0 / sdur
            f = 1.0 - (1.0 - lai_min) * np.minimum(1.0, self._senesc_stage)

        # update self.LAIdecid and total LAI
        self._LAIdecid = self._LAIdecid_max * f
        self.LAI = self._LAIconif + self._LAIdecid
        return f

    def dry_canopy_et(self, D, Qp, AE, Ta, Ra=25.0, Ras=250.0, CO2=380.0, f=1.0, Rew=1.0, beta=1.0, fPheno=1.0):
        """
        Computes ET from 2-layer canopy in absense of intercepted preciptiation, i.e. in dry-canopy conditions
        IN:
           self - object
           D - vpd in kPa
           Qp - PAR in Wm-2
           AE - available energy in Wm-2
           Ta - air temperature degC
           Ra - aerodynamic resistance (s/m)
           Ras - soil aerodynamic resistance (s/m)
           CO2 - atm. CO2 mixing ratio (ppm)
           f - franction of local Rnet available for evaporation at ground [-]
           Rew - relative extractable water [-]
           beta - term for soil evaporation resistance (Wliq/FC) [-]
           fPheno - phenology modifier [-]
        Args:
           Tr - transpiration rate (mm s-1)
           Efloor - forest floor evaporation rate (mm s-1)
           Gc - canopy conductance (integrated stomatal conductance)  (m s-1)
        SOURCES:
        Launiainen et al. (2016). Do the energy fluxes and surface conductance
        of boreal coniferous forests in Europe scale with leaf area?
        Global Change Biol.
        Modified from: Leuning et al. 2008. A Simple surface conductance model
        to estimate regional evaporation using MODIS leaf area index and the
        Penman-Montheith equation. Water. Resources. Res., 44, W10419
        Original idea Kelliher et al. (1995). Maximum conductances for
        evaporation from global vegetation types. Agric. For. Met 85, 135-147

        Samuli Launiainen, Luke
        Last edit: 12 / 2017
        """

        # --- gsref as LAI -weighted average of conifers and decid.
        gsref = 1./self.LAI * (self._LAIconif * self.physpara['gsref_conif']
                + self._LAIdecid *self.physpara['gsref_decid'])

        kp = self.physpara['kp']  # (-) attenuation coefficient for PAR
        q50 = self.physpara['q50']  # Wm-2, half-sat. of leaf light response
        rw = self.physpara['rw']  # rew parameter
        rwmin = self.physpara['rwmin']  # rew parameter

        tau = np.exp(-kp * self.LAI)  # fraction of Qp at ground relative to canopy top

        """--- canopy conductance Gc (integrated stomatal conductance)----- """

        # fQ: Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
        fQ = 1./ kp * np.log((kp*Qp + q50) / (kp*Qp*np.exp(-kp * self.LAI) + q50 + eps) )

        # the next formulation is from Leuning et al., 2008 WRR for daily Gc; they refer to 
        # Kelliher et al. 1995 AFM but the resulting equation is not exact integral of K95.        
        # fQ = 1./ kp * np.log((Qp + q50) / (Qp*np.exp(-kp*self.LAI) + q50))

        # VPD -response
        fD = 1.0 / (np.sqrt(D) + eps)

        # soil moisture response: Lagergren & Lindroth, xxxx"""
        fRew = np.minimum(1.0, np.maximum(Rew / rw, rwmin))
        # fRew = 1.0

        # CO2 -response of canopy conductance, derived from APES-simulations
        # (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
        fCO2 = 1.0 - 0.387 * np.log(CO2 / 380.0)

        Gc = gsref * fQ * fD * fRew * fCO2 * fPheno
        Gc[np.isnan(Gc)] = eps

        """ --- transpiration rate --- """

        Tr = cu.penman_monteith((1.-tau)*AE, 1e3*D, Ta, Gc, 1./Ra, units='mm')
        Tr[Tr < 0] = 0.0

        """--- forest floor evaporation rate--- """
        # soil conductance is function of relative water availability
        gcs = 1. / self.soilrp * beta**2.0

        # beta = Wliq / FC; Best et al., 2011 Geosci. Model. Dev. JULES

        Efloor = cu.penman_monteith(tau * f * AE, 1e3*D, Ta, gcs, 1./Ras, units='mm')
        Efloor[f==0] = 0.0  # snow on ground restricts Efloor
        
        return Tr, Efloor, Gc  # , fQ, fD, fRew


    def canopy_water_snow(self, dt, T, Prec, AE, D, Ra=25.0, U=2.0):
        """
        Calculates canopy water interception and SWE during timestep dt
        Args: 
            self - object
            dt - timestep [s]
            T - air temperature (degC)
            Prec - precipitation rate during (mm d-1)
            AE - available energy (~net radiation) (Wm-2)
            D - vapor pressure deficit (kPa)
            Ra - canopy aerodynamic resistance (s m-1)
        Returns:
            self - updated state W, Wf, SWE, SWEi, SWEl
            Infil - potential infiltration to soil profile (mm)
            Evap - evaporation / sublimation from canopy store (mm)
            MBE - mass balance error (mm)
        Samuli Launiainen & Ari Laurén 2014 - 2017
        Last edit 12 / 2017
        """

        # quality of precipitation
        Tmin = 0.0  # 'C, below all is snow
        Tmax = 1.0  # 'C, above all is water
        Tmelt = 0.0  # 'C, T when melting starts

        # storage capacities mm
        Wmax = self.wmax * self.LAI
        Wmaxsnow = self.wmaxsnow * self.LAI

        # melting/freezing coefficients mm/s
        Kmelt = self.Kmelt - 1.64 * self.cf / (24 * 3600)  # Kuusisto E, 'Lumi Suomessa' --- oli /dt
        Kfreeze = self.Kfreeze

        kp = self.physpara['kp']
        tau = np.exp(-kp*self.LAI)  # fraction of Rn at ground

        # inputs to arrays, needed for indexing later in the code
        gridshape = np.shape(self.LAI)  # rows, cols
    
        if np.shape(T) != gridshape:
            T = np.ones(gridshape) * T
            Prec = np.ones(gridshape) * Prec
            AE = np.ones(gridshape) * AE
            D = np.ones(gridshape) * D
            Ra = np.ones(gridshape) * Ra

        Prec = Prec * dt  # mm
        # latent heat of vaporization (Lv) and sublimation (Ls) J kg-1
        Lv = 1e3 * (3147.5 - 2.37 * (T + 273.15))
        Ls = Lv + 3.3e5

        # compute 'potential' evaporation / sublimation rates for each grid cell
        erate = np.zeros(gridshape)
        ixs = np.where((Prec == 0) & (T <= Tmin))
        ixr = np.where((Prec == 0) & (T > Tmin))
        Ga = 1. / Ra  # aerodynamic conductance

        # resistance for snow sublimation adopted from:
        # Pomeroy et al. 1998 Hydrol proc; Essery et al. 2003 J. Climate;
        # Best et al. 2011 Geosci. Mod. Dev.
        # ri = (2/3*rhoi*r**2/Dw) / (Ce*Sh*W) == 7.68 / (Ce*Sh*W

        Ce = 0.01*((self.W + eps) / Wmaxsnow)**(-0.4)  # exposure coeff (-)
        Sh = (1.79 + 3.0*U**0.5)  # Sherwood numbner (-)
        gi = Sh*self.W*Ce / 7.68 + eps # m s-1
        # print ixs
        # print('ixs', np.shape(ixs), 'gi', np.shape(gi[ixs]), 'ga', np.shape(Ga[ixs]),
        #      'T', np.shape(T[ixs]), 'AE', np.shape(AE[ixs]), 'tau', np.shape(tau[ixs]))
        erate[ixs] = dt / Ls[ixs] * cu.penman_monteith((1.0 - tau[ixs])*AE[ixs], 1e3*D[ixs], T[ixs], Ga[ixs], gi[ixs], units='W')
#        print('gi', gi, 'Ce', Ce, 'Sh', Sh)

        # evaporation of intercepted water, mm
        gs = 1e6
        erate[ixr] = dt / Lv[ixr] * cu.penman_monteith((1.0 - tau[ixr])*AE[ixr], 1e3*D[ixr], T[ixr], Ga[ixr], gs, units='W')
        
        # print('erate', erate)

        # ---state of precipitation [as water (fW) or as snow(fS)]
        fW = np.zeros(gridshape)
        fS = np.zeros(gridshape)

        fW[T >= Tmax] = 1.0
        fS[T <= Tmin] = 1.0

        ix = np.where((T > Tmin) & (T < Tmax))
        fW[ix] = (T[ix] - Tmin) / (Tmax - Tmin)
        fS[ix] = 1.0 - fW[ix]
        del ix

        # --- Local fluxes (mm)
        Unload = np.zeros(gridshape)  # snow unloading
        Interc = np.zeros(gridshape)  # interception
        Melt = np.zeros(gridshape)   # melting
        Freeze = np.zeros(gridshape)  # freezing
        Evap = np.zeros(gridshape)

        """ --- initial conditions for calculating mass balance error --"""
        Wo = self.W  # canopy storage
        SWEo = self.SWE  # Snow water equivalent mm

        """ --------- Canopy water storage change -----"""
        # snow unloading from canopy, ensures also that seasonal LAI development does
        # not mess up computations
        ix = (T >= Tmax)
        Unload[ix] = np.maximum(self.W[ix] - Wmax[ix], 0.0)
        self.W = self.W - Unload
        del ix
        dW = self.W - Wo

        # Interception of rain or snow: asymptotic approach of saturation.
        # Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
        # Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
        ix = (T < Tmin)
        Interc[ix] = (Wmaxsnow[ix] - self.W[ix]) \
                    * (1.0 - np.exp(-(self.cf[ix] / Wmaxsnow[ix]) * Prec[ix]))
        del ix
        
        # above Tmin, interception capacity equals that of liquid precip
        ix = (T >= Tmin)
        Interc[ix] = np.maximum(0.0, (Wmax[ix] - self.W[ix]))\
                    * (1.0 - np.exp(-(self.cf[ix] / Wmax[ix]) * Prec[ix]))
        del ix
        self.W = self.W + Interc  # new canopy storage, mm

        Trfall = Prec + Unload - Interc  # Throughfall to field layer or snowpack

        # evaporate from canopy and update storage
        Evap = np.minimum(erate, self.W)  # mm
        self.W = self.W - Evap

        """ Snowpack (in case no snow, all Trfall routed to floor) """
        ix = np.where(T >= Tmelt)
        Melt[ix] = np.minimum(self.SWEi[ix], Kmelt[ix] * dt * (T[ix] - Tmelt))  # mm
        del ix
        ix = np.where(T < Tmelt)
        Freeze[ix] = np.minimum(self.SWEl[ix], Kfreeze * dt * (Tmelt - T[ix]))  # mm
        del ix

        # amount of water as ice and liquid in snowpack
        Sice = np.maximum(0.0, self.SWEi + fS * Trfall + Freeze - Melt)
        Sliq = np.maximum(0.0, self.SWEl + fW * Trfall - Freeze + Melt)

        PotInf = np.maximum(0.0, Sliq - Sice * self.R)  # mm
        Sliq = np.maximum(0.0, Sliq - PotInf)  # mm, liquid water in snow

        # update Snowpack state variables
        self.SWEl = Sliq
        self.SWEi = Sice
        self.SWE = self.SWEl + self.SWEi
        
        # mass-balance error mm
        MBE = (self.W + self.SWE) - (Wo + SWEo) - (Prec - Evap - PotInf)

        MBEsnow = self.SWE - SWEo - (Trfall - PotInf)
#        if np.nanmax(MBEsnow) > 0.1:
#        print('MBEsnow', MBEsnow)
#        print('Freeze', Freeze)
#        print('melt', Melt)
#        print('Potinf', PotInf)
#        print('Trfall', Trfall)
#        print('SWE', self.SWE)
#        print('SWo', SWEo)
#        print('dSWE', self.SWE-SWEo)
#        print('T', T)
#        print('fS', fS)
#        print('fW', fW)
#        print('Kmelt', Kmelt)
#        MBEcan = (self.W - Wo) - (Interc - Evap - Unload)
#        if np.nanmax(MBEcan) > 0.1:
#            print('Mbecan', MBEcan)
#            print('Unload', Unload)
#            print('dW', dW)

        return PotInf, Trfall, Evap, Interc, MBE, erate, Unload, fS + fW

